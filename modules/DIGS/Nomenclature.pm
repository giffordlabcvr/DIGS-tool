#!usr/bin/perl -w
############################################################################
# Module:      Nomenclature.pm   
# Description: Nomenclature functions in the DIGS tool
# History:     December  2017: Created by Robert Gifford 
############################################################################
package Test;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::FileIO;
use Base::Console;
use Base::DevTools;

# Program components
use DIGS::ScreenBuilder; # Functions to set up screen

############################################################################
# Globals
############################################################################

# Base objects
my $fileio    = FileIO->new();
my $console   = Console->new();
my $devtools  = DevTools->new();

# DIGS test database connection globals
my $server   = 'localhost';
my $user     = 'root';
my $password = 'blenat2';
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create new Test 'object'
#***************************************************************************
sub new {

	my ($invocant, $digs_obj) = @_;
	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {

		# Global settings
		process_id             => $digs_obj->{process_id},
		program_version        => $digs_obj->{program_version},
		
		# DIGS tool object
		digs_obj => $digs_obj,

	};
	
	bless ($self, $class);
	return $self;
}


############################################################################
# INTERNAL FUNCTIONS: nomenclature
############################################################################

#***************************************************************************
# Subroutine:  create_standard_locus_ids
# Description: apply standard ids to a set of loci from one or more genomes
#***************************************************************************
sub create_standard_locus_ids {

	my ($self, $infile) = @_;

	# Create nomenclature tables if they don't exist already
	my $digs_obj = $self->{digs_obj};
	my $db_ref = $digs_obj->{db};
	my $dbh = $db_ref->{dbh};
	my $nomenclature_exists = $db_ref->does_table_exist('nomenclature');
	unless ($nomenclature_exists) {
		$db_ref->create_nomenclature_table($dbh);
	}
	my $nom_tracks_exists = $db_ref->does_table_exist('nomenclature_tracks');
	unless ($nom_tracks_exists) {
		$db_ref->create_nomenclature_track_table($dbh);
	}
	my $nom_chains_exists = $db_ref->does_table_exist('nomenclature_chains');
	unless ($nom_chains_exists) {
		$db_ref->create_nomenclature_chains_table($dbh);
	}

	# Load nomenclature tables
	$db_ref->load_nomenclature_tracks_table($dbh);
	$db_ref->load_nomenclature_chains_table($dbh);
	$db_ref->load_nomenclature_table($dbh);
	my $tracks_table = $db_ref->{nomenclature_tracks_table};
	my $chains_table = $db_ref->{nomenclature_chains_table};
	my $nom_table    = $db_ref->{nomenclature_table};
	unless ($nom_table and $tracks_table and $chains_table) { die; }

	# Check whether to flush the table
	my $question = "\n\n\t  Flush the tables before uploading tracks?";
	my $flush = $console->ask_yes_no_question($question);
	if ($flush eq 'y') { 
		$tracks_table->flush();
		$chains_table->flush();
		$nom_table->flush();
	}

	# Load tracks into table in a DIGS locus format
	$self->load_nomenclature_tracks();

	# Cluster tracks
	$self->create_nomenclature_track();

	# Apply standard names to locus clusters
	my $organism_code = $self->{organism_code};
	my $locus_class   = $self->{locus_class};
	unless ($organism_code and $locus_class) { die; } # Sanity checking
	$self->apply_standard_names_to_clusters($locus_class, $organism_code);

}

#***************************************************************************
# Subroutine:  create_nomenclature_track
# Description: cluster nomenclature tracks
#***************************************************************************
sub create_nomenclature_track {

	my ($self) = @_;
	
	# Get sorted tracks from nomenclature table
	my @sorted;
	$self->get_sorted_nomenclature_tracks(\@sorted);
	#$devtools->print_array(\@sorted); die;
	my $total_loci = scalar @sorted;
	print "\n\n\t # $total_loci rows in the nomenclature table";

	# Compose clusters of related sequences
	my %settings;
	$settings{range} = '0';
	$settings{start} = 'extract_start';
	$settings{end}   = 'extract_end';
	my %clusters;
	my $digs_obj = $self->{digs_obj};
	$digs_obj->compose_clusters(\%clusters, \@sorted, \%settings);
	$devtools->print_hash(\%clusters); die;

	# Cluster IDs
	my @cluster_ids  = keys %clusters;
	my $num_clusters = scalar @cluster_ids;
	print "\n\t # $num_clusters locus groups in total";
	$self->{nomenclature_clusters} = \%clusters;	

	# Set the organism and version fields (a convenience), & configure namespace
	my %namespace;
	my $organism       = $self->{nomenclature_organism};
	my $target_version = $self->{nomenclature_version};
	foreach my $cluster_id (@cluster_ids) {
		
		my $cluster_ref = $clusters{$cluster_id};
		#$devtools->print_array($cluster_ref); exit;
		foreach my $locus_ref (@$cluster_ref) {		
			
			# Set organism and target version
			$locus_ref->{organism}       = $organism;
			$locus_ref->{target_version} = $target_version;			

			# Check namespace
			my $assigned_name = $locus_ref->{assigned_name};
			my $namespace_id  = $locus_ref->{namespace_id};
			if ($namespace_id ne 'NULL') {
				
				if ($namespace{$assigned_name}) {			
					my $taxon_namespace_ref = $namespace{$assigned_name};
					$taxon_namespace_ref->{$namespace_id} = 1;
				}
				else{
					my %taxon_namespace;
					$taxon_namespace{$namespace_id} = 1;
					$namespace{$assigned_name} = \%taxon_namespace;
				}				
			}
		}
	}
	$self->{namespace} = \%namespace;
}

#***************************************************************************
# Subroutine:  apply_standard_names_to_clusters
# Description: apply standard names to clusters  
#***************************************************************************
sub apply_standard_names_to_clusters {

	my ($self, $locus_class, $organism_code) = @_;

	my $db_ref        = $self->{db};
	my $chains_table  = $db_ref->{nomenclature_chains_table};
	my $nom_table     = $db_ref->{nomenclature_table};
	unless ($nom_table and $chains_table) { die; }

	# Load translations
	my %translations;
	$self->load_translations(\%translations);

	# Iterate through the clusters
	my %counter;
	my @nomenclature;
	my $clusters_ref = $self->{nomenclature_clusters};	
	my @cluster_ids  = keys %$clusters_ref;
	#$self->show_clusters($clusters_ref);
	foreach my $cluster_id (@cluster_ids) {

		my $mixed;
		my $skip;
		my $lowest;
		my $highest;
		my $taxname_final;
		my $scaffold_final;
		my $orientation_final;
		my $namespace_id_final = undef;
		
		# Get the array of loci
		my $cluster_ref = $clusters_ref->{$cluster_id};
		my %last_locus;
		my %composite;
		foreach my $locus_ref (@$cluster_ref) {
		
			#$devtools->print_hash($locus_ref);
			my $start        = $locus_ref->{extract_start};
			my $end          = $locus_ref->{extract_end};
			my $track        = $locus_ref->{track_name};
			my $taxname      = $locus_ref->{assigned_name};
			my $namespace_id = $locus_ref->{namespace_id};
			my $orientation  = $locus_ref->{orientation};
			my $scaffold     = $locus_ref->{scaffold};
			if ($namespace_id ne 'NULL') {	
				$namespace_id_final = $namespace_id;
			}

			# Set coordinates
			my $last_start   = $last_locus{extract_start};
			my $last_end     = $last_locus{extract_end};
			
			if ($last_start) {
				if ($start < $last_start) { $lowest = $start; }
				if ($end > $last_end)     { $highest = $end;  }
			}
			else {
				$lowest = $start;
				$highest = $end;
			}
			
			$taxname_final = $taxname;
			$scaffold_final = $scaffold;
			$orientation_final = $orientation;
		}
		
		# Create numeric ID
		$composite{track_name}    = 'MASTER';
		$composite{assigned_name} = $taxname_final;
		$composite{scaffold}      = $scaffold_final;
		$composite{extract_start} = $lowest;
		$composite{extract_end}   = $highest;
		$composite{orientation}   = $orientation_final;
		$composite{locus_class}   = $locus_class;
		$composite{organism_code} = $organism_code;
		$composite{namespace_id}  = $namespace_id_final;
		my $numeric_id = $self->create_numeric_id(\%composite, \%counter);

		# Create the locus ID		
		my @id;
		push (@id, $locus_class);
		push (@id, $taxname_final);
		push (@id, $numeric_id);
		push (@id, $organism_code);	 
		my $id = join('.', @id);
		print "\n\t # ID: $id";

		# Update the nomenclature table
		$composite{full_id} = $id;
		$composite{namespace_id} = $numeric_id;
		my $nom_id = $nom_table->insert_row(\%composite);
	
		# Update chains table
		foreach my $locus_ref (@$cluster_ref) {			
			my $locus_id = $locus_ref->{record_id};
			my %data;
			$data{track_id}             = $locus_id;
			$data{nomenclature_locus_id} = $nom_id;
			$chains_table->insert_row(\%data);
		}
	}	
}

#***************************************************************************
# Subroutine:  create_numeric_id
# Description: create a unique ID for a locus (tracking an ID namespace)
#***************************************************************************
sub create_numeric_id {

	my ($self, $locus_ref, $counter_ref) = @_;

	# Get data structures and values
	my $taxname       = $locus_ref->{assigned_name};
	my $namespace_id  = $locus_ref->{namespace_id};
	my $namespace_ref = $self->{namespace};

	# Create numeric ID
	my $numeric_id;			
	if ($namespace_id) {	
		$numeric_id = $namespace_id;
	}
	else {			

		# Get namespace for this taxon
		my $taxon_namespace = $namespace_ref->{$taxname};
		
		# Create new numeric ID that does not infringe the current namespace
		my $unique = undef;
		
		if ($counter_ref->{$taxname}) {
			$numeric_id = $counter_ref->{$taxname};	
		}
		else {
			$numeric_id = '0';
		}
						
		do { # Increment until we get a number that doesn't infringe the namespace
			$numeric_id++;
			unless ($taxon_namespace->{$numeric_id}) {
				$unique = 'true';
			}
		} until ($unique);
			
		$counter_ref->{$taxname} = $numeric_id;
	}	
	return $numeric_id;
}

#***************************************************************************
# Subroutine:  load_translations
# Description: load translation tables
#***************************************************************************
sub load_translations {

	my ($self, $translations_ref) = @_;

	# Read translation from file
	my $translations_path = $self->{translation_path};
	unless ($translations_path) { die; }
	my @file;
	$fileio->read_file($translations_path, \@file);
	my $header = shift @file;
	chomp $header;
	my @header = split("\t", $header); 
	my %levels;
	my $i = 0;
	foreach my $element (@header) {
		$i++;
		$levels{$i} = $element;
	}

	# Set up the translations
	foreach my $line (@file) {

		chomp $line;
		my @line  = split("\t", $line);
		my $j = 0;
		my %taxonomy;
		foreach my $value (@line) {
			$j++;
			my $level = $levels{$j};
			unless ($level) { die; }		
			$taxonomy{$level} = $value;			
		}
		my $id = shift @line;
		$translations_ref->{$id} = \%taxonomy;		
	}
}

#***************************************************************************
# Subroutine:  load_nomenclature_tracks
# Description: load input tracks into table, in a DIGS locus format
#***************************************************************************
sub load_nomenclature_tracks {

	my ($self) = @_;
	
	# Load nomenclature table
	my $db_ref = $self->{db};
	my $nom_table = $db_ref->{nomenclature_tracks_table};
	unless ($nom_table) { die; }

	# Read tracks from file path
	my @new_track;
	my $new_track_path = $self->{new_track_path};
	unless ($new_track_path) { die; }
	$fileio->read_file($new_track_path, \@new_track);
	
	# Load tracks into table
	foreach my $line (@new_track) {
	
		#print $line;
		chomp $line;
		my @line = split("\t", $line);
		my $track_name    = shift @line;
		my $assigned_name = shift @line;
		my $scaffold      = shift @line;
		my $extract_start = shift @line;
		my $extract_end   = shift @line;
		my $assigned_gene = shift @line;
		my $namespace_id  = shift @line;		
		my $span;
		my $orientation;

		#
		if ($extract_end > $extract_start) {
			$orientation = '+';
			$span = $extract_end - $extract_start;
		}
		else {
			$orientation = '-';
			$span = $extract_start - $extract_end;		
			my $start = $extract_start;
			$extract_start = $extract_end;
			$extract_end   = $start;
		}

		my %data;
		$data{track_name}      = $track_name;
		$data{assigned_name}   = $assigned_name;
		$data{scaffold}        = $scaffold;
		$data{extract_start}   = $extract_start;
		$data{extract_end}     = $extract_end;
		$data{sequence_length} = $span;
		$data{orientation}     = $orientation;
		$data{assigned_gene}   = $assigned_gene;
		unless ($namespace_id) { $namespace_id = 'NULL'; }
		$data{namespace_id}    = $namespace_id;		
		$nom_table->insert_row(\%data);	

	}
}

#***************************************************************************
# Subroutine:  get_sorted_nomenclature_tracks
# Description: get nomenclature set rows, sorted by scaffold, in order of location
#***************************************************************************
sub get_sorted_nomenclature_tracks {

	my ($self, $data_ref, $where) = @_;;

	# Set SQL 'where' clause to sort the rows
	my $sort  = " ORDER BY scaffold, extract_start ";
	if ($where) { $where .= $sort; }
	else        { $where  = $sort; }

	# Get database tables
	my $digs_obj = $self->{digs_obj};
	my $db = $digs_obj->{db};
	my $nomenclature_table = $db->{nomenclature_tracks_table};
	unless ($nomenclature_table) {  
		$devtools->print_hash($db); die; 
	}
		
	# Set the fields to get values for
	my @fields = qw [ record_id track_name
	                  assigned_name assigned_gene
	                  scaffold orientation namespace_id
                      extract_start extract_end sequence_length ];
	$nomenclature_table->select_rows(\@fields, $data_ref, $where);
}

#***************************************************************************
# Subroutine:  show_translations
# Description:  
#***************************************************************************
sub show_translations {

	my ($self, $translations_ref) = @_;

	# Validate
	#show_translations(\%taxonomy); #die;
		
	my @keys = keys %$translations_ref;
	foreach my $key (@keys) {
		my $data_ref = $translations_ref->{$key};
		my @keys2 = keys %$data_ref;
		foreach my $key2 (@keys2) {
			my $value = $data_ref->{$key2};
			print "\n\t $key	$key2	$value";
		}
	}
}

############################################################################
# EOF
############################################################################