#!usr/bin/perl -w
############################################################################
# Module:      Nomenclature.pm   
# Description: Nomenclature functions in the DIGS tool
# History:     December  2017: Created by Robert Gifford 
############################################################################
package Nomenclature;

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
# TOP LEVEL FUNCTION
############################################################################

#***************************************************************************
# Subroutine:  create_standard_locus_ids
# Description: apply standard ids to a set of loci from one or more genomes
#***************************************************************************
sub create_standard_locus_ids {

	my ($self, $infile) = @_;

 	# Show title
	my $digs_obj = $self->{digs_obj};
	$digs_obj->show_title();  

	# Set up for this process
	$self->initialise_nomenclature_process($infile);

	# Get the locus data and created clustered annotations
	$self->cluster_annotations();

	# Apply standard names to locus clusters
	$self->apply_standard_names_to_clusters();

}

############################################################################
# INITIALISATION: top level handlers
############################################################################

#***************************************************************************
# Subroutine:  initialise_nomenclature_process
# Description: set up the program to execute the ID allocation process 
#***************************************************************************
sub initialise_nomenclature_process {

	my ($self, $infile) = @_;

	# Set the 'defragment_mode' (determines rules for clustering loci)
	my $digs_obj = $self->{digs_obj};
	$digs_obj->{defragment_mode} = 'consolidate';
	
	# Initialise database
	$self->initialise_nomenclature_db($infile);

	# Do console dialogue to get input parameters
	$self->do_console_setup_dialogue();

}

#***************************************************************************
# Subroutine:  initialise_nomenclature_process
# Description: set up the program to execute the ID allocation process 
#***************************************************************************
sub do_console_setup_dialogue {

	my ($self) = @_;

	print "\n\n\t  #≈#~# SETUP: overview\n";
	print "\n\t  # Please define the following:";
	print "\n\t  # ";
	print "\n\t  #   1.   Locus class (e.g. ERV)";
	print "\n\t  #   2.   A source of tracks (required)";
	print "\n\t  #   3.   A taxonomy translation table for the locus class (optional)";
	print "\n\t  #   4.   An alternative name translation table for gene names (optional)";
	print "\n\t  #   5.   An organism ID, or organism name=>ID translation table";

	# Get the locus class
	print "\n\n\t  #≈#~# SET LOCUS CLASS\n";
	my $question2 = "\n\t  What locus class name to use?";
	my $locus_class = $console->ask_question($question2);	
	#$self->{locus_class} = $locus_class;
	$self->{locus_class} = 'ERV';	

	# Nomenclature tracks
	$self->load_nomenclature_tracks();
	
	# Set the rules for using translation tables (options)
	$self->load_taxonomy_tables();

	# Set the rules for handling gene names (options)
	
	# Get the organism ID
	my $question1 = "\n\t  What is the organism code to use?";
	#my $organism_code = $console->ask_yes_no_question($question1);	
	#$self->{organism_code} = $organism_code;
	$self->{organism_code} = 'HoSa';
}

############################################################################
# INITIALISATION: core database tables for ID allocation
############################################################################

#***************************************************************************
# Subroutine:  initialise_nomenclature_db
# Description: load core nomenclature tables, create if they don't exist
#***************************************************************************
sub initialise_nomenclature_db {

	my ($self, $infile) = @_;
	
	# Parse control file and connect to DB
	unless($infile) {
		die "\n\t This option requires an infile\n\n";
	}

	# Parse control file and connect to DB
	$self->parse_ctl_file_and_connect_to_db($infile);

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
		$db_ref->create_nomenclature_tracks_table($dbh);
	}
	my $nom_chains_exists = $db_ref->does_table_exist('nomenclature_chains');
	unless ($nom_chains_exists) {
		$db_ref->create_nomenclature_chains_table($dbh);
	}

	# Load nomenclature tables
	$db_ref->load_nomenclature_tracks_table($dbh);
	$db_ref->load_nomenclature_chains_table($dbh);
	$db_ref->load_nomenclature_table($dbh);	
}

#***************************************************************************
# Subroutine:  parse_ctl_file_and_connect_to_db
# Description: connect to a DIGS screening DB by parsing a DIGS control file
#***************************************************************************
sub parse_ctl_file_and_connect_to_db {

	my ($self, $infile) = @_;

	my $digs_obj = $self->{digs_obj};
	
	# Try opening control file
	my @ctl_file;
	my $valid = $fileio->read_file($infile, \@ctl_file);
	unless ($valid) {  # Exit if we can't open the file
		die "\n\t ### Couldn't open control file '$infile'\n\n\n ";
	}
	
	# If control file looks OK, store the path and parse the file
	$self->{ctl_file} = $infile;
	my $loader_obj = ScreenBuilder->new($digs_obj);
	$loader_obj->parse_control_file($infile, $digs_obj);

	# Store the ScreenBuilder object (used later)
	$self->{loader_obj} = $loader_obj;

	# Load/create the screening database
	my $db_name = $loader_obj->{db_name};
	unless ($db_name) { die "\n\t Error: no DB name defined \n\n\n"; }
	$digs_obj->initialise_screening_db($db_name);
}

#***************************************************************************
# Subroutine:  do_flush_tables_dialogue
# Description: ask user whether or not to flush the core nomenclature tables
#***************************************************************************
sub do_flush_tables_dialogue {

	my ($self) = @_;

	my $digs_obj = $self->{digs_obj};
	my $db_ref = $digs_obj->{db};
	my $dbh = $db_ref->{dbh};

	# Check whether to flush the table
	my $tracks_table = $db_ref->{nomenclature_tracks_table};
	my $chains_table = $db_ref->{nomenclature_chains_table};
	my $nom_table    = $db_ref->{nomenclature_table};
	unless ($nom_table and $tracks_table and $chains_table) { die; }
	my $question1 = "\n\n\t  Flush the tables before uploading tracks?";
	my $flush = $console->ask_yes_no_question($question1);
	if ($flush eq 'y') { 
		$tracks_table->flush();
		$chains_table->flush();
		$nom_table->flush();
	}	
}

############################################################################
# INITIALISATION: core database tables for ID allocation
############################################################################

#***************************************************************************
# Subroutine:  load_nomenclature_tracks
# Description: load input tracks into table, in a DIGS locus format
#***************************************************************************
sub load_nomenclature_tracks {

	my ($self) = @_;

	# Option to load a taxon translation table
	# Reason: tables allow IDs to be based on alternative taxonomic levels
	print "\n\n\t  #≈#~# SET SOURCE OF ANNOTATION TRACKS\n";
	my $question1 = "\n\t  Load tracks from file, or use a screening DB table";
	my @choices = qw [ f t ];
	my $load_from_file = $console->ask_simple_choice_question($question1, \@choices);
	if ($load_from_file eq 'f') {
		$self->load_nomenclature_tracks_from_file();
	}
	elsif ($load_from_file eq 't') {

	}
}

#***************************************************************************
# Subroutine:  load_nomenclature_tracks_from_file
# Description: load input tracks into table, in a DIGS locus format
#***************************************************************************
sub load_nomenclature_tracks_from_file {

	my ($self) = @_;
	
	# Load nomenclature table
	my $digs_obj = $self->{digs_obj};
	my $db_ref = $digs_obj->{db};

	$self->do_flush_tables_dialogue();

	my $nom_table = $db_ref->{nomenclature_tracks_table};
	unless ($nom_table) { die; }

	# Read tracks from file path
	print "\n\n\t #### WARNING: This function expects a tab-delimited data table with column headers!";
	my $question1 = "\n\n\t Please enter the path to the file with the table data and column headings\n\n\t";
	my $infile = $console->ask_question($question1);
	unless ($infile) { die; }
	my @new_track;
	$fileio->read_file($infile, \@new_track);

	# Load tracks into table
	my $line_number = 0;
	foreach my $line (@new_track) {

		$line_number++;
		if     ($line =~ /^\s*$/)  { next; } # discard blank line
		elsif  ($line =~ /^\s*#/)  { next; } # discard comment line 
		unless ($line =~ /\t/)     { print "\n\t Incorrect formatting at line '$line_number'"; die; }
	
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
# Subroutine:  load_taxonomy_tables
# Description: load a table that captures taxonomic definitions
#***************************************************************************
sub load_taxonomy_tables {

	my ($self) = @_;

	# Option to load a taxon translation table
	# Reason: tables allow IDs to be based on alternative taxonomic levels
	print "\n\n\t  #≈#~# LOAD TAXONOMY TRANSLATION TABLES\n";
	my $question1 = "\n\t  Load a taxonomy table from file?";
	my $load_from_file = $console->ask_yes_no_question($question1);
	if ($load_from_file eq 'y') {
		my %translations;
		$self->load_translations_from_file(\%translations);
		$self->{translations} = \%translations;
	}
	else {
		my $question2 = "\n\t  Use a table in the screening database?";
		my $load_from_db_table = $console->ask_yes_no_question($question2);
		if ($load_from_db_table eq 'y') {
			my %translations;
			$self->load_translations_from_table(\%translations);
			$self->{translations} = \%translations;
		}
	}
}

#***************************************************************************
# Subroutine:  load_translations_from_file
# Description: load translation tables
#***************************************************************************
sub load_translations_from_file {

	my ($self, $translations_ref) = @_;

	# Read translation from file
	#my $translations_path = $self->{translation_path};
	my $translations_path =  '../local/human/translations.txt';
	die;
	
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
			unless ($level) { 
				$level = 'wtf';
				#die; 
			}		
			$taxonomy{$level} = $value;			
		}
		my $id = shift @line;
		$translations_ref->{$id} = \%taxonomy;		
	}
}

#***************************************************************************
# Subroutine:  load_gene_name_translation_tables
# Description: load a table that captures alternative names of homologous genes
#***************************************************************************
sub load_gene_name_translation_tables {

	my ($self) = @_;

	# Option to load a gene look-up table (resolve to equivalents)
	# TODO
	print "\n\n\t #### GENE NAME TRANSLATION TABLES";
	my $question1 = "\n\t  Load a gene name translation table from file?";
	my $load_from_file = $console->ask_yes_no_question($question1);
	if ($load_from_file eq 'y') {	
		my %gene_names;
		die;
		#$self->load_gene_name_translations_from_file(\%translations);
		#$self->{gene_names} = \%gene_names;
	}
	else {
		my $question2 = "\n\t  Use a gene name translation table in the screening database?";
		my $load_from_db_table = $console->ask_yes_no_question($question2);
		if ($load_from_db_table eq 'y') {	
			my %gene_names;
			#$self->load_translations_from_table(\%gene_names);
			#$self->{gene_names} = \%gene_names;
		}
	}
}

############################################################################
# ANNOTATION CLUSTERING
############################################################################

#***************************************************************************
# Subroutine:  cluster_annotations
# Description: cluster nomenclature tracks
#***************************************************************************
sub cluster_annotations {

	my ($self) = @_;

	# Get sorted tracks from nomenclature table
	$self->get_sorted_nomenclature_tracks();
		
	# Compose clusters of related sequences
	my %settings;
	$settings{range} = '0';
	$settings{start} = 'extract_start';
	$settings{end}   = 'extract_end';
	my %clusters;
	my $digs_obj = $self->{digs_obj};
	my $sorted_ref = $self->{sorted_loci};
	$digs_obj->compose_clusters(\%clusters, $sorted_ref, \%settings);
	#$devtools->print_hash(\%clusters); die;

	# Cluster IDs
	my @cluster_ids  = keys %clusters;
	my $num_clusters = scalar @cluster_ids;
	print "\n\t # $num_clusters locus groups in total";
	$self->{nomenclature_clusters} = \%clusters;	
	
}

#***************************************************************************
# Subroutine:  get_sorted_nomenclature_tracks
# Description: get nomenclature set rows, sorted by scaffold, in order of location
#***************************************************************************
sub get_sorted_nomenclature_tracks {

	my ($self) = @_;;

	# Set SQL 'where' clause to sort the rows
	my $where  = " ORDER BY scaffold, extract_start ";

	# Get database tables
	my $digs_obj = $self->{digs_obj};
	my $db = $digs_obj->{db};
	my $nomenclature_table = $db->{nomenclature_tracks_table};
	unless ($nomenclature_table) {  
		$devtools->print_hash($db); die; 
	}
		
	# Set the fields to get values for
	my @sorted;
	my @fields = qw [ record_id track_name
	                  assigned_name assigned_gene
	                  scaffold orientation namespace_id
                      extract_start extract_end sequence_length ];
	$nomenclature_table->select_rows(\@fields, \@sorted, $where);

	my $total_loci = scalar @sorted;
	print "\n\n\t # $total_loci rows in the nomenclature table";
	#$devtools->print_array(\@sorted); die;
	$self->{sorted_loci} = \@sorted;	
}

############################################################################
# STANDARD ID ALLOCATION
############################################################################

#***************************************************************************
# Subroutine:  apply_standard_names_to_clusters
# Description: apply standard names to clusters  
#***************************************************************************
sub apply_standard_names_to_clusters {

	my ($self) = @_;

	my $digs_obj      = $self->{digs_obj};
	my $db_ref        = $digs_obj->{db};
	my $chains_table  = $db_ref->{nomenclature_chains_table};
	my $nom_table     = $db_ref->{nomenclature_table};
	unless ($nom_table and $chains_table) { die; }

	my $organism_code = $self->{organism_code};
	my $locus_class   = $self->{locus_class};
	unless ($organism_code and $locus_class) { die; } # Sanity checking

	# Configure namespace and preprocess
	$self->configure_namespace();

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
# Subroutine:  configure_namespace
# Description: configure the namespace for named loci  (also preprocess for locus assign)
#***************************************************************************
sub configure_namespace {

	my ($self) = @_;

	# Set the organism and version fields (a convenience), & configure namespace
	my $clusters_ref = $self->{nomenclature_clusters};	
	my @cluster_ids  = keys %$clusters_ref;

	my $organism       = $self->{nomenclature_organism};
	my $target_version = $self->{nomenclature_version};
	my %namespace;
	foreach my $cluster_id (@cluster_ids) {
		
		my $cluster_ref = $clusters_ref->{$cluster_id};
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
# Subroutine:  translate_taxaname
# Description: translate taxa name
#***************************************************************************
sub translate_taxaname {

	my ($self, $original_name, $taxonomic_level) = @_;

	# Do translation
	my $translations_ref = $self->{translations};
	my $taxonomy_ref = $translations_ref->{$original_name};	
	my $translated_name = $taxonomy_ref->{$taxonomic_level};	
	return $translated_name;
}

############################################################################
# UTILITY / DEV
############################################################################

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