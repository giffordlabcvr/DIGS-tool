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

	# Set the 'defragment_mode' (determines rules for clustering loci)
	$digs_obj->{defragment_mode} = 'consolidate';
	
	# Initialise database
	$self->initialise_nomenclature_db($infile);

	# Do console dialogue to get input parameters
	$self->do_console_setup_dialogue();

	# Get sorted tracks from nomenclature table
	$self->get_sorted_tracks();
		
	# Get the locus data and created clustered annotations
	$self->cluster_annotations();
	
	# Apply standard names to locus clusters
	$self->apply_standard_names_to_clusters();

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

	# Show the results of clustering
	my @cluster_ids  = keys %clusters;
	my $num_clusters = scalar @cluster_ids;
	print "\n\t # $num_clusters locus groups in total";
	$self->{nomenclature_clusters} = \%clusters;	
	
}

#***************************************************************************
# Subroutine:  get_sorted_tracks
# Description: get nomenclature set rows, sorted by scaffold, in order of location
#***************************************************************************
sub get_sorted_tracks {

	my ($self) = @_;;

	# Set SQL 'where' clause to sort the rows
	my $where  = " ORDER BY scaffold, start_position ";

	# Get database tables
	my $digs_obj = $self->{digs_obj};
	my $db_ref = $digs_obj->{db};	
	my $nomenclature_table = $db_ref->{nomenclature_tracks_table};
	unless ($nomenclature_table) {  
		$devtools->print_hash($db_ref); die; 
	}
		
	# Set the fields to get values for
	my @sorted;
	my @fields = qw [ record_id track_name
	                  organism_code locus_class
	                  taxon gene scaffold orientation
	                  start_position end_position
	                  namespace_id  ];
	$nomenclature_table->select_rows(\@fields, \@sorted, $where);

	my $total_annotations = scalar @sorted;
	print "\n\n\t # $total_annotations rows in the track table";
	#$devtools->print_array(\@sorted); die;

	# Format for DIGS clustering	
	my $start;
	my $last_start;
	my $scaffold;
	my $last_scaffold;
	
	my $count = 0;
	foreach my $row (@sorted) {
		$count++;
		$row->{assigned_name} = $row->{taxon};
		$row->{extract_start} = $row->{start_position};
		$row->{extract_end}   = $row->{end_position};
		$row->{assigned_gene} = $row->{gene};
	}
	
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
		#print "\n\t # ID: $id";

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
# INITIALISATION: top level handlers
############################################################################

#***************************************************************************
# Subroutine:  initialise_nomenclature_process
# Description: set up the program to execute the ID allocation process 
#***************************************************************************
sub do_console_setup_dialogue {

	my ($self) = @_;

	print "\n\n\t  #≈#~# SETUP: overview\n";
	#print "\n\t  # Please define the following:\n";
	#print "\n\t  #   1.   Locus class (e.g. ERV)";
	#print "\n\t  #   2.   A source of tracks (required)";
	#print "\n\t  #   3.   A taxonomy translation table for the locus class (optional)";
	#print "\n\t  #   4.   An alternative name translation table for gene names (optional)";
	#print "\n\t  #   5.   An organism ID, or a translation table linking organism name with an ID";

	# Get the locus class
	print "\n\n\t  #≈#~# 1: SET LOCUS CLASS\n";
	my $question2 = "\n\t  What locus class name to use?";
	#my $locus_class = $console->ask_question($question2);	
	my $locus_class = 'ERV';	
	print "\n\t      - locus class set to '$locus_class'";
	$self->{locus_class} = $locus_class;	

	# Set the source for nomenclature tracks
	my $track_table = $self->set_nomenclature_tracks();
	print "\n\t      - annotations sourced from table '$track_table'";
	
	# Set the rules for using translation tables (options)
	my $taxonomy_table = $self->load_taxonomy_table();
	print "\n\t      - taxonomy translations sourced from table '$taxonomy_table'";

	# Set the rules for handling gene names (options)
	my $gene_name_table = $self->load_gene_name_translation_table();
	print "\n\t      - gene name translations sourced from table '$gene_name_table'";
	
	# Get the organism ID
	print "\n\n\t  #≈#~# 5: SET ORGANISM CODE\n";
	my $question1 = "\n\t  What is the organism code to use?";
	#my $organism_code = $console->ask_yes_no_question($question1);	
	my $organism_code = 'HoSa';
	$self->{organism_code} = $organism_code;
	print "\n\t      - organism code set to '$organism_code'";

}

############################################################################
# INITIALISATION: core database tables for ID allocation
############################################################################

#***************************************************************************
# Subroutine:  set_nomenclature_tracks
# Description: set the source of annotations for ID allocation
#***************************************************************************
sub set_nomenclature_tracks {

	my ($self) = @_;

	# Get DIGS object and DB
	my $digs_obj  = $self->{digs_obj};
	my $db        = $digs_obj->{db};  
	my $dbh       = $db->{dbh};  

	# Option to load a taxon translation table
	# Reason: tables allow IDs to be based on alternative taxonomic levels
	print "\n\n\t  #≈#~# 2: SET SOURCE OF ANNOTATION TRACKS\n";
	my $table_name;
	my $question1 = "\n\t  Load tracks from file, or use a screening DB table";
	my @choices = qw [ l t ];
	my $load_from_file = $console->ask_simple_choice_question($question1, \@choices);

	#my $load_from_file = 't';
	if ($load_from_file eq 'l') {
		$table_name = $self->do_load_tracks_to_tables_dialogue();
	}
	elsif ($load_from_file eq 't') {
		$table_name = 'herv_input';
		$db->load_tracks_table($dbh, $table_name);
	}

	# Set table 
	$self->{input_track_table_name} = $table_name;
	return $table_name;

}

#***************************************************************************
# Subroutine:  do_load_tracks_to_tables_dialogue
# Description: 
#***************************************************************************
sub do_load_tracks_to_tables_dialogue {

	my ($self, $track_data_ref) = @_;

	# Get DIGS object and DB
	my $digs_obj  = $self->{digs_obj};
	my $db        = $digs_obj->{db};  
	my $dbh       = $db->{dbh};  

	# Get the file path
	print "\n\n\t #### WARNING: This function expects a tab-delimited data table with column headers!";
	my $question = "\n\n\t Please enter the path to the file with the table data and column headings\n\n\t";
	#my $infile = $self->ask_question($question);
	my $infile = '../local/human/hg19_combined_tracks.txt';
	unless ($infile) { die; }

	# Read in the data from a tab delimited file
	my @data;
	my %fields;
	my @fields;
	$console->do_read_tabdelim_dialogue($infile, \@data, \@fields, \%fields);

	# Get a reference to a table object for the ancillary table
	my $table_name = $db->do_track_table_dialogue(\@data, \@fields, \%fields);		
	return $table_name;

}

#***************************************************************************
# Subroutine:  choose_nomenclature_table
# Description: 
#***************************************************************************
sub choose_nomenclature_table {

	my ($self, $answer) = @_;

	# Get DIGS object and DB
	my $digs_obj = $self->{digs_obj};
	my $db = $digs_obj->{db};  


	my %fields;
	my @fields;

	my %extra_tables;
	my @extra_tables = qw [ digs_results ];

	# Create new table if option 1 selected
	my $table_name;
	if ($answer eq '1') {	
		my $table_name_question = "\n\t What is the name of the new table?";
		$table_name = $console->ask_question($table_name_question);
		$db->create_ancillary_table($table_name, \@fields, \%fields);	

	}
	# or choose one of the ancillary tables already in the DB
	else {

		# Get the ancillary tables in this DB
		$db->get_ancillary_table_names(\@extra_tables);
		
		my $table_num = 0;
		foreach my $extra_table_name (@extra_tables) {
		
			if ($extra_table_name eq 'digs_results'
			or  $extra_table_name eq 'nomenclature') {
				next;
			}
			else {
				$table_num++;
				$extra_tables{$table_num} = $extra_table_name;
				print "\n\t\t Table $table_num: '$extra_table_name'";
			}
		}
		my @table_choices = keys %extra_tables;

		my $question5 = "\n\n\t Apply to which of the above tables?";
		my $answer5   = $console->ask_simple_choice_question($question5, \@table_choices);
		$table_name = $extra_tables{$answer5};
		unless ($table_name) { die; }
	}
	return $table_name;
}

#***************************************************************************
# Subroutine:  set_table_as_track_source
# Description: set a DIGS screening DB table as a track source
#***************************************************************************
sub set_table_as_track_source {

	my ($self) = @_;


	die;

	# Get DIGS object and DB
	my $digs_obj = $self->{digs_obj};
	my $db = $digs_obj->{db};  
	
	# Get the ancillary tables in this DB
	my %extra_tables;
	my @extra_tables = qw [ digs_results ];
	$db->get_ancillary_table_names(\@extra_tables);
		
	my $table_num = 0;
	foreach my $table_name (@extra_tables) {
		$table_num++;
		$extra_tables{$table_num} = $table_name;
		print "\n\t\t Table $table_num: '$table_name'";
	}

	my $table_to_use;	
	my $question = "\n\n\t Use which table as a source";
	my $answer   = $console->ask_list_question($question, $table_num);
	$table_to_use = $extra_tables{$answer};
	unless ($table_to_use) { die; }
}

############################################################################
# INITIALISATION: load taxonomy & gene name translation tables
############################################################################

#***************************************************************************
# Subroutine:  load_taxonomy_table
# Description: load a table that captures taxonomic definitions
# Reason: tables allow IDs to be based on alternative taxonomic levels
#***************************************************************************
sub load_taxonomy_table {

	my ($self) = @_;

	my $table_name;
	print "\n\n\t  #≈#~# 3: LOAD TAXONOMY TRANSLATION TABLE\n";
	$table_name = 'taxonomy_translations';
	$self->{taxonomy_table_name} = $table_name;
	my %translations;
	#$self->load_taxonomy_translations_from_table(\%translations);
	#$self->{taxonomy_translations} = \%translations;

	return $table_name;
}

#***************************************************************************
# Subroutine:  load_gene_name_translation_tables
# Description: load a table that captures alternative names of homologous genes
#***************************************************************************
sub load_gene_name_translation_table {

	my ($self) = @_;

	my $table_name;
	print "\n\n\t  #≈#~# 4: LOAD GENE NAME TRANSLATION TABLE\n";
	$table_name = 'gene_translations';
	$self->{gene_names_table_name} = $table_name;
	my %translations;
	#$self->load_genename_translations_from_table(\%translations);
	#$self->{gene_name_translations} = \%translations;

	return $table_name;

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
	my $nom_chains_exists = $db_ref->does_table_exist('nomenclature_chains');
	unless ($nom_chains_exists) {
		$db_ref->create_nomenclature_chains_table($dbh);
	}

	# Load nomenclature tables
	print "\n\n\t  #        Loading nomenclature table";
	$db_ref->load_nomenclature_table($dbh);	
	print "\n\t  #        Loading nomenclature chains table";
	$db_ref->load_nomenclature_chains_table($dbh);

	# Flush nomenclature table
	my $nom_table    = $db_ref->{nomenclature_table};
	print "\n\t  #        Flushing nomenclature table";
	$nom_table->flush();

	# Flush nomenclature chains table
	my $chains_table = $db_ref->{nomenclature_chains_table};
	print "\n\t  #        Flushing nomenclature chains table";
	$chains_table->flush();

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


############################################################################
# UTILITY / DEV
############################################################################

#***************************************************************************
# Subroutine:  show_nomenclature_chains
# Description: Show nomenclature chains for all consolidated annotations
#***************************************************************************
sub show_nomenclature_chains {
	
	my ($self, $digs_obj) = @_;

	# Get relevant variables and objects
	unless ($digs_obj) { die; } # Sanity checking
	my $db = $digs_obj->{db};
	unless ($db) { die; } # Sanity checking

	# Get relevant variables and objects
	my $dbh = $db->{dbh};
	$db->load_nomenclature_table($dbh);
	$db->load_nomenclature_chains_table($dbh);
	my $nomenclature_table = $db->{nomenclature_table}; 
	my $chains_table  = $db->{nomenclature_chains_table};
	unless ($nomenclature_table and $chains_table) { die; } # Sanity checking

	# Get all named loci
	my $nom_where = " ORDER BY record_id ";
	my @loci;
	my @fields = qw [ record_id full_id ];
	$nomenclature_table->select_rows(\@fields, \@loci, $nom_where);	 
	
	# Iterate through loci
	foreach my $locus_ref (@loci) {

		my $locus_id = $locus_ref->{record_id};
		print "\n\t ### Chain $locus_id: ";	
		my $chain_where = " WHERE nomenclature_locus_id = $locus_id ";
		my @results;
		@fields = qw [ record_id track_id nomenclature_locus_id ];
		$chains_table->select_rows(\@fields, \@results, $chain_where);	 
		
		foreach my $result_ref (@results) {

			my @digs_results;
			my $track_entry_id = $result_ref->{track_id};
			my @fields = qw [ assigned_name assigned_gene ];
			my $where  = " WHERE record_id = $track_entry_id ";
			$nomenclature_table->select_rows(\@fields, \@digs_results, $where);
			
			foreach my $result_ref (@digs_results) {
		
				my $assigned_name = $result_ref->{assigned_name};
				my $assigned_gene = $result_ref->{assigned_gene};
				print " $track_entry_id ";
	
			}
		}
	}
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