#!usr/bin/perl -w
############################################################################
# Module:      DIGS.pm   database-integrated genome screening (DIGS)
# Description: Functions for implementing DIGS
# History:     December  2013: Created by Robert Gifford 
############################################################################
package DIGS;

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

# Maximum range for defragment
my $maximum   = 100000000;

# Flags
my $verbose   = undef;
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create new DIGS 'object'
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {
		
		# Flags
		process_id             => $parameter_ref->{process_id},
		program_version        => $parameter_ref->{program_version},
		
		# Paths and member variables
		genome_use_path        => $parameter_ref->{genome_use_path},
		output_path            => $parameter_ref->{output_path},

		# Database variables
		db_name                => '',   # Obtained from control file
		server                 => '',   # Obtained from control file
		username               => '',   # Obtained from control file
		password               => '',   # Obtained from control file
	
		# Member classes 
		blast_obj              => $parameter_ref->{blast_obj},
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# MAIN LOOP
############################################################################

#***************************************************************************
# Subroutine:  run_digs_process
# Description: handler for main DIGS functions 
#***************************************************************************
sub run_digs_process {

	my ($self, $option, $ctl_file) = @_;

 	# Show title
	$self->show_title();  

	# Hand off to functions 
	if ($option eq 1) { # Format targets files for BLAST searching
		my $target_db_obj = TargetDB->new($self);
		$target_db_obj->format_targets_for_blast();
	}
	else {
	
		# Initialise (an infile must be defined)
		unless ($ctl_file) { die "\n\t Option '$option' requires an infile\n\n"; }
		$self->initialise($ctl_file);
	
		# Load/create the screening
		$self->load_screening_db($ctl_file);

		# Hand off to DIGS functions
		if ($option eq 2) { # Screen;
			$self->setup_digs();
			$self->do_digs();	
		}
		elsif ($option eq 3) { # Reassign data in digs_results table
			my @digs_results;
			$self->initialise_reassign(\@digs_results); # Set up 
			$self->reassign(\@digs_results);	
		}
		elsif ($option eq 4) { # Interactively defragment results 	
			$self->interactive_defragment();	
		}
		elsif ($option eq 5) { # Consolidate results into higher level structures 
			$self->consolidate_loci();
		}
		elsif ($option eq 6) { # Standardised locus naming
			# Initialise (an infile must be defined)
			$self->create_standard_locus_ids();
		}
		else {
			print "\n\t  Unrecognized option '-m=$option'\n";
		}
	}
}

############################################################################
# PRIMARY FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  do_digs
# Description: do the core database-integrated genome screening processes
#***************************************************************************
sub do_digs {

	my ($self, $mode) = @_;

	# Get relevant member variables and objects
	my $queries_ref = $self->{queries};
	my $loader_obj  = $self->{loader_obj};
	unless ($queries_ref and $loader_obj) { die; }  # Sanity checking

	# Iterate through and excute screens
	print "\n\t  Starting database-integrated genome screening\n";
	my @probes = keys %$queries_ref;

	my $completed = 0;
	my %crossmatching;
	$self->{crossmatching} = \%crossmatching;
	foreach my $probe_name (@probes) {
		
		# Get the array of queries for this target file
		my $probe_queries = $queries_ref->{$probe_name};
		foreach my $query_ref (@$probe_queries) {  
	
			# Increment query count
			$completed++;		
			$self->{completed} = $completed;

			# Do the 1st BLAST (probe vs target)
			$self->search($query_ref);
		
			# Compress DB
			my @new_hits;
			$self->compress_active_set($query_ref, \@new_hits);

			# Extract newly identified or extended sequences
			my @extracted;
			my $target_path = $query_ref->{target_path};
			$self->extract($target_path, \@new_hits, \@extracted);	
	
			# Do the 2nd BLAST (hits from 1st BLAST vs reference library)
			my $assigned_count   = 0;
			my $crossmatch_count = 0;
			foreach my $hit_ref (@extracted) { # Iterate through the matches
				# Execute the 'reverse' BLAST (2nd BLAST in a round of paired BLAST)				
				my $assigned = $self->do_blast_genotyping($hit_ref);
				if ($assigned) { $assigned_count++; }
				# Get the unique key for this probe
				my $probe_name  = $query_ref->{probe_name};
				my $probe_gene  = $query_ref->{probe_gene};
				my $probe_key = $probe_name . '_' . $probe_gene; 		
				# Record cross-matching
				if ($probe_key ne $assigned) {
					$crossmatch_count++;
					$self->update_cross_matching($probe_key, $assigned);
				}
			}	
			if ($assigned_count > 0) {
				print "\n\t\t # $assigned_count extracted sequences matched to reference library";
				print "\n\t\t # $crossmatch_count cross-matched to something other than the probe";
			}
			
			# Update DB
			$self->update_db(\@extracted);
			
			# Show progress
			$self->show_digs_progress();
		}	
	}
	
	# Cleanup
	my $output_dir = $loader_obj->{report_dir};
	my $command1 = "rm -rf $output_dir";
	system $command1;

	# Show cross matching at end if verbose output setting is on
	if ($verbose) { $self->show_cross_matching(); }

	# Print finished message
	print "\n\n\t ### SCREEN COMPLETE ~ + ~ + ~";
}

#***************************************************************************
# Subroutine:  reassign
# Description: reassign sequences in the digs_results_table (for use after
#              the reference library has been updated)
#***************************************************************************
sub reassign {
	
	my ($self, $extracted_seqs_ref) = @_;

	# Get data structures and variables from self
	my $blast_obj       = $self->{blast_obj};
	my $result_path     = $self->{report_dir};
	my $db              = $self->{db};
	my $digs_results_table = $db->{digs_results_table};
	unless ($digs_results_table) { die; }
	
	# Iterate through the matches
	print "\n\n\t  Reassigning hits in the digs_results table\n";
	my $count = 0;
	my %reassign_matrix;
	my %unique_keys;
	foreach my $hit_ref (@$extracted_seqs_ref) {

		# Set the linking to the BLAST result table
		my $record_id  = $hit_ref->{record_id};
		
		my $extract_start   = $hit_ref->{extract_start};
		my $extract_end     = $hit_ref->{extract_end};
		$hit_ref->{subject_start} = $extract_start;
		$hit_ref->{subject_end}   = $extract_end;
		delete $hit_ref->{extract_start};
		delete $hit_ref->{extract_end};
	
		$hit_ref->{organism} = $hit_ref->{organism} ;
	
		# Execute the 'reverse' BLAST (2nd BLAST in a round of paired BLAST)	
		my $previous_assign = $hit_ref->{assigned_name};
		my $previous_gene   = $hit_ref->{assigned_gene};
		#print "\n\t Redoing assign for record ID $blast_id assigned to $previous_assign";
		#print "\n\t coordinates: $extract_start-$extract_end";
		$self->do_blast_genotyping($hit_ref);
		
		$count++;
		if (($count % 100) eq 0) {
			print "\n\t  Checked $count rows";
		}

		my $assigned_name = $hit_ref->{assigned_name};
		my $assigned_gene = $hit_ref->{assigned_gene};
		if ($assigned_name ne $previous_assign or  $assigned_gene ne $previous_gene) {
			
			# Show the results
			print "\n\t ##### Reassigned $record_id from $previous_assign ($previous_gene)";
			print " to $assigned_name ($assigned_gene)";
			my $where = " WHERE record_id = $record_id ";
			
			# Update the matrix
			my $previous_key = $previous_assign . '_' . $previous_gene;
			my $assigned_key = $assigned_name . '_' . $assigned_gene;
			$unique_keys{$assigned_key} = 1;
			$unique_keys{$previous_key} = 1;
			if ($reassign_matrix{$previous_key}) {
				my $hash_ref = $reassign_matrix{$previous_key};
				if ($hash_ref->{$assigned_key}) {
					$hash_ref->{$assigned_key}++;
				}
				else {
					$hash_ref->{$assigned_key} = 1;
				}
			}
			else {
				my %hash;
				$hash{$assigned_key} = 1;
				$reassign_matrix{$previous_key} = \%hash; 
			}
			# Insert the data
			delete $hit_ref->{record_id}; 
			delete $hit_ref->{organism}; 
			$digs_results_table->update($hit_ref, $where);
		}
	}
	
	# Write out the cross-matching matrix
	$self->write_crossmatching_data_to_file(\%reassign_matrix, \%unique_keys);

	# Cleanup
	my $output_dir = $self->{report_dir};
	my $command1 = "rm -rf $output_dir";
	system $command1;
}

#***************************************************************************
# Subroutine:  interactive_defragment 
# Description: console driven menu options for interactive defragment
#***************************************************************************
sub interactive_defragment {

	my ($self) = @_;

	# Display current settings
	my $redundancy_mode  = $self->{redundancy_mode};
	my $defragment_range = $self->{defragment_range};
	my $genome_use_path  = $self->{genome_use_path};
	unless ($redundancy_mode and  $defragment_range ) { die; } 
	unless ($genome_use_path) { die; } 
	print "\n\n\t\t Current settings (based on control file)";
	print "\n\t\t redundancy mode:  $redundancy_mode";
	print "\n\t\t defragment_range: $defragment_range";

	# Get a list of all the target files from the screening DB
	my $db = $self->{db};
	my $digs_results_table = $db->{digs_results_table};
	my @fields = qw [ organism target_datatype target_version target_name ];
	my @targets;
	$digs_results_table->select_distinct(\@fields, \@targets);

	# Settings
	my $choice;
	my %cluster_info;
	$cluster_info{total_hits} = '0';
	$cluster_info{total_clusters} = '0';
	my $t_range;
	
	# Question loop
	do {
		my $question1 = "\n\n\t # Set the range for merging hits";
		$t_range = $console->ask_int_with_bounds_question($question1, $defragment_range, $maximum);		

		# Apply the settings
		$self->interactive_defragment_loop(\@targets, \%cluster_info, $t_range);
		my $total_hits     = $cluster_info{total_hits};
		my $total_clusters = $cluster_info{total_clusters};
	
		# Prompt for what to do next
		print "\n\t\t\t TOTAL HITS:     $total_hits";
		print "\n\t\t\t TOTAL CLUSTERS: $total_clusters ";

		print "\n\n\t\t Option 1: preview new parameters";
		print "\n\t\t Option 2: apply these parameters";
		print "\n\t\t Option 3: exit";
		my $list_question = "\n\n\t # Choose an option:";
		$choice = $console->ask_list_question($list_question, 3);
		
	} until ($choice > 1);
	
	if ($choice eq 2) { # Apply the changes

		# Create a backup of the current digs_results
		$db->backup_digs_results_table();

		foreach my $target_ref (@targets) {

			my $organism        = $target_ref->{organism};
			my $target_name     = $target_ref->{target_name};
			my $target_datatype = $target_ref->{target_datatype};
			my $target_version  = $target_ref->{target_version};
			my $target_path  = $genome_use_path;
			my $target_path .= "~/Genomes/Mammalia/$organism/$target_datatype/$target_version/$target_name";
		
			# Create the relevant set of previously extracted loci
			my @combined;
			my $where  = " WHERE organism      = '$organism' ";
			$where    .= " AND target_datatype = '$target_datatype' ";
			$where    .= " AND target_version  = '$target_version' ";
			$where    .= " AND target_name     = '$target_name' "; 

			$self->get_sorted_digs_results(\@combined, $where);
			my $num_hits = scalar @combined;
		
			# Compose clusters of overlapping/adjacent BLAST hits and extracted loci
			my %settings;
			my %target_defragmented;
			$settings{range} = $t_range;
			$settings{start} = 'extract_start';
			$settings{end}   = 'extract_end';
			$self->compose_clusters(\%target_defragmented, \@combined, \%settings);
			#$self->show_clusters(\%defragmented);  # Show clusters
		
			# Defragment the results using these settings
			my @defragmented;
			$self->compress(\%settings, \%target_defragmented, \@combined, \@defragmented);
		
			# Extract newly identified or extended sequences
			my @extracted;
			$self->extract($target_path, \@defragmented, \@extracted);	
			#$devtools->print_array(\@extracted); die;
				
			# Do the 2nd BLAST (hits from 1st BLAST vs reference library)
			my $assigned_count   = 0;
			my $crossmatch_count = 0;
			foreach my $hit_ref (@extracted) { # Iterate through the matches
				# Execute the 'reverse' BLAST (2nd BLAST in a round of paired BLAST)				
				my $assigned = $self->do_blast_genotyping($hit_ref);
				if ($assigned) { $assigned_count++; }
			}
			
			# Update DB                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
			$self->update_db(\@extracted);
		}
	}
	elsif ($choice eq 3) { print "\n"; exit; }
}

#***************************************************************************
# Subroutine:  extend_screening_db
# Description: console managemant of ancillary tables in the screening database
#***************************************************************************
sub extend_screening_db {

	my ($self) = @_;

	# Get database handle
	my $db = $self->{db};
	unless ($db) { die; }
	my $dbh = $db->{dbh};
	unless ($dbh) { die "\n\t Couldn't retrieve database handle \n\n"; }

	# Declare the variables & data structures we need
	my %extra_tables;
	my @extra_tables;
	my $table_to_use;
	my %fields;
	my @fields;
	my $anc_table;
	my $table_name;

	# Show the options
	my @choices = qw [ 1 2 3 4 ];
	print "\n\n\t\t 1. Create new ancillary table";
	print "\n\t\t 2. Append data to existing ancillary table";
	print "\n\t\t 3. Flush existing ancillary table and upload fresh data";
	print "\n\t\t 4. Drop an ancillary table\n";
	my $question4 = "\n\t Choose an option";
	my $answer4   = $console->ask_simple_choice_question($question4, \@choices);

	# Create new table
	if ($answer4 == '1') {	
		my $table_name_question = "\n\t What is the name of the new table?";
		$table_name = $console->ask_question($table_name_question);
	}
	# or choose one of the ancillary tables already in the DB
	else {

		# Get the ancillary tables in this DB
		$db->get_ancillary_table_names(\@extra_tables);
		
		my $table_num = 0;
		foreach my $table_name (@extra_tables) {
			$table_num++;
			$extra_tables{$table_num} = $table_name;
			print "\n\t\t Table $table_num: '$table_name'";
		}
		my @table_choices = keys %extra_tables;

		my $question5 = "\n\n\t Apply to which of the above tables?";
		my $answer5   = $console->ask_simple_choice_question($question5, \@table_choices);
		$table_to_use = $extra_tables{$answer5};
		unless ($table_to_use) { die; }
	}
		
	# Upload data to table
	my @data;
	unless ($answer4 eq 4) {

		# Try to read the tab-delimited infile
		print "\n\n\t #### WARNING: This function expects a tab-delimited data table with column headers!";
		my $question1 = "\n\n\t Please enter the path to the file with the table data and column headings\n\n\t";
		my $infile = $console->ask_question($question1);
		unless ($infile) { die; }
		my @infile;
		$fileio->read_file($infile, \@infile);

		my $line_number = 0;
		foreach my $line (@infile) {
			$line_number++;
			if     ($line =~ /^\s*$/)  { next; } # discard blank line
			elsif  ($line =~ /^\s*#/)  { next; } # discard comment line 
			unless ($line =~ /\t/)     { print "\n\t Incorrect formatting at line '$line_number'"; die; }
			push (@data, $line);
		}
		my $data = scalar @data;
		unless ($data) {
			die "\n\t Couldn't read input file\n\n";
		}
		
		my $header_row = shift @data;
		my @header_row = split ("\t", $header_row);		
		print "\n\n\t The following column headers (i.e. table fields) were obtained\n";
		my $i;
		foreach my $element (@header_row) {
			chomp $element;
			$i++;
			$element =~ s/\s+/_/g;
			if ($element eq '') { $element = 'EMPTY_COLUMN_' . $i; } 
			print "\n\t\t Column $i: '$element'";
			push (@fields, $element);
			$fields{$element} = "varchar";
		}
		
		# Prompt user - did we read the file correctly?
		my $question3 = "\n\n\t Is this correct?";
		my $answer3 = $console->ask_yes_no_question($question3);
		if ($answer3 eq 'n') { # Exit if theres a problem with the infile
			print "\n\t\t Aborted!\n\n\n"; exit;
		}
	}


	# Create table if first time
	if ($answer4 eq 1) {
		$table_to_use = $db->create_ancillary_table($table_name, \@fields, \%fields);	
	}

	# Get a reference to a table object for the ancillary table
	$anc_table = MySQLtable->new($table_to_use, $dbh, \%fields);
	$db->{$table_to_use} = $anc_table;
		
	if ($answer4 == '4') {	# Drop the ancillary table
		$db->drop_ancillary_table($table_to_use);
		return;
	}
	if ($answer4 eq 3)   {  # Flush the table if requested
		$anc_table->flush();
		$anc_table->reset_primary_keys();
	}
	
	my $row_count = 0;
	foreach my $line (@data) { # Add data to the table
		$row_count++;
		chomp $line;
		my %insert;
		my @elements = split ("\t", $line);
		my $column_num = 0;
		foreach my $field (@fields) {
			my $value = $elements[$column_num];
			$column_num++;
			my $type  = $fields{$column_num};
			if ($verbose) {
				print "\n\t Row count $row_count: uploading value '$value' to field '$field'";
			}
			unless ($value) { 
				$value = 'NULL';
			}
			$insert{$field} = $value;
		}
		$anc_table->insert_row(\%insert);
	}
}

#***************************************************************************
# Subroutine:  consolidate_loci
# Description: assemble digs_results rows into higher-order loci 
#***************************************************************************
sub consolidate_loci {

	my ($self) = @_;

	print "\n\n\t  ### Consolidating assigned extracted sequences into loci \n";

    # Set up for consolidation
	my @sorted;

	# DEBUG OPTIONS
	#my $where = " WHERE orientation = '+'";
	#my $where = " WHERE assigned_name = 'HERV-E'";
	#$self->get_sorted_digs_results(\@sorted, $where);

	# Get the 
	$self->get_sorted_digs_results(\@sorted);
	my $total_hits = scalar @sorted;
	print "\n\t...$total_hits individual hits";
	#$devtools->print_array(\@sorted); die; # Show loci 	

	# Create tables if they don't exist already
	my $db_ref = $self->{db};
	my $dbh = $db_ref->{dbh};
	my $loci_exists = $db_ref->does_table_exist('loci');
	unless ($loci_exists) {
		$db_ref->create_loci_table($dbh);
	}
	my $loci_chains_exists = $db_ref->does_table_exist('loci_chains');
	unless ($loci_chains_exists) {
		$db_ref->create_loci_chains_table($dbh);
	}
	
	$db_ref->load_loci_table($dbh);
	$db_ref->load_loci_chains_table($dbh);

	# Set up for consolidate
	my %settings;
	my %consolidated;
	my $range = $self->{consolidate_range};
	my $d_range = $self->{defragment_range};
	unless ($range) { 
		my $question1 = "\n\n\t # Set the range for consolidating digs results";
		$range = $console->ask_int_with_bounds_question($question1, $d_range, $maximum);		
	}
	
	# Settings
	$settings{range} = $range;
	$settings{start} = 'extract_start';
	$settings{end}   = 'extract_end';

	# Compose clusters of overlapping/adjacent BLAST hits and extracted loci
	$self->compose_clusters(\%consolidated, \@sorted, \%settings);
	my @cluster_ids  = keys %consolidated;
	my $num_clusters = scalar @cluster_ids;
	if ($total_hits > $num_clusters) {
		print "\n\t...compressed to $num_clusters overlapping/contiguous clusters";
	}
	
	#$self->show_clusters(\%consolidated);  die; # DEBUG- Show clusters 	
	my $loci_table        = $db_ref->{loci_table};
	my $loci_chains_table = $db_ref->{loci_chains_table};
	foreach my $cluster_id (@cluster_ids) {
	
		my $hits_ref = $consolidated{$cluster_id};
		my $cluster_size = scalar @$hits_ref;
		if ($cluster_size > 1) {
			#$self->show_cluster($hits_ref, $cluster_id);
		}	
		
		my @locus_structure;
		my $initialised = undef;
		my $lowest;
		my $highest;
		my $scaffold;
		my $version;
		my $orientation;
		my $datatype;
		my $organism;
		my $target_name;
		my $assigned_name;

		foreach my $hit_ref (@$hits_ref) {
			
			my $feature     = $hit_ref->{assigned_gene};
			my $start       = $hit_ref->{extract_start};
			my $end         = $hit_ref->{extract_end};
			$version        = $hit_ref->{target_version};
			$datatype       = $hit_ref->{target_datatype};
			$orientation    = $hit_ref->{orientation};
			$scaffold       = $hit_ref->{scaffold};
			$organism       = $hit_ref->{organism};
			$target_name    = $hit_ref->{target_name};

			unless ($feature and $orientation) { die; }
			if ($orientation eq '+') {
				push(@locus_structure, $feature);
			}
			elsif ($orientation eq '-') {
				unshift(@locus_structure, $feature);			
			}

			if ($initialised) {
				if ($end > $highest) {
					$highest = $end;
				}
				if ($start < $lowest) {
					$lowest = $start;					
				}
			}
			else {
				$highest = $end;
				$lowest = $start;
				$initialised = 'true';								
			}
		}
		
		my $locus_structure = join('-', @locus_structure);
		my %locus;
		$locus{organism}        = $organism;
		$locus{target_version}  = $version;
		$locus{target_name}     = $target_name;
		$locus{target_datatype} = $datatype;
		$locus{scaffold}        = $scaffold;
		$locus{orientation}     = $orientation;
		$locus{extract_start}   = $lowest;
		$locus{extract_end}     = $highest;
		$locus{assigned_name}   = $assigned_name;
		$locus{locus_structure} = $locus_structure;
		my $locus_id  = $loci_table->insert_row(\%locus);
		
		foreach my $hit_ref (@$hits_ref) {
			my $digs_result_id = $hit_ref->{record_id};
			my %chain_data;
			$chain_data{digs_result_id} = $digs_result_id;
			$chain_data{locus_id}       = $locus_id;
			$loci_chains_table->insert_row(\%chain_data);
		}		
	}
}

#***************************************************************************
# Subroutine:  create_standard_locus_ids
# Description: apply standard ids to a set of loci from one or more genomes
#***************************************************************************
sub create_standard_locus_ids {

	my ($self, $infile) = @_;

	# Create nomenclature tables if they don't exist already
	my $db_ref = $self->{db};
	my $dbh = $db_ref->{dbh};
	my $nomenclature_exists = $db_ref->does_table_exist('nomenclature');
	unless ($nomenclature_exists) {
		$db_ref->create_nomenclature_table($dbh);
	}
	# Load nomenclature table
	$db_ref->load_nomenclature_table($dbh);

	# Load tracks into table
	$self->load_nomenclature_tracks();
	
	# Load translation table
	my %translations;
	$self->load_translations(\%translations);

	# Cluster annotations
	#$self->create_annotation_clusters(\%params, \%output, \@new_track, \%translations);
	die;

	# Create the report
	#$self->write_report(\%params, \%output);
}

############################################################################
# INTERNAL FUNCTIONS: MAIN DIGS SCREENING LOOP
############################################################################

#***************************************************************************
# Subroutine:  search
# Description: execute a similarity search and parse the results
#***************************************************************************
sub search {
	
	my ($self, $query_ref) = @_;

	# Get relevant member variables and objects
	my $blast_obj    = $self->{blast_obj};
	my $tmp_path     = $self->{tmp_path};
	my $min_length   = $self->{seq_length_minimum};
	my $min_score    = $self->{bitscore_minimum};

	# Sanity checking
	unless ($blast_obj)       { die; } 
	unless ($tmp_path)        { die; } 

	# Get screening database table objects
	my $db_ref           = $self->{db};
	my $searches_table   = $db_ref->{searches_table};
	my $active_set_table = $db_ref->{active_set_table};
	unless ($db_ref)          { die; } 

	# Get query details
	my $probe_id        = $query_ref->{probe_id};
	my $probe_name      = $query_ref->{probe_name};
	my $probe_gene      = $query_ref->{probe_gene};
	my $probe_type      = $query_ref->{probe_type};
	my $probe_path      = $query_ref->{probe_path};
	my $organism        = $query_ref->{organism};
	my $version         = $query_ref->{target_version};
	my $datatype        = $query_ref->{target_datatype};
	my $target_name     = $query_ref->{target_name};
	my $target_path     = $query_ref->{target_path};
	my $blast_alg       = $query_ref->{blast_alg};
	my $result_file     = $tmp_path . "/$probe_id" . "_$target_name.blast_result.tmp";
	#unless ($probe_id, $probe_type) { die; }

	# Do BLAST similarity search
	my $completed = $self->{completed};	
	print "\n\n\t  $blast_alg: $completed: '$organism' ($version, $datatype)";
	print   "\n\t  target: '$target_name'";
	print   "\n\t  probe:  '$probe_id'";   
	$blast_obj->blast($blast_alg, $target_path, $probe_path, $result_file);
	
	# Extract the results from tabular format BLAST output
	my @hits;
	$blast_obj->parse_tab_format_results($result_file, \@hits);
	my $rm_command = "rm $result_file";
	system $rm_command; # Remove the result file
	# TODO: catch error from BLAST and don't update "Searches_performed" table	

	# Summarise raw results of BLAST search
	my $num_hits = scalar @hits;
	if ($num_hits > 0) {
		print "\n\t\t # $num_hits matches to probe: $probe_name, $probe_gene";
	}
	
	# Apply filters & store results
	my $num_retained_hits = 0;
	my $score_exclude_count = '0';
	my $length_exclude_count = '0';
	foreach my $hit_ref (@hits) {

		my $skip = undef;

		# Apply length cutoff
		if ($min_length) { # Skip sequences that are too short
			my $start  = $hit_ref->{aln_start};
			my $end    = $hit_ref->{aln_stop};
			if ($end - $start < $min_length) {  
				$skip = 'true';
				$length_exclude_count++; 				
			}
		}

		# Apply bitscore cutoff
		if ($min_score) { # Skip sequences that have too low bit scores
			my $query_score = $hit_ref->{bitscore};
			if ($query_score < $min_score) {  
				$skip = 'true';
				$score_exclude_count++;
			}
		}
	
		unless ($skip) {		
			# Insert values into 'active_set' table
			$hit_ref->{digs_result_id}  = 0;
			$hit_ref->{organism}        = $organism;
			$hit_ref->{target_version}  = $version;
			$hit_ref->{target_datatype} = $datatype;
			$hit_ref->{target_name}     = $target_name;
			$hit_ref->{probe_id}        = $probe_id;
			$hit_ref->{probe_name}      = $probe_name;
			$hit_ref->{probe_gene}      = $probe_gene;
			$hit_ref->{probe_type}      = $probe_type;
			$hit_ref->{align_len}       = $hit_ref->{align_len};
			$hit_ref->{subject_start}   = $hit_ref->{aln_start}; 	# Rename to match DB
			$hit_ref->{subject_end}     = $hit_ref->{aln_stop}; 	# Rename to match DB
			$hit_ref->{query_end}       = $hit_ref->{query_stop}; # Rename to match DB
			$active_set_table->insert_row($hit_ref);
			$num_retained_hits++;			
		}
	} 

	# Show summary of BLAST results after filtering
	if ($score_exclude_count or $length_exclude_count) {
		print "\n\t\t # $num_retained_hits matches above threshold ";
		print "(excluded $length_exclude_count < length; $score_exclude_count < bitscore)";
	}

	# Update the searches table, to indicate this search has been performed
	$searches_table->insert_row($query_ref);
	
	return $num_hits;
}

#***************************************************************************
# Subroutine:  compress_active_set
# Description: determine what to extract based on current results
#***************************************************************************
sub compress_active_set {
	
	my ($self, $query_ref, $to_extract_ref) = @_;

	# Create the relevant set of previously extracted loci
	my @loci;
	my $where = $self->set_redundancy($query_ref);
	$self->get_sorted_digs_results(\@loci, $where);
	my $num_loci = scalar @loci;
	print "\n\t\t # $num_loci previously extracted loci";
		
	# Add relevant set of  previously extracted loci to the active set of BLAST results
	$self->create_combined_active_set(\@loci);

	# Get new BLAST results & previously extracted loci in a sorted list
	my @combined;
	$self->get_sorted_active_set(\@combined);

	my %settings;
	my %defragmented;
	$settings{range} = $self->{defragment_range};
	$settings{start} = 'subject_start';
	$settings{end}   = 'subject_end';
	$self->compress(\%settings, \%defragmented, \@combined, $to_extract_ref);
}

#***************************************************************************
# Subroutine:  compress
# Description: 
#***************************************************************************
sub compress {

	my ($self, $settings_ref, $defragmented_ref, $combined_ref, $to_extract_ref) = @_;

	# Compose clusters of overlapping/adjacent BLAST hits and extracted loci
	my $total_hits = scalar @$combined_ref;
	if ($total_hits > 0) {
		print "\n\t\t # $total_hits new hits & previously extracted loci ";
	}

	# Compose clusters of related sequences
	$self->compose_clusters($defragmented_ref, $combined_ref, $settings_ref);
	my @cluster_ids  = keys %$defragmented_ref;
	my $num_clusters = scalar @$combined_ref;
	if ($total_hits > $num_clusters) {
		print "...compressed to $num_clusters overlapping/contiguous clusters";
	}

	# Determine what to extract, and extract it
	# DEBUG $self->show_clusters($defragmented_ref);  # Show clusters
	$self->merge_clustered_loci($defragmented_ref, $to_extract_ref);
	my $num_new = scalar @$to_extract_ref;
	if ($num_new){
		print "\n\t\t # $num_new newly extracted sequences for assign/reassign";
	}	
}

#***************************************************************************
# Subroutine:  extract_sequences
# Description: extract sequences from target databases
#***************************************************************************
sub extract {

	my ($self, $target_path, $hits_ref, $extracted_ref) = @_;

	# Get paths, objects, data structures and variables from self
	my $blast_obj   = $self->{blast_obj};
	my $buffer      = $self->{extract_buffer};

	# Iterate through the list of sequences to extract
	my $new_hits = scalar @$hits_ref;
	my $i;
	foreach my $hit_ref (@$hits_ref) {
				
		# Add any buffer 
		$i++;
		#print "\n\t Extract $i";
		my $orientation   = $hit_ref->{orientation};
		if ($buffer) {
			if ($orientation eq '-') {
				$hit_ref->{start} = $hit_ref->{start} + $buffer;
				$hit_ref->{end}   = $hit_ref->{end} - $buffer;
				if ($hit_ref->{end} < 1) { # Don't allow negative coordinates
					$hit_ref->{end} = 1;
				}	
			}
			else {
				$hit_ref->{start} = $hit_ref->{start} - $buffer;
				if ($hit_ref->{start} < 1) { # Don't allow negative coordinates
					$hit_ref->{start} = 1;
				}	
				$hit_ref->{end}   = $hit_ref->{end} + $buffer;
			}
		}
		
		# Extract the sequence
		my $sequence = $blast_obj->extract_sequence($target_path, $hit_ref);
		if ($sequence) {	
			if ($verbose) {
				print "\t ........extracting";
			}
			my $seq_length = length $sequence; # Set sequence length
			$hit_ref->{extract_start}   = $hit_ref->{start};
			$hit_ref->{extract_end}     = $hit_ref->{end};
			$hit_ref->{sequence}        = $sequence;
			$hit_ref->{sequence_length} = $seq_length;
			#print "\n\t Sequence $sequence\n\n";	die;
			push (@$extracted_ref, $hit_ref);
		}
		else {
			if ($verbose) {
				print "\n\t Sequence extraction failed";
			}
		}
	}	
	return $new_hits;
}

#***************************************************************************
# Subroutine:  update_db
# Description: update the screening DB based on a completed round of DIGS
#***************************************************************************
sub update_db {

	my ($self, $extracted_ref) = @_;
		
	# Get parameters from self
	my $db_ref = $self->{db};
	my $digs_results_table  = $db_ref->{digs_results_table}; 
	my $active_set_table    = $db_ref->{active_set_table}; 
	my $blast_chains_table  = $db_ref->{blast_chains_table}; 

	# Flush the active set table
	#print "\n # Flushing 'active_set' table";
	$active_set_table->flush();
	

	# Iterate through the extracted sequences
	foreach my $hit_ref (@$extracted_ref) {

		# Insert the data to the Extracted_sequences table
		$hit_ref->{organism} = $hit_ref->{organism}; # Translate field name
		my $digs_result_id = $digs_results_table->insert_row($hit_ref);
		#$devtools->print_hash($hit_ref); die;
		
		# Insert the data to the BLAST_chains table
		my $blast_chains = $hit_ref->{blast_chains};
		#$devtools->print_hash($blast_chains);

		if ($blast_chains) {		
			my @blast_ids = keys %$blast_chains;
			foreach my $blast_id (@blast_ids) {							
				my $data_ref = $blast_chains->{$blast_id};
				$data_ref->{digs_result_id} = $digs_result_id;	
				#$devtools->print_hash($data_ref); exit;
				$blast_chains_table->insert_row($data_ref);
			}
		}

		# Delete superfluous data from the
		my $digs_result_ids_ref = $hit_ref->{digs_result_ids};
		foreach my $old_digs_result_id (@$digs_result_ids_ref) {			
			
			# Delete superfluous extract rows
			my $extracted_where = " WHERE record_id = $old_digs_result_id ";	
			$digs_results_table->delete_rows($extracted_where);
			# Update extract IDs			
			my $chains_where = " WHERE digs_result_id = $old_digs_result_id ";
			my %new_id;
			$new_id{digs_result_id} = $digs_result_id;	
			$blast_chains_table->update(\%new_id, $chains_where);
			#$devtools->print_hash($data_ref); exit;	
		}
	}
}

#***************************************************************************
# Subroutine:  show_digs_progress
# Description: show progress in DIGS screening
#***************************************************************************
sub show_digs_progress {

	my ($self) = @_;

	# Get the counts
	my $total_queries   = $self->{total_queries};
	my $completed       = $self->{completed};	
	unless ($completed and $total_queries) { die; } # Sanity checking
	
	# Calculate percentage progress
	my $percent_prog    = ($completed / $total_queries) * 100;
	my $f_percent_prog  = sprintf("%.2f", $percent_prog);
	print "\n\t\t  %$f_percent_prog completed";
}

#***************************************************************************
# Subroutine:  do_blast_genotyping
# Description: genotype a nucleotide sequence
#***************************************************************************
sub do_blast_genotyping {

	my ($self, $hit_ref) = @_;
	
	# Get paths and objects from self
	my $result_path   = $self->{tmp_path};
	my $blast_obj     = $self->{blast_obj};
	unless ($result_path and $blast_obj) { die; }
	
	# Get required data about the hit, prior to performing reverse BLAST
	my $sequence      = $hit_ref->{sequence};
	my $organism      = $hit_ref->{organism};
	my $probe_type    = $hit_ref->{probe_type};
	
	# Sanity checking
	unless ($organism)   { die; }
	unless ($probe_type) { die; }
	unless ($sequence)   { die "\n\t # ERROR: No sequence found for reverse BLAST"; } 

	# Make a FASTA query file for the reverse BLAST procedure
	$sequence =~ s/-//g;   # Remove any gaps that might happen to be there
	$sequence =~ s/~//g;   # Remove any gaps that might happen to be there
	$sequence =~ s/\s+//g; # Remove any gaps that might happen to be there
	my $fasta      = ">TEMPORARY\n$sequence";
	my $query_file = $result_path . 'TEMPORARY.fas';
	$fileio->write_text_to_file($query_file, $fasta);
	my $result_file = $result_path . 'TEMPORARY.blast_result';
		
	# Do the BLAST according to the type of sequence (AA or NA)
	my $blast_alg;
	my $lib_path;
	if ($probe_type eq 'UTR') {
		$lib_path  = $self->{blast_utr_lib_path};
		$blast_alg = 'blastn';
		unless ($lib_path) { die "\n\t NO UTR LIBRARY defined"; }
	}
	elsif ($probe_type eq 'ORF') {
		$lib_path  = $self->{blast_orf_lib_path};
		unless ($lib_path) {  die "\n\t NO ORF LIBRARY defined"; }
		$blast_alg = 'blastx';
	}
	else { die; }
	unless ($lib_path) { return; }

	# Execute the 'reverse' BLAST (2nd BLAST in a round of paired BLAST)	
	$blast_obj->blast($blast_alg, $lib_path, $query_file, $result_file);
	my @results;
	$blast_obj->parse_tab_format_results($result_file, \@results);

	# Define some variables for capturing the result
	my $top_match = shift @results;
	my $query_start   = $top_match->{query_start};
	my $query_end     = $top_match->{query_stop};
	my $subject_start = $top_match->{aln_start};
	my $subject_end   = $top_match->{aln_stop};
	my $assigned_name = $top_match->{scaffold};	
	my $assigned;

	# Deal with a query that matched nothing in the 2nd BLAST search
	unless ($assigned_name) {	
		$hit_ref->{assigned_name}    = 'Unassigned';
		$hit_ref->{assigned_gene}    = 'Unassigned';
		$hit_ref->{identity}         = 0;
		$hit_ref->{bitscore}        = 0;
		$hit_ref->{evalue_exp}      = 0;
		$hit_ref->{evalue_num}      = 0;
		$hit_ref->{mismatches}       = 0;
		$hit_ref->{align_len}        = 0;
		$hit_ref->{gap_openings}     = 0;
		$hit_ref->{query_end}        = 0;
		$hit_ref->{query_start}      = 0;
		$hit_ref->{subject_end}      = 0;
		$hit_ref->{subject_start}    = 0;
		$assigned = undef;
	}
	else {	# Assign the extracted sequence based on matches from 2nd BLAST search

		# Split assigned to into (i) refseq match (ii) refseq description (e.g. gene)	
		my @assigned_name = split('_', $assigned_name);
		my $assigned_gene = pop @assigned_name;
		$assigned_name = join ('_', @assigned_name);
		$hit_ref->{assigned_name}    = $assigned_name;
		$hit_ref->{assigned_gene}    = $assigned_gene;
		$hit_ref->{identity}         = $top_match->{identity};
		$hit_ref->{bitscore}        = $top_match->{bitscore};
		$hit_ref->{evalue_exp}      = $top_match->{evalue_exp};
		$hit_ref->{evalue_num}      = $top_match->{evalue_num};
		$hit_ref->{mismatches}       = $top_match->{mismatches};
		$hit_ref->{align_len}        = $top_match->{align_len};
		$hit_ref->{gap_openings}     = $top_match->{gap_openings};
		$hit_ref->{query_end}        = $query_end;
		$hit_ref->{query_start}      = $query_start;
		$hit_ref->{subject_end}      = $subject_end;
		$hit_ref->{subject_start}    = $subject_start;
		$assigned = $assigned_name . '_' . $assigned_gene;
	}

	# Clean up
	my $command1 = "rm $query_file";
	my $command2 = "rm $result_file";
	system $command1;
	system $command2;

	unless ($assigned) { $assigned = 'Unassigned'; }
	return $assigned;
}

#***************************************************************************
# Subroutine:  set_redundancy
# Description: compose SQL WHERE statement based on redundancy settings
#***************************************************************************
sub set_redundancy {
	
	my ($self, $query_ref) = @_;
	
	# Get query details
	my $probe_name      = $query_ref->{probe_name};
	my $probe_gene      = $query_ref->{probe_gene};
	my $probe_type      = $query_ref->{probe_type};
	my $target_name     = $query_ref->{target_name};
	my $organism        = $query_ref->{organism};
	my $target_datatype = $query_ref->{target_datatype};

	# Compose the WHERE statement based on redundancy settings
	my $redundancy_mode = $self->{redundancy_mode};
	unless ($redundancy_mode) { die; } 
	my $where  = " WHERE organism = '$organism' ";
	   $where .= " AND target_name = '$target_name' "; # Always limit by target
	if ($redundancy_mode eq 2 or $redundancy_mode eq 3) {
		# Mode 2 & 3, select all Extracted table rows with same value for 'assigned_gene'
		$where .= " AND assigned_gene = '$probe_gene' ";
	}
	if ($redundancy_mode eq 3) {
		# Mode 3, select all rows with same value for 'assigned_name' and 'assigned_gene'
		$where .= " AND assigned_name = '$probe_name' ";
	}

	return $where;
}

############################################################################
# INTERNAL FUNCTIONS: comparing and merging overlapping/adjacent hits
############################################################################

#***************************************************************************
# Subroutine:  compose_clusters 
# Description:
#***************************************************************************
sub compose_clusters {

	my ($self, $defragmented_ref, $hits_ref, $settings_ref) = @_;
	
	# Get settings
	my $range       = $settings_ref->{range};
	my $start_token = $settings_ref->{start};
	my $end_token   = $settings_ref->{end};
	
	# Iterate through consolidating as we go
	my $j = 1;
	my %last_hit;
	my %name_counts;
	my $initialised = undef;
	foreach my $hit_ref (@$hits_ref)  {

		# Get hit values
		my $record_id     = $hit_ref->{record_id};
		my $scaffold      = $hit_ref->{scaffold};
		my $target_name   = $hit_ref->{target_name};
		my $assigned_name = $hit_ref->{assigned_name};
		my $orientation   = $hit_ref->{orientation};
		my $start         = $hit_ref->{$start_token};
		my $end           = $hit_ref->{$end_token};
		#print "\n\t\t Range $start-$end on scaffold $scaffold";
		
		# Get last hit values
		my $last_record_id     = $last_hit{record_id};
		my $last_scaffold      = $last_hit{scaffold};
		my $last_target_name   = $last_hit{target_name};
		my $last_assigned_name = $last_hit{assigned_name};
		my $last_orientation   = $last_hit{orientation};
		my $last_start         = $last_hit{$start_token};
		my $last_end           = $last_hit{$end_token};		

	    if ($initialised) {
			# Sanity checking - are sequences in sorted order for this scaffold?
			if ( $scaffold eq $last_scaffold) {
				unless ($start >= $last_start) { 
					print "\n\t\t $start is less than $last_start";
					print " (end $end , last end $last_end) on scaffold $scaffold";
					die;
				}
			}			
            
            # Work out wether to merge this hit with the last
            my $new = $self->compare_adjacent_hits($hit_ref, \%last_hit, $settings_ref);
            
            if ($new) {
                 # Increment the count
                $j++;       
                # Initialise record
                $self->initialise_cluster($defragmented_ref, $hit_ref, $j);
            }
            else {             
                # Extend record
                $self->extend_cluster($defragmented_ref, $hit_ref, $j);
            }
        }
		else {
            $initialised = 'true'; # Set flag - we have seen at least one   
            $self->initialise_cluster($defragmented_ref, $hit_ref, $j);
		}

		# Update last hit data
		$last_hit{record_id}     = $record_id;
		$last_hit{assigned_name} = $assigned_name;
		$last_hit{scaffold}      = $scaffold;
		$last_hit{orientation}   = $orientation;
		$last_hit{$start_token}  = $start;
		$last_hit{$end_token}    = $end;
	}	
}

#***************************************************************************
# Subroutine:  compare_adjacent_hits
# Description: compare two hits 
#***************************************************************************
sub compare_adjacent_hits {

	my ($self, $hit_ref, $last_hit_ref, $settings_ref) = @_;

	# Get settings
	my $range       = $settings_ref->{range};
	my $start_token = $settings_ref->{start};
	my $end_token   = $settings_ref->{end};

	# Get the current hit values
	my $name           = $hit_ref->{assigned_name};
	my $scaffold       = $hit_ref->{scaffold};	
	my $start          = $hit_ref->{$start_token};
	my $end            = $hit_ref->{$end_token};
	my $orientation    = $hit_ref->{orientation};			

	# Get the last hit values
	my $last_name        = $last_hit_ref->{assigned_name};
	my $last_scaffold    = $last_hit_ref->{scaffold};	
	my $last_start       = $last_hit_ref->{$start_token};
	my $last_end         = $last_hit_ref->{$end_token};
	my $last_orientation = $last_hit_ref->{orientation};			
	
	if ($scaffold ne $last_scaffold) {
		return 1;
	}
	
	if ($orientation ne $last_orientation) {
		return 1;
	}
	
	my $gap = $start - $last_end;		
	#print "\n\t #\t CALC: '$scaffold': '$start'-'$last_end' = $gap";

	# Test whether to combine this pair into a set
	if ($gap < $range) {  # Combine
		#print "\n\t #\t COMBINE";
        return 0;
	
	}
	else { # Don't combine
		return 1;
	}
}

#***************************************************************************
# Subroutine:  initialise_cluster
# Description: 
#***************************************************************************
sub initialise_cluster {

	my ($self, $defragmented_ref, $hit_ref, $count) = @_;

    # Get the current hit values
	#print "\n\t New record ($count)";
    my @array;
    my %hit = %$hit_ref;
    push (@array, \%hit);
    $defragmented_ref->{$count} = \@array;
}

#***************************************************************************
# Subroutine:  extend_cluster 
# Description: 
#***************************************************************************
sub extend_cluster {

	my ($self, $defragmented_ref, $hit_ref, $count) = @_;

    # Get the current hit values
	#print "\n\t Extending record ($count) ";
    #$devtools->print_hash($hit_ref); die;
    my $array_ref = $defragmented_ref->{$count};
    push (@$array_ref, $hit_ref);

}

#***************************************************************************
# Subroutine:  merge_clustered_loci
# Description: 
#***************************************************************************
sub merge_clustered_loci {
	
	my ($self, $defragmented_ref, $to_extract_ref) = @_;

	#$devtools->print_hash($defragmented_ref); die;

	# Get screening database table objects
	my $db_ref           = $self->{db};
	my $searches_table   = $db_ref->{searches_table};
	my $active_set_table = $db_ref->{active_set_table};	

	my $extend_count = '0';
	my @cluster_ids = keys %$defragmented_ref;
	foreach my $cluster_id (@cluster_ids) {

		# Get data for this cluster
		my $cluster_ref = $defragmented_ref->{$cluster_id};
		unless ($cluster_ref) { die; }
		
		# Determine what to extract for this cluster
		my %new_blast_chains;
		my %previous_digs_result_ids;
		my $highest_end   = undef;
		my $lowest_start  = undef;
		my $previous_digs_result_id = undef;
		my $target_name;		
		my $version;
		my $target_datatype;		
		my $scaffold;			
		my $orientation;
		my $organism;
		my $probe_type;
		foreach my $hit_ref (@$cluster_ref) {
					
			my $record_id    = $hit_ref->{record_id};					
			my $digs_result_id   = $hit_ref->{digs_result_id};					
			my $start        = $hit_ref->{extract_start};			
			my $end          = $hit_ref->{extract_end};
			unless ($start) {
		        $start = $hit_ref->{subject_start};
		    }
			unless ($end) {
		        $end   = $hit_ref->{subject_end};
		    }

			$target_name     = $hit_ref->{target_name};					
			$target_datatype = $hit_ref->{target_datatype};			
			$version         = $hit_ref->{target_version};
			$scaffold        = $hit_ref->{scaffold};			
			$orientation     = $hit_ref->{orientation};
			$organism        = $hit_ref->{organism};
			unless ($organism) {
		        $organism = $hit_ref->{organism};
		    }
			$probe_type      = $hit_ref->{probe_type};
			#print "\n\t\t : $organism $target_name $scaffold $orientation $start";
			unless ($target_name) {
				$devtools->print_hash($hit_ref); die;
			}
							
			# Check if this is a a previously extracted locus
			if ($digs_result_id) {
				if ($previous_digs_result_id) { # If this represents a merge of multiple extracted loci, record the ids
					unless ($digs_result_id eq $previous_digs_result_id) {
					}
				}
				$previous_digs_result_ids{$digs_result_id} = $hit_ref;
				$previous_digs_result_id = $digs_result_id;		
			}
			# Store this BLAST result as part of a chain if it is new
			else {
				my %data = %$hit_ref;
				$new_blast_chains{$record_id} = \%data; 				
			}
			
			#print "\n\t\t RECORD ID:\t $scaffold $extract_start-$extract_end";
			if ($lowest_start and $highest_end ) {		
				if ($start > $lowest_start) { $lowest_start = $start; }
				if ($end > $highest_end)    { $highest_end  = $end;   }
			}
			else {
				$highest_end  = $end;
				$lowest_start = $start;
			}									
		}

		# Determine whether or not to extract
		my $extract = undef;		

		# Is this cluster composed entirely of new loci?
		unless ($previous_digs_result_id) {
			$extract = 'true';				
		}
		# Is this a merge of multiple previously extracted loci?
		my @previous_digs_result_ids = keys %previous_digs_result_ids;
		my $num_previously_extracted_loci_in_cluster = scalar @previous_digs_result_ids;
		if ($num_previously_extracted_loci_in_cluster > 1) {
			my $combined = join (',', @previous_digs_result_ids);
			print "\n\t\t # Merging previously extracted loci: ($combined) ";
			$extract = 'true';							
		}
		# If it includes a single extracted locus, does this locus need to be extended?
		elsif ($previous_digs_result_id) {
		
			# get the indexed query 
			my $data_ref = $previous_digs_result_ids{$previous_digs_result_id};
			my $extract_start = $data_ref->{subject_start};
			my $extract_end   = $data_ref->{subject_end};
			unless ($lowest_start >= $extract_start and $highest_end <= $extract_end) {	
				$extend_count++;
				$extract = 'true';							
				#print "\n\t\t # Extending extracted sequence: $extract_start, $extract_end: ($lowest_start-$highest_end) ";
			}			
		}
		if ($extract) {
		
			# Set the extract params for this cluster
			my %extract;
			$extract{digs_result_id}      = $previous_digs_result_id;
			$extract{target_name}     = $target_name;
			$extract{target_datatype} = $target_datatype;
			$extract{target_version}  = $version;
			$extract{organism} = $organism;
			$extract{probe_type}      = $probe_type;
			$extract{digs_result_id}      = $previous_digs_result_id;
			$extract{target_name}     = $target_name;
			$extract{scaffold}        = $scaffold;
			$extract{start}           = $lowest_start;	
			$extract{end}             = $highest_end;
			$extract{orientation}     = $orientation;
			$extract{digs_result_ids}   = \@previous_digs_result_ids;
			my $num_chains = scalar keys %new_blast_chains;
			if ($num_chains) {
				$extract{blast_chains} = \%new_blast_chains;
			}
			push (@$to_extract_ref, \%extract);	
		}
	}
	print "\n\t\t # $extend_count extensions to previously extracted sequences ";
	
}

#***************************************************************************
# Subroutine:  show_clusters
# Description: 
#***************************************************************************
sub show_clusters {

	my ($self, $defragmented_ref) = @_;

	#$devtools->print_hash($defragmented_ref); die;

	my @cluster_ids = keys %$defragmented_ref;
	my $cluster_count;
	foreach my $id (@cluster_ids) {
		$cluster_count++;
		my $hits_ref = $defragmented_ref->{$id};
		my $cluster_size = scalar @$hits_ref;
		if ($cluster_size > 1) {
			$self->show_cluster($hits_ref, $cluster_count);
		}	
	}
}

#***************************************************************************
# Subroutine:  show_cluster
# Description: 
#***************************************************************************
sub show_cluster {

	my ($self, $hits_ref, $cluster_id) = @_;

 	#$devtools->print_array($hits_ref); die;	
	#print "\n";
	
	foreach my $hit_ref (@$hits_ref) {
   		
   		#$devtools->print_hash($hit_ref); die;	
		my $organism      = $hit_ref->{organism};			
		my $assigned_name = $hit_ref->{probe_name};			
		my $assigned_gene = $hit_ref->{probe_gene};			
		my $orientation   = $hit_ref->{orientation};			

		unless ($assigned_name) {
			$assigned_name = $hit_ref->{assigned_name};			
		}
		unless ($assigned_gene) {
			$assigned_gene = $hit_ref->{assigned_gene};			
		}

		my $scaffold      = $hit_ref->{scaffold};			
		my $start         = $hit_ref->{extract_start};			
		my $end           = $hit_ref->{extract_end};
		unless ($start) {
			$start = $hit_ref->{subject_start};			
		}
		unless ($end) {
			$end = $hit_ref->{subject_end};			
		}

		my $digs_result_id    = $hit_ref->{digs_result_id};

		print "\n\t\t CLUSTER $cluster_id $organism: ";
		print "$assigned_name: $assigned_gene: $scaffold $start-$end ($orientation)";
		if ($digs_result_id) {
			print " (extract ID: $digs_result_id)";				
		}
	}			
}

#***************************************************************************
# Subroutine:  interactive_defragment_loop 
# Description: 
#***************************************************************************
sub interactive_defragment_loop {

    my ($self, $targets_ref, $cluster_params, $t_range) = @_;
   
	# Apply the settings
	my $total_hits     = $cluster_params->{total_hits};
	my $total_clusters = $cluster_params->{total_clusters};
   
	foreach my $target_ref (@$targets_ref) {

		my $organism        = $target_ref->{organism};
		my $target_name     = $target_ref->{target_name};
		my $target_datatype = $target_ref->{target_datatype};
		my $target_version  = $target_ref->{target_version};
			
		# Create the relevant set of previously extracted loci
		my @loci;
		my $where  = " WHERE organism      = '$organism' ";
		$where    .= " AND target_datatype = '$target_datatype' ";
		$where    .= " AND target_version  = '$target_version' ";
		$where    .= " AND target_name     = '$target_name' "; 

		$self->get_sorted_digs_results(\@loci, $where);
		my $num_hits = scalar @loci;
		
		# Compose clusters of overlapping/adjacent BLAST hits and extracted loci
		my %settings;
		my %target_defragmented;
		$settings{range} = $t_range;
		$settings{start} = 'extract_start';
		$settings{end}   = 'extract_end';
		$self->compose_clusters(\%target_defragmented, \@loci, \%settings);
		
		# Show clusters
		my @cluster_ids  = keys %target_defragmented;
		my $num_clusters = scalar @cluster_ids;
		if ($verbose) {
			print "\n\n\t\t $num_hits hits in target $target_name";
			if ($num_hits > $num_clusters) {
				print "\n\t\t   > $num_clusters overlapping/contiguous clusters";
				$self->show_clusters(\%target_defragmented);
			}
		}
		$total_hits = $total_hits + $num_hits;
		$total_clusters = $total_clusters + $num_clusters;
	}

	$cluster_params->{total_hits}     = $total_hits;
	$cluster_params->{total_clusters} = $total_clusters;
		
}

############################################################################
# INTERNAL FUNCTIONS: recording cross-matching during DIGS
###########################################################################

#***************************************************************************
# Subroutine:  update_cross_matching
# Description: update a hash to record cross-matches
#***************************************************************************
sub update_cross_matching {

	my ($self, $probe_key, $assigned) = @_;
	
	my $crossmatch_ref = $self->{crossmatching};
	
	if ($crossmatch_ref->{$probe_key}) {
		my $cross_matches_ref = $crossmatch_ref->{$probe_key};
		if ($cross_matches_ref->{$assigned}) {
			$cross_matches_ref->{$assigned}++;
		}
		else {
			$cross_matches_ref->{$assigned} = 1;
		}
	}
	else {
		my %crossmatch;
		$crossmatch{$assigned} = 1;
		$crossmatch_ref->{$probe_key} = \%crossmatch;
	}
}

#***************************************************************************
# Subroutine:  show_cross_matching
# Description: show contents of hash that records cross-matches
#***************************************************************************
sub show_cross_matching {

	my ($self) = @_;

	print "\n\n\t  Summary of cross-matching";   
	my $crossmatch_ref = $self->{crossmatching};
	my @probe_names = keys 	%$crossmatch_ref;
	foreach my $probe_name (@probe_names) {
		
		my $cross_matches_ref = $crossmatch_ref->{$probe_name};
		my @cross_matches = keys %$cross_matches_ref;
		foreach my $cross_match (@cross_matches) {
			my $count = $cross_matches_ref->{$cross_match};
			print "\n\t\t #   $count x $probe_name to $cross_match";
		}
	}
}

#***************************************************************************
# Subroutine:  write_crossmatching_data_to_file 
# Description: write cross-matching results as a matrix
#***************************************************************************
sub write_crossmatching_data_to_file {

    my ($self, $matrix_ref, $keys_ref) = @_;

    my @matrix;
    my @keys = sort keys %$keys_ref;
	
	# Reformat the keys - split into the two separate elements
	my %reformatted;
    my @horizontal;
	foreach my $key (@keys) {
    
		my @key = split('_', $key);
		my $second_part = pop @key;
		my $first_part = shift @key;
		my $reformatted_key = $first_part . " ($second_part)";
		$reformatted{$key} = $reformatted_key;
		push (@horizontal, $reformatted_key);
	}

	# Create the header row for the matrix
    my $horizontal = join("\t", @horizontal);
	push (@matrix, "\t$horizontal\n");

	foreach my $key (@keys) {
		
		# Write the name
		my @line;
		my $f_key = $reformatted{$key};
		push (@line, $f_key);
		
		# Write the numbers (internal part of the matrix)
		foreach my $comparison_key (@keys) {
			
			my $count;
			my $hash_ref = $matrix_ref->{$key};
			$count = $hash_ref->{$comparison_key};
			unless ($count) { $count = '0'; }
			push (@line, $count);
		}
		# Create matrix row
    	my $line = join("\t", @line);
		push (@matrix, "$line\n");
	}

	my $file = 'reassign_matrix.txt';
	my $file_path = $self->{tmp_path} . $file;
	$fileio->write_file($file_path, \@matrix);
	#$devtools->print_hash($matrix_ref);
}

############################################################################
# INTERNAL FUNCTIONS: interacting with screening DB (indexing, sorting)
############################################################################

#***************************************************************************
# Subroutine:  index_previously_executed_searches 
# Description: index BLAST searches that have previously been executed
#***************************************************************************
sub index_previously_executed_searches {
	
	my ($self, $done_ref) = @_;

	my $db = $self->{db};	
	my $searches_table = $db->{searches_table};
	unless ($searches_table) { die "\n\t Searches_performed table not loaded\n\n"; }
	my @data;
	my @fields = qw [ record_id
                      probe_name probe_gene 
	                  organism target_datatype target_version target_name ];
	my $where = " ORDER BY record_id ";
	$searches_table->select_rows(\@fields, \@data, $where);
	
	# Index the executed searches
	foreach my $data_ref (@data) {
		
		# Get the query parameters
		my $organism        = $data_ref->{organism};
		my $target_datatype = $data_ref->{target_datatype};
		my $version         = $data_ref->{target_version};
		my $target_name     = $data_ref->{target_name};
		my $probe_name      = $data_ref->{probe_name};
		my $probe_gene      = $data_ref->{probe_gene};
	
		# Sanity checking
		unless ( $organism )        { die; }
        unless ( $target_datatype ) { die; }
        unless ( $version )         { die; }
        unless ( $target_name )     { die; }
        unless ( $probe_name )      { die; }
        unless ( $probe_gene )      { die; }
		
		# Create the unique key for this search
		my @genome = ( $organism , $target_datatype, $version );
		my $target_id = join ('|', @genome);
		my $probe_id  = $probe_name . '_' .  $probe_gene;
		my @key = ( $target_id, $target_name, $probe_id );
		my $key = join ('|', @key);

		# Record the query, indexed by it's unique key
		$done_ref->{$key} = $data_ref;		
	}

	#$devtools->print_hash($done_ref); die; # DEBUG
}

#***************************************************************************
# Subroutine:  get sorted digs results
# Description: 
#***************************************************************************
sub get_sorted_digs_results {

	my ($self, $data_ref, $where) = @_;;

	# Set statement to sort loci
	my $sort  = " ORDER BY target_name, scaffold, extract_start ";
	if ($where) { $where .= $sort; }
	else        { $where  = $sort; }

	# Get database tables
	my $db = $self->{db};
	my $digs_results_table  = $db->{digs_results_table};
		
	# Set the fields to get values for
	my @extract_fields = qw [ record_id organism 
	                          target_name target_version target_datatype
	                          assigned_name assigned_gene probe_type
	                          scaffold orientation
	                          bitscore gap_openings
	                          query_start query_end 
	                          mismatches align_len
                              evalue_num evalue_exp identity 
                              extract_start extract_end sequence_length ];
	$digs_results_table->select_rows(\@extract_fields, $data_ref, $where);
}

#***************************************************************************
# Subroutine:  create_combined_active_set 
# Description: get all extracted loci for this target file (& probe)
#***************************************************************************
sub create_combined_active_set {
	
	my ($self, $data_ref) = @_;;

	# Get database tables
	my $db = $self->{db};
	my $active_set_table  = $db->{active_set_table};

	# Enter all relevant extracted loci into 'active_set' table 
	my $num_loci = scalar @$data_ref;
	foreach my $locus_ref (@$data_ref) {

		#print "\n\t\t # inserting extract ID $digs_result_id";
		#$devtools->print_hash($locus_ref);
		my $digs_result_id = $locus_ref->{record_id};
		
		# Translations
		$locus_ref->{digs_result_id}       = $digs_result_id;
		$locus_ref->{probe_name}       = $locus_ref->{assigned_name};
		$locus_ref->{probe_gene}       = $locus_ref->{assigned_gene};
		$locus_ref->{subject_start}    = $locus_ref->{extract_start};
		$locus_ref->{subject_end}      = $locus_ref->{extract_end};
		$locus_ref->{organism}  = $locus_ref->{organism};
		$active_set_table->insert_row($locus_ref);
	}
}

#***************************************************************************
# Subroutine:  get sorted active set 
# Description: 
#***************************************************************************
sub get_sorted_active_set {
	
	my ($self, $data_ref, $where) = @_;;

	# Set statement to sort loci
	if ($where) { $where .= " ORDER BY scaffold, subject_start "; }
	else        { $where  = " ORDER BY scaffold, subject_start "; }

	# Get database tables
	my $db = $self->{db};
	my $active_set_table    = $db->{active_set_table};
	
	# Get sorted, combined extracted loci and new blast results	
	my @blast_fields = qw [ record_id digs_result_id  
	                        organism target_datatype target_version target_name
	                        probe_name probe_gene probe_type
	                        bitscore gap_openings
	                        query_start query_end 
	                        align_len mismatches 
                            evalue_num evalue_exp identity 
	                        scaffold orientation
	                        subject_start subject_end ];

	$active_set_table->select_rows(\@blast_fields, $data_ref, $where);

}

#***************************************************************************
# Subroutine:  show_blast_chains
# Description: Show BLAST chains for all extracted sequences
#***************************************************************************
sub show_blast_chains {
	
	my ($self) = @_;

	# Get relevant variables and objects
	my $db = $self->{db};
	unless ($db) { die; } # Sanity checking
	my $digs_results_table = $db->{digs_results_table}; 
	my $blast_chains_table = $db->{blast_chains_table};
	my $extract_where = " ORDER BY record_id ";
	my @extracted_ids;
	my @fields = qw [ record_id assigned_name assigned_gene ];
	$digs_results_table->select_rows(\@fields, \@extracted_ids, $extract_where);	 
	
	foreach my $hit_ref (@extracted_ids) {
		my $digs_result_id = $hit_ref->{record_id};
		my @chain;
		my $assigned_name = $hit_ref->{assigned_name};
		my $assigned_gene = $hit_ref->{assigned_gene};
		my @chain_fields = qw [ record_id probe_name probe_gene 
		                        organism target_name 
		                        scaffold subject_start subject_end
		                        bitscore identity align_len ];
		my $blast_where  = " WHERE Extract_ID = $digs_result_id ";
		   $blast_where .= " ORDER BY subject_start";
		$blast_chains_table->select_rows(\@chain_fields, \@chain, $blast_where);	
		print "\n\t ### BLAST hit chain for extracted locus $digs_result_id";
		print " ($assigned_name, $assigned_gene):";
		foreach my $hit_ref (@chain) {
			my $blast_id    = $hit_ref->{record_id};
			my $probe_name  = $hit_ref->{probe_name};
			my $probe_gene  = $hit_ref->{probe_gene};
			my $organism    = $hit_ref->{organism};
			my $scaffold    = $hit_ref->{scaffold};
			my $start       = $hit_ref->{subject_start};
			my $end         = $hit_ref->{subject_end};
			my $bitscore    = $hit_ref->{bitscore};
			my $identity    = $hit_ref->{identity};
			my $f_identity  = sprintf("%.2f", $identity);
			my $align_len   = $hit_ref->{align_len};
			print "\n\t\t $blast_id:\t Score: $bitscore, \%$f_identity identity ";
			print "across $align_len aa ($start-$end) to:\t $probe_name ($probe_gene) ";
		}
	}
}

#***************************************************************************
# Subroutine:  show_locus_chains
# Description: Show digs_result chains for all consolidated loci
#***************************************************************************
sub show_locus_chains {
	
	my ($self) = @_;

	# Get relevant variables and objects
	my $db = $self->{db};
	unless ($db) { die; } # Sanity checking
	my $dbh = $db->{dbh};
	$db->load_loci_table($dbh);
	$db->load_loci_chains_table($dbh);
	#$devtools->print_hash($db); die;

	my $digs_results_table = $db->{digs_results_table}; 
	my $loci_table         = $db->{loci_table};
	my $loci_chains_table  = $db->{loci_chains_table};
	unless ($digs_results_table and $loci_chains_table) {
		print "\n\t # Locus tables not found - run consolidate first\n\n\n";
		exit;
	}

	# Get all loci
	my $loci_where = " ORDER BY record_id ";
	my @loci;
	my @fields = qw [ record_id locus_structure ];
	$loci_table->select_rows(\@fields, \@loci, $loci_where);	 
	
	# Iterate through loci
	foreach my $locus_ref (@loci) {

		my $locus_id = $locus_ref->{record_id};
		print "\n\t ### Chain $locus_id ";	
		my $chain_where = " WHERE locus_id = $locus_id ";
		my @results;
		my @fields = qw [ record_id locus_id digs_result_id ];
		$loci_chains_table->select_rows(\@fields, \@results, $chain_where);	 
		
		foreach my $result_ref (@results) {

			my @digs_results;
			my $digs_result_id = $result_ref->{digs_result_id};
			my @result_fields = qw [ assigned_name assigned_gene 
			                         organism target_name 
			                         scaffold extract_start extract_end
			                         bitscore identity align_len ];
			my $where  = " WHERE record_id = $digs_result_id ";
			$digs_results_table->select_rows(\@result_fields, \@digs_results, $where);			

			foreach my $result_ref (@digs_results) {
		
				my $assigned_name = $result_ref->{assigned_name};
				my $assigned_gene = $result_ref->{assigned_gene};
				my $organism      = $result_ref->{organism};	
				my $scaffold      = $result_ref->{scaffold};
				my $start         = $result_ref->{extract_start};
				my $end           = $result_ref->{extract_end};
				my $bitscore      = $result_ref->{bitscore};
				my $identity      = $result_ref->{identity};
				my $align_len     = $result_ref->{align_len};
				my $f_identity    = sprintf("%.2f", $identity);
				print "\n\t\t $digs_result_id:\t Score: $bitscore, \%$f_identity identity ";
				print "across $align_len aa ($start-$end) to:\t $assigned_name ($assigned_gene) ";
	
			}
		}
	}
}

############################################################################
# INTERNAL FUNCTIONS: nomenclature
############################################################################

#***************************************************************************
# Subroutine:  extract_values_from_track_line  
# Description: 
#***************************************************************************
sub extract_values_from_track_line {

    my ($self, $line, $hit_ref) = @_;

    # Chomp and split the line on tabs
    chomp $line;
	my @line = split("\t", $line);
    
    # Get values
    my $track     = shift @line;
    my $name      = shift @line;
    my $scaffold  = shift @line;
    my $start     = shift @line;
    my $end       = shift @line;
    my $gene      = shift @line;
    my $master_id = shift @line;
    my $orientation = get_orientation($start, $end);

    unless ($track) { die; }
    unless ($name)  { die; }
    unless ($start) { 
    	if ($start eq '0') {
    		$start = 1;
    	}
    	else { print "\n\n$line\n\n"; die; } 
    
    }
    unless ($end)   { die; }
    unless ($gene)  { print "\n\n$line\n\n"; die; }
    unless ($orientation)  { die; }

    # Extract line
    $hit_ref->{track}       = $track;
    $hit_ref->{name}        = $name;
    $hit_ref->{scaffold}    = $scaffold;
    $hit_ref->{start}       = $start;
    $hit_ref->{end}         = $end;
    $hit_ref->{orientation} = $orientation;
    $hit_ref->{gene}        = $gene;
    $hit_ref->{master_id}   = $master_id;
    $hit_ref->{hit_id}      = $track . '_' . $name;
    $hit_ref->{line}        = $line;

}

#***************************************************************************
# Subroutine:  load_nomenclature_tracks  
# Description: 
#***************************************************************************
sub load_nomenclature_tracks {

    my ($self) = @_;
    
	# Get new track(s) into table
	my @new_track;
	my $new_track_path = $self->{new_track_path};
	unless ($new_track_path) { die; }
	$fileio->read_file($new_track_path, \@new_track);
	
}

#***************************************************************************
# Subroutine:  load_translations
# Description: load translation tables
#***************************************************************************
sub load_translations {

	my ($self, $translations_ref) = @_;

	# Read translation from file
	my $translations_path = $self->{translations_path};
	unless ($translations_path) { die; }
	my @file;
	read_file($translations_path, \@file);
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

############################################################################
# INTERNAL FUNCTIONS: initialisation
############################################################################

#***************************************************************************
# Subroutine:  initialise 
# Description: initialise module for interacting with screening database
#              and perform basic validation of options and input file 
#***************************************************************************
sub initialise {

	my ($self, $ctl_file) = @_;

	# Try opening control file
	my @ctl_file;
	my $valid = $fileio->read_file($ctl_file, \@ctl_file);
	unless ($valid) {  # Exit if we can't open the file
		die "\n\t ### Couldn't open control file '$ctl_file'\n\n\n ";
	}

	# If control file looks OK, store the path and parse the file
	$self->{ctl_file}   = $ctl_file;
	my $loader_obj = ScreenBuilder->new($self);
	$loader_obj->parse_control_file($ctl_file, $self);

	# Update numthreads setting in BLAST object
	my $num_threads = $loader_obj->{num_threads};
	unless ($num_threads) { $num_threads = 1; }  # Default setting
	$self->{blast_obj}->{num_threads} = $num_threads;

	# Store the ScreenBuilder object (used later)
	$self->{loader_obj} = $loader_obj;
	
}

#***************************************************************************
# Subroutine:  initialise_reassign 
# Description: set up for reassign process
#***************************************************************************
sub initialise_reassign {

	my ($self, $extracted_seqs_ref) = @_;

	# Create a unique ID and report directory for this run
	my $output_path = $self->{output_path};
	my $process_id  = $self->{process_id};
	my $db          = $self->{db};
	my $db_name     = $db->{db_name};
	unless ($db and $db_name and $process_id and $output_path) { die; }
	#$devtools->print_hash($self); die;
	
	# Create report directory
	my $loader_obj = $self->{loader_obj};
	$loader_obj->create_output_directories($self);

	# Set up the reference library
	$loader_obj->setup_reference_library($self);
	my $where;
	if ($self->{blast_utr_lib_path}) {
		$where = " WHERE probe_type = 'UTR'";
	}
	elsif ($self->{blast_orf_lib_path}) {
		$where = " WHERE probe_type = 'ORF'";
	}

	# Get the assigned data
	my $digs_results_table = $db->{digs_results_table};
	my @fields  = qw [ record_id probe_type 
	                   assigned_name assigned_gene
	                   organism 
	                   extract_start extract_end sequence ];
	$digs_results_table->select_rows(\@fields, $extracted_seqs_ref, $where);

}

#***************************************************************************
# Subroutine:  setup_digs
# Description: prepare database and DIGS query list to commence screening
#***************************************************************************
sub setup_digs {

	my ($self) = @_;

	# Flush active set
	my $db  = $self->{db};
	unless ($db) { die "\n\t Error: no DB defined \n\n\n"; }

	#print "\n\t  Flushing 'active_set' table\n";
	my $active_set_table = $db->{active_set_table};
	$active_set_table->flush();
	
	# Index previously executed searches
	my %done;
	$self->index_previously_executed_searches(\%done);
	
	# Set up the screening queries
	my %queries;
	my $loader_obj = $self->{loader_obj};
	unless ($loader_obj) { die; }  # Sanity checking
	#$devtools->print_hash(\%done);
	$loader_obj->{previously_executed_searches} = \%done;

	my $total_queries = $loader_obj->set_up_screen($self, \%queries);
	unless ($total_queries)  { 
		print "\n\t  Exiting without screening.\n\n";	
		exit;
	}
		
	# Record queries 
	$self->{queries}       = \%queries;
	$self->{total_queries} = $total_queries;
}

#***************************************************************************
# Subroutine:  load_screening_db
# Description: load a DIGS screening database, create if doesn't exist 
#***************************************************************************
sub load_screening_db {

	my ($self, $ctl_file) = @_;

	# Get required objects and info from self, check everything looks OK
	my $loader_obj = $self->{loader_obj};
	unless ($loader_obj) { die; } 
	my $db_name = $loader_obj->{db_name};
	unless ($db_name) { die "\n\t Error: no DB name defined \n\n\n"; }

	# Create the screening DB object
	my $db_obj = ScreeningDB->new($loader_obj);

	# Check if this screening DB exists, if not then create it
	my $db_exists = $db_obj->does_db_exist($db_name);
	unless ($db_exists) {
		$db_obj->create_screening_db($db_name);	
	}
	
	# Load the table handles into screening database object
	$db_obj->load_screening_db($db_name);	
	$self->{db} = $db_obj; # Store the database object reference 
}

############################################################################
# INTERNAL FUNCTIONS: title and help display
############################################################################

#***************************************************************************
# Subroutine:  show_title
# Description: show command line title blurb 
#***************************************************************************
sub show_title {

	my ($self) = @_;

	my $version_num =  $self->{program_version};
	unless ($version_num) {
		$version_num = 'version undefined (use with caution)';
	}
	$console->refresh();
	my $title       = 'DIGS';
	my $description = 'Database-Integrated Genome Screening';
	my $author      = 'Robert J. Gifford';
	my $contact	    = '<robert.gifford@glasgow.ac.uk>';
	$console->show_about_box($title, $version_num, $description, $author, $contact);
}

#***************************************************************************
# Subroutine:  show_help_page
# Description: show help page information
#***************************************************************************
sub show_help_page {

	my ($self) = @_;

	# Initialise usage statement to print if usage is incorrect
	my ($HELP)  = "\n\t Usage: $0 -m=[option] -i=[control file]\n";

        $HELP  .= "\n\t ### Main functions"; 
		$HELP  .= "\n\t -m=1  Index target files under path '$ENV{DIGS_GENOMES}'"; 
		$HELP  .= "\n\t -m=2  Screen"; 
		$HELP  .= "\n\t -m=3  Reassign"; 
		$HELP  .= "\n\t -m=4  Defragment"; 
		$HELP  .= "\n\t -m=5  Consolidate"; 
		$HELP  .= "\n\t -m=6  Create standard locus IDs\n"; 
		$HELP  .= "\n\t Run  $0 -e to see information on utility functions \n\n"; 

	print $HELP;
}

############################################################################
# UTILITY FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  show_utility_help_page
# Description: show help page information for utility functions
#***************************************************************************
sub show_utility_help_page {

	my ($self) = @_;

	# Initialise usage statement to print if usage is incorrect
	my ($HELP)  = "\n\t Usage: $0 -m=[option] -i=[control file]\n";

        $HELP  .= "\n\t ### Utility functions"; 
		$HELP  .= "\n\t -u=1  Add extra tables to screening DB"; 
		$HELP  .= "\n\t -u=2  Flush screening DB"; 
		$HELP  .= "\n\t -u=3  Drop screening DB"; 
		$HELP  .= "\n\t -u=4  Show BLAST chains"; 
		$HELP  .= "\n\t -u=5  Show locus chains"; 
		$HELP  .= "\n\t -u=6  Summarise genomes (short, by species)";
		$HELP  .= "\n\t -u=7  Summarise genomes (long, by target file)";
		$HELP  .= "\n\t -u=8  Translate DB schema\n\n"; 
		#$HELP  .= "\n\t -u=6  Extract sequences using track";

	print $HELP;
}

#***************************************************************************
# Subroutine:  run_utility_process
# Description: handler for DIGS tool utility functions 
#***************************************************************************
sub run_utility_process {

	my ($self, $option, $infile) = @_;

 	# Show title
	$self->show_title();  

	# Initialise (an infile must be defined)
	unless ($infile) { die "\n\t Option '$option' requires an infile\n\n"; }
	$self->initialise($infile);
	
	# Load/create the screening
	$self->load_screening_db($infile);

	# Hand off to functions 
	if ($option eq 1) { # Add a table of data to the screening database
		$self->extend_screening_db();
	}
	elsif ($option eq 2) { # Flush screening DB
		my $db = $self->{db};
		$db->flush_screening_db();
	}
	elsif ($option eq 3) { # Drop screening DB 
		my $db = $self->{db};
		$db->drop_screening_db();    
	}
	elsif ($option eq 4) { # Show the BLAST chains for each extracted locus
		$self->show_blast_chains();
	}
	elsif ($option eq 5) { # Consolidate DIGS results into higher level loci 
		$self->show_locus_chains();
	}

	elsif ($option eq 6)    { # Summarise target genome directory (short)
		my $target_db_obj = TargetDB->new($self);
		$target_db_obj->summarise_genomes_short();
	}
	elsif ($option eq 7) { # Summarise target genome directory (long)
		my $target_db_obj = TargetDB->new($self);
		$target_db_obj->summarise_genomes_long();
	}
	elsif ($option eq 8) { # DB schema translation
		my $db_obj = $self->{db};
		$db_obj->translate_schema();
	}
	elsif ($option eq 9) {
		unless ($infile) {  die "\n\t Option '$option' requires an infile\n\n"; }
		my $loader_obj = ScreenBuilder->new($self);
		my @extracted;
		$loader_obj->extract_track_sequences(\@extracted, $infile);
	}
	else {
		print "\n\t  Unrecognized option '-u=$option'\n";
	}
}

############################################################################
# VALIDATE FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  validate
# Description: 
#***************************************************************************
sub validate {

	my ($self) = @_;

	# Build the test set

}

############################################################################
# EOF
############################################################################
