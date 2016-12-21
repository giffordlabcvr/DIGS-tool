#!usr/bin/perl -w
############################################################################
# Module:      DIGS.pm
# Description: Exploring genomes using BLAST and a relational database
# History:     December 2013: Created by Robert Gifford 
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
my $devtools   = DevTools->new();

# Flags
my $verbose      = undef;
1;

############################################################################
# LIFECYCLE & TOP LEVEL HANDLER FUNCTIONS
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
		blast_bin_path         => $parameter_ref->{blast_bin_path},
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

#***************************************************************************
# Subroutine:  run_digs_process
# Description: handler for main DIGS functions 
#***************************************************************************
sub run_digs_process {

	my ($self, $option, $ctl_file) = @_;

 	# Show title
	$self->show_title();  

	# Do initial set up and sanity checking for options that require it
	# An infile must be defined
	unless ($ctl_file) {
		unless ($option eq 8) {  die "\n\t Option '$option' requires an infile\n\n"; }
	}
	else {
		$self->initialise($ctl_file);
	}
	
	# If it exists, load the screening database specified in the control file
	if ( $option > 1 and $option < 8) {
		my $db_name = $self->{db_name};
		unless ($db_name) { die "\n\t Error: no DB name defined \n\n\n"; }
		my $db_obj = ScreeningDB->new($self);
		$db_obj->load_screening_db($db_name);	
		$self->{db} = $db_obj; # Store the database object reference 
	}

	# Hand off to functions 
	if ($option eq 1)    { # Create a screening DB 
		$self->create_screening_db($ctl_file);
	}
	elsif ($option eq 2) { # Screenne;
		$self->index_previously_executed_searches();
		$self->set_up_screen();
		$self->screen();	
	}
	elsif ($option eq 3) { # Reassign data in Exracted table
		my @extracted_seqs;
		$self->initialise_reassign(\@extracted_seqs); # Set up 
		$self->reassign(\@extracted_seqs);	
	}
	elsif ($option eq 4) { # Reassign data in Exracted table
		$self->interactive_defragment();	
	}
	elsif ($option eq 5) { # Flush screening DB
		my $db = $self->{db};
		$db->flush_screening_db();
	}
	elsif ($option eq 6) { # Drop screening DB 
		my $db = $self->{db};
		$db->drop_screening_db();    
	}
	elsif ($option eq 7) { # Add a table of data to the screening database
		$self->extend_screening_db();
	}
	elsif ($option eq 8) { # Add a table of data to the screening database
		my $target_db_obj = TargetDB->new($self);
		$target_db_obj->refresh_genomes();
	}
	else {
		print "\n\t  Unrecognized option '-m=$option'\n";
	}
}

#***************************************************************************
# Subroutine:  run_utility_process
# Description: handler for DIGS tool utility functions 
#***************************************************************************
sub run_utility_process {

	my ($self, $option, $infile) = @_;

 	# Show title
	$self->show_title();  

	# Hand off to functions 
	if ($option eq 1)    { # Summarise target genome directory (short)
		my $target_db_obj = TargetDB->new($self);
		$target_db_obj->summarise_genomes_short();
	}
	elsif ($option eq 2) { # Summarise target genome directory (long)
		my $target_db_obj = TargetDB->new($self);
		$target_db_obj->summarise_genomes_long();
	}
	elsif ($option eq 3) {
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
# INITIALISING FUNCTIONS
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
# Description: set up for reassigning the sequences in the Extracted_sequences table
#***************************************************************************
sub initialise_reassign {

	my ($self, $extracted_seqs_ref) = @_;

	# Create a unique ID and report directory for this run
	my $output_path = $self->{output_path};
	my $process_id  = $self->{process_id};
	my $db          = $self->{db};
	my $db_name     = $db->{db_name};
	unless ($db and $db_name and $process_id and $output_path) { die; }
	
	# Create report directory
	my $loader_obj = $self->{loader_obj};
	$loader_obj->create_output_directories($self);

	# Get the assigned data
	my $extracted_table = $db->{extracted_table};
	my @fields  = qw [ record_id probe_type assigned_name assigned_gene 
	                       extract_start extract_end sequence organism ];
	$extracted_table->select_rows(\@fields, $extracted_seqs_ref);

	# Set up the reference library
	$loader_obj->setup_reference_library($self);

}

#***************************************************************************
# Subroutine:  set_up_screen
# Description: set up files and directories for screening
#***************************************************************************
sub set_up_screen {

	my ($self) = @_;
	
	# Set up the screening queries
	my %queries;
	my $loader_obj = $self->{loader_obj};
	unless ($loader_obj) { die; }  # Sanity checking
	my $total_queries = $loader_obj->set_up_screen($self, \%queries);
	unless ($total_queries)  { 
		print "\n\t  Exiting without screening.\n\n";	
		exit;
	}
	$self->{queries}       = \%queries;
	$self->{total_queries} = $total_queries;
}


############################################################################
# MAIN FUNCTIONS - TOP LEVEL PROCESS FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  create_screening_db
# Description: create a screening datbase 
#***************************************************************************
sub create_screening_db {

	my ($self, $ctl_file) = @_;

	# Get parameters
	my $loader_obj = $self->{loader_obj};
	unless ($loader_obj) { die; }  # Sanity checking
	my $db_name = $loader_obj->{db_name};
	unless ($db_name)  { die; } 
		
	my $db_obj = ScreeningDB->new($loader_obj);
	$db_obj->create_screening_db($db_name);	
	$self->{db} = $db_obj; # Store the database object reference 
}

#***************************************************************************
# Subroutine:  screen
# Description: do the core database-integrated genome screening processes
#***************************************************************************
sub screen {

	my ($self, $mode) = @_;

	# Get relevant member variables and objects
	my $db_ref      = $self->{db};
	my $queries_ref = $self->{queries};
	my $loader_obj  = $self->{loader_obj};
	unless ($db_ref and $queries_ref and $loader_obj) { die; }  # Sanity checking

	# Iterate through and excute screens
	print "\n\t  Starting database-integrated genome screening\n";
	my @probes = keys %$queries_ref;

	my $num_queries = 0;
	my %crossmatching;
	$self->{crossmatching} = \%crossmatching;
	foreach my $probe_name (@probes) {
		
		# Get the array of queries for this target file
		my $probe_queries = $queries_ref->{$probe_name};
		foreach my $query_ref (@$probe_queries) {  
	
			# Increment query count
			$num_queries++;		
			$self->{num_queries} = $num_queries;

			# Do the 1st BLAST (probe vs target)
			$self->search($query_ref);
			
			# Update DB and extract sequences to assign
			my @new_hits;
			$self->update_db($query_ref, \@new_hits);

			# Extract newly identified or extended sequences
			my @extracted;
			$self->extract($query_ref,\@extracted);	
	
			# Do the 2nd BLAST (hits from 1st BLAST vs reference library)
			$self->assign($query_ref, \@extracted);
			die;	
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
# Description: reassign sequences in the extracted_table (for use after
#              the reference library has been updated)
#***************************************************************************
sub reassign {
	
	my ($self, $extracted_seqs_ref) = @_;

	# Get data structures and variables from self
	my $blast_obj       = $self->{blast_obj};
	my $result_path     = $self->{report_dir};
	my $db              = $self->{db};
	my $extracted_table = $db->{extracted_table};
	unless ($extracted_table) { die; }
	
	# Iterate through the matches
	print "\n\n\t  Reassigning Extracted_sequences table\n";
	my $count = 0;
	my %reassign_matrix;
	my %unique_keys;
	foreach my $hit_ref (@$extracted_seqs_ref) {

		# Set the linking to the BLAST result table
		my $blast_id  = $hit_ref->{record_id};
		$hit_ref->{blast_id} = $blast_id;
		delete $hit_ref->{record_id};
		my $extract_start   = $hit_ref->{extract_start};
		my $extract_end     = $hit_ref->{extract_end};
		$hit_ref->{subject_start} = $extract_start;
		$hit_ref->{subject_end}   = $extract_end;
		delete $hit_ref->{extract_start};
		delete $hit_ref->{extract_end};
	
		# Execute the 'reverse' BLAST (2nd BLAST in a round of paired BLAST)	
		my $previous_assign = $hit_ref->{assigned_name};
		my $previous_gene   = $hit_ref->{assigned_gene};
		#print "\n\t Redoing assign for record ID $blast_id assigned to $previous_assign";
		#print "\n\t coordinates: $extract_start-$extract_end";
		$self->do_reverse_blast($hit_ref);
		
		$count++;
		if (($count % 100) eq 0) {
			print "\n\t  Checked $count rows";
		}

		my $assigned_name = $hit_ref->{assigned_name};
		my $assigned_gene = $hit_ref->{assigned_gene};
		if ($assigned_name ne $previous_assign or  $assigned_gene ne $previous_gene) {
			
			# Show the results
			print "\n\t ##### Reassigned $blast_id from $previous_assign ($previous_gene)";
			print " to $assigned_name ($assigned_gene)";
			my $where = " WHERE Record_id = $blast_id ";
			
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
			$extracted_table->update($hit_ref, $where);
		}
	}
	
	# Write out the matrix
	$self->write_matrix(\%reassign_matrix, \%unique_keys);

	# Cleanup
	my $output_dir = $self->{report_dir};
	my $command1 = "rm -rf $output_dir";
	system $command1;
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


############################################################################
# MAIN SCREENING LOOP
############################################################################

#***************************************************************************
# Subroutine:  search
# Description: execute a search (i.e. run a BLAST query)
#***************************************************************************
sub search {
	
	my ($self, $query_ref) = @_;

	# Get relevant member variables and objects
	my $db_ref          = $self->{db};
	my $blast_obj       = $self->{blast_obj};
	my $tmp_path        = $self->{tmp_path};
	my $min_length      = $self->{seq_length_minimum};

	# Show screening process
	my $total_queries   = $self->{total_queries};
	my $num_queries     = $self->{num_queries};	
	unless ($num_queries and $total_queries) { die; }
	my $percent_prog    = ($num_queries / $total_queries) * 100;
	my $f_percent_prog  = sprintf("%.2f", $percent_prog);

	# Sanity checking
	unless ($blast_obj)       { die; } 
	unless ($tmp_path)        { die; } 
	unless ($db_ref)          { die; } 

	# Get query details
	my $probe_id     = $query_ref->{probe_id};
	my $probe_name   = $query_ref->{probe_name};
	my $probe_gene   = $query_ref->{probe_gene};
	my $probe_path   = $query_ref->{probe_path};
	my $organism     = $query_ref->{organism};
	my $version      = $query_ref->{version};
	my $data_type    = $query_ref->{data_type};
	my $target_name  = $query_ref->{target_name};
	my $target_path  = $query_ref->{target_path};
	my $blast_alg    = $query_ref->{blast_alg};
	my $cutoff       = $query_ref->{bitscore_cutoff};
	my $result_file  = $tmp_path . "/$probe_id" . "_$target_name.blast_result.tmp";
	unless ($cutoff) {die;}

	# Do the BLAST search
	print "\n\t  $blast_alg screen: $num_queries (%$f_percent_prog done): '$organism' file '$target_name' with probe $probe_id";   
	$blast_obj->blast($blast_alg, $target_path, $probe_path, $result_file);
	
	# Parse out the alignment
	my @hits;
	$blast_obj->parse_tab_format_results($result_file, \@hits, $cutoff);
	# TODO: catch error from BLAST and don't update "Searches_performed" table	

	# Clean up - remove the result file
	my $rm_command = "rm $result_file";
	system $rm_command;

	# Store new new BLAST hits that meet conditions
	my $blast_results_table = $db_ref->{blast_results_table};
	my $num_hits = scalar @hits;
	if ($num_hits > 0) {
		print "\n\t\t # $num_hits matches to probe: $probe_name, $probe_gene";
	}
	foreach my $hit_ref (@hits) {
		
		#$devtools->print_hash($hit_ref ); die;
		
		if ($min_length) { # Skip sequences that are too short
			my $start  = $hit_ref->{aln_start};
			my $end    = $hit_ref->{aln_stop};
			if ($end - $start < $min_length) {  
				next;
			}
		}
	}
	
	foreach my $hit_ref (@hits) {
		
		if ($min_length) { # Skip sequences that are too short
			my $start  = $hit_ref->{aln_start};
			my $end    = $hit_ref->{aln_stop};
			if ($end - $start < $min_length) {  next; }
		}
			
		# Record hit in BLAST results table
		$hit_ref->{organism}      = $organism;
		$hit_ref->{version}       = $version;
		$hit_ref->{data_type}     = $data_type;
		$hit_ref->{target_name}   = $target_name;
		$hit_ref->{probe_id}      = $probe_id;
		$hit_ref->{probe_name}    = $probe_name;
		$hit_ref->{probe_gene}    = $probe_gene;
		$hit_ref->{probe_type}    = $query_ref->{probe_type};
		$hit_ref->{hit_length}    = $hit_ref->{align_len};
		$hit_ref->{subject_start} = $hit_ref->{aln_start}; 	# Rename coordinate fields to match DB
		$hit_ref->{subject_end}   = $hit_ref->{aln_stop}; 	# Rename coordinate fields to match DB
		$hit_ref->{query_end}     = $hit_ref->{query_stop}; # Rename coordinate fields to match DB

		$blast_results_table->insert_row($hit_ref);
	} 

	# Update the searches table
	my $searches_table = $db_ref->{searches_table};
	$searches_table->insert_row($query_ref);
}

#***************************************************************************
# Subroutine:  update_db
# Description: 
#***************************************************************************
sub update_db {
	
	my ($self, $query_ref, $extracted_ref) = @_;

	# Get the information for this query
	my $target_name = $query_ref->{target_name};
	my $probe_name  = $query_ref->{probe_name};
	my $probe_gene  = $query_ref->{probe_gene};

	# Get data structures and variables from self
	my $db = $self->{db};
	my $blast_results_table = $db->{blast_results_table};

	# Index BLAST chains
	my %extracted;
	my $where = " WHERE target_name = '$target_name' ";
	$self->index_blast_chains(\%extracted, $where);

	# Get all BLAST results for this target file
	my @hits;
	my $where = " WHERE Target_name = '$target_name'
                  AND probe_name = '$probe_name' 
                  AND probe_gene = '$probe_gene'
	              ORDER BY scaffold, subject_start";
	my @fields = qw [ record_id 
	               organism data_type version 
                   probe_name probe_gene probe_type
				   orientation scaffold target_name
                   subject_start subject_end 
		           query_start query_end ];
	$blast_results_table->select_rows(\@fields, $hits_ref, $where); 
	my $num_hits = scalar @hits;
	if ($verbose) { print "\n\t ### There are $num_hits hits to extract"; }

	# Merge overlapping, redundant and fragmented BLAST hits
	# Get the set of hits (rows in BLAST_results table) to look at
	my @hits;

	my ($self, $hits_ref, $target_name, $probe_name, $probe_gene) = @_;
	
	# Get relevant member variables and objects

	my $redundancy_mode     = $self->{redundancy_mode};

	unless ($redundancy_mode) { die; } 
	
	# Set the fields to get values for
	my @fields = qw [ record_id scaffold orientation
	                  subject_start subject_end
                      query_start query_end ];

	# Build an SQL "where" statement to control what hits are selected
	my $where  = " WHERE target_name = '$target_name'";

	# Handle different modes
	if ($redundancy_mode eq 1) {
		return; # Do nothing else
	}
	if ($self->{redundancy_mode} eq 2) {
		# In mode 2, we select all BLAST table rows with the same value for 'probe_gene'
		$where .= " AND probe_gene = '$probe_gene' ";
	}
	elsif ($self->{redundancy_mode} eq 3) {
		# In mode 3, we select all BLAST table rows with the same value for both 'probe_gene' and 'probe_gene'
		$where .= " AND probe_name = '$probe_name'
                    AND probe_gene = '$probe_gene' ";
	}

	# Order by ascending start coordinates within each scaffold in the target file
	$where .= "ORDER BY scaffold, subject_start ";
	
	# Get the table rows for these hits
	#$blast_results_table->select_rows(\@fields, $hits_ref, $where);



	
	# Apply consolidation rules to overlapping and/or redundant hits
	my %merged;
	my %retained;
	$self->do_consolidation($query_ref, \@hits, \%merged, \%retained);



	# Update database
	$self->update_db_loci($table, \@hits, \%merged, \%retained);

		# Skip previously extracted hits
		#my $record_id   = $hit_ref->{record_id};
		#if ($verbose) {
		#	print "\n\t Checking $record_id in extracted table";
		#}
		#if ($extracted{$record_id}) { 
		#	if ($verbose) {
		#		print "\t ALREADY EXTRACTED";
		#	}
		#	next;
		#}	
}

#***************************************************************************
# Subroutine:  extract_sequences
# Description: extract sequences from target databases
#***************************************************************************
sub extract {

	my ($self, $query_ref, $hits_ref, $extracted_ref) = @_;

	# Get paths, objects, data structures and variables from self
	my $blast_obj   = $self->{blast_obj};
	my $buffer      = $self->{extract_buffer};
	my $target_path = $query_ref->{target_path};

	# Iterate through the list of sequences to extract
	my $new_hits = scalar @$hits_ref;
	foreach my $hit_ref (@$hits_ref) {
				
		# Add any buffer 
		my $orientation   = $hit_ref->{orientation};
		if ($buffer) {
			if ($orientation eq '-') {
				$hit_ref->{subject_start} = $hit_ref->{subject_start} + $buffer;
				$hit_ref->{subject_end}   = $hit_ref->{subject_end} - $buffer;
				if ($hit_ref->{subject_end} < 1) { # Don't allow negative coordinates
					$hit_ref->{subject_end} = 1;
				}	
			}
			else {
				$hit_ref->{subject_start} = $hit_ref->{subject_start} - $buffer;
				if ($hit_ref->{subject_start} < 1) { # Don't allow negative coordinates
					$hit_ref->{subject_start} = 1;
				}	
				$hit_ref->{subject_end}   = $hit_ref->{subject_end} + $buffer;
			}
		}
		
		# Extract the sequence
		my $sequence = $blast_obj->extract_sequence($target_path, $hit_ref, $buffer);
		if ($sequence) {	
			if ($verbose) {
				print "\t ........extracting";
			}
			my $seq_length = length $sequence; # Set sequence length
			$hit_ref->{sequence_length} = $seq_length;
			$hit_ref->{sequence} = $sequence;
			push (@$extracted_ref, $hit_ref);
		}
		else {
			if ($verbose) {
				print "\n\t Sequence extraction failed";
			}
		}
	}
	#my $num_extracted = scalar @$extracted_ref;	
	#print "\n\t Num extracted $num_extracted";	

	return $new_hits;
}

#***************************************************************************
# Subroutine:  assign
# Description: assign sequences that matched probes in a BLAST search 
#***************************************************************************
sub assign {
	
	my ($self, $query_ref, $extracted_ref) = @_;
	
	# Get parameters from self
	my $db_ref = $self->{db};
	my $extracted_table    = $db_ref->{extracted_table}; 
	my $blast_chains_table = $db_ref->{blast_chains_table}; 
	
	# Iterate through the matches
	my $assigned_count   = 0;
	my $crossmatch_count = 0;
	my $new_sequences    = 0;
	foreach my $hit_ref (@$extracted_ref) {

		# Set the linking to the BLAST result table
		my $blast_id  = $hit_ref->{record_id};
		$hit_ref->{blast_id} = $blast_id;
		
		# Execute the 'reverse' BLAST (2nd BLAST in a round of paired BLAST)	
		my %data = %$hit_ref; # Make a copy
		my $assigned = $self->do_reverse_blast(\%data);
		if ($assigned) { $assigned_count++; }
		my $probe_name  = $query_ref->{probe_name};
		my $probe_gene  = $query_ref->{probe_gene};
		my $probe_key = $probe_name . '_' . $probe_gene; 
		
		# Record cross-matching
		if ($probe_key ne $assigned) {
			$crossmatch_count++;
			$self->update_cross_matching($probe_key, $assigned);
		}

		# Insert the data to the Extracted_sequences table
		my $extract_id = $extracted_table->insert_row(\%data);
		#$devtools->print_hash(\%data); die;

		# Insert the data to the BLAST_chains table
		$data{extract_id} = $extract_id;
		$blast_chains_table->insert_row(\%data);

	}
	print "\n\t\t # $new_sequences newly identified hits";
	
	if ($new_sequences > 0) {
		print "\n\t\t # $assigned_count extracted sequences matched to reference library";
		print "\n\t\t # $crossmatch_count cross-matched to something other than the probe";
	}
}

#***************************************************************************
# Subroutine:  do_reverse_blast
# Description: Execute the 2nd BLAST in a round of paired BLAST
#***************************************************************************
sub do_reverse_blast {

	my ($self, $hit_ref) = @_;
	
	# Get paths and objects from self
	my $result_path   = $self->{tmp_path};
	my $blast_obj     = $self->{blast_obj};
	unless ($result_path and $blast_obj) { die; }
	
	if ($hit_ref->{subject_start} and $hit_ref->{subject_end}) { 
		# If we are coming direct from the 1st BLAST, do this translation
		$hit_ref->{extract_start} = $hit_ref->{subject_start};
		$hit_ref->{extract_end}   = $hit_ref->{subject_end};
	}

	# Get required data about the hit, prior to performing reverse BLAST
	my $blast_id      = $hit_ref->{blast_id};
	my $sequence      = $hit_ref->{sequence};
	my $organism      = $hit_ref->{organism};
	my $probe_type    = $hit_ref->{probe_type};
	
	# Sanity checking
	unless ($probe_type and  $blast_id and $organism) { die; }
	unless ($sequence) {  die "\n\t # ERROR: No sequence found for reverse BLAST"; } 
	
	# Make a FASTA query file for the reverse BLAST procedure
	$sequence =~ s/-//g;   # Remove any gaps that might happen to be there
	$sequence =~ s/~//g;   # Remove any gaps that might happen to be there
	$sequence =~ s/\s+//g; # Remove any gaps that might happen to be there
	my $fasta      = ">$blast_id\n$sequence";
	my $query_file = $result_path . $blast_id . '.fas';
	$fileio->write_text_to_file($query_file, $fasta);
	my $result_file = $result_path . $blast_id . '.blast_result';
		
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
		$hit_ref->{bit_score}        = 0;
		$hit_ref->{e_value_exp}      = 0;
		$hit_ref->{e_value_num}      = 0;
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
		$hit_ref->{bit_score}        = $top_match->{bit_score};
		$hit_ref->{e_value_exp}      = $top_match->{e_value_exp};
		$hit_ref->{e_value_num}      = $top_match->{e_value_num};
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


############################################################################
# DEFRAGMENTING - comparing and merging overlapping/adjacent hits
###########################################################################	

#***************************************************************************
# Subroutine:  do_consolidation
# Description: apply consolidation rules to overlapping/redundant hits
#***************************************************************************
sub do_consolidation {
	
	my ($self, $query_ref, $hits_ref, $retained_ref, $merged_ref) = @_;
	
	# Iterate through consolidating as we go
	my $i;
	my %last_hit;
	my $merged_count = 0;
	foreach my $hit_ref (@$hits_ref)  {

		# Get current hit values
		$i++;
		my $record_id          = $hit_ref->{record_id};
		my $scaffold           = $hit_ref->{scaffold};
		my $orientation        = $hit_ref->{orientation};
		my $subject_start      = $hit_ref->{subject_start};
		
		# Get last hit values
		my $last_record_id     = $last_hit{record_id};
		my $last_scaffold      = $last_hit{scaffold};
		my $last_orientation   = $last_hit{orientation};
		my $merged;

		if ($verbose) {
			print "\n\t $subject_start: RECORD ID $record_id";
		}

		# Keep if first in this loop process
		if ($i eq 1) {
			$retained_ref->{$record_id} = 1;
		}
		# ...or if first hit on scaffold
		elsif ($scaffold ne $last_scaffold) {
			$retained_ref->{$record_id} = 1;
		}
		# ...or if in opposite orientation to last hit
		elsif ($orientation ne $last_orientation) {
			$retained_ref->{$record_id} = 1;
		}
		else { # If we get this far we have two hits on the same scaffold
			
			# Check whether to consolidate hit on the same target scaffold
			$merged = $self->inspect_adjacent_hits(\%last_hit, $hit_ref);

			# Keep track of the outcome
			if ($merged) {
				$merged_count++;
				my %hit = %last_hit; # Make a copy
				$hit{hit_length} = ($hit{subject_end} - $hit{subject_start}) + 1;
				$merged_ref->{$record_id} = \%hit;
				$retained_ref->{$last_record_id} = 1;
			}
			else { $retained_ref->{$record_id} = 1; }
		}
		
		# Update the 'last hit' to the current one before exiting this iteration
		unless ($merged) {
			$last_hit{record_id}     = $record_id;
			$last_hit{scaffold}      = $scaffold;
			$last_hit{orientation}   = $orientation;
			$last_hit{subject_start} = $hit_ref->{subject_start};
			$last_hit{subject_end}   = $hit_ref->{subject_end};
			$last_hit{query_start}   = $hit_ref->{query_start};
			$last_hit{query_end}     = $hit_ref->{query_end};
		}
	}
}

#***************************************************************************
# Subroutine:  inspect_adjacent_hits
# Description: determine whether two adjacent hits should be joined
#***************************************************************************
sub inspect_adjacent_hits {
	
	my ($self, $last_hit_ref, $hit_ref) = @_;

	# Get parameters for consolidating hits from self
	my $defragment_range   = $self->{defragment_range};
	unless ($defragment_range) { die "\n\t Defragment range is not set\n\n\n"; } 

	# Get data for 1st hit (i.e. most 'leftward' in the target sequence)
	my $last_record_id     = $last_hit_ref->{record_id};
	my $last_subject_start = $last_hit_ref->{subject_start};
	my $last_subject_end   = $last_hit_ref->{subject_end};
	my $last_query_start   = $last_hit_ref->{query_start};
	my $last_query_end     = $last_hit_ref->{query_end};
	
	# Get data for 2nd hit (i.e. most 'rightward' in the target sequence)
	my $record_id          = $hit_ref->{record_id};
	my $scaffold           = $hit_ref->{scaffold};
	my $orientation        = $hit_ref->{orientation};
	my $subject_start      = $hit_ref->{subject_start};
	my $subject_end        = $hit_ref->{subject_end};
	my $query_start        = $hit_ref->{query_start};
	my $query_end          = $hit_ref->{query_end};

	# Calculate the gap between these sequence hits
	my $subject_gap = $subject_start - $last_subject_end;

	print "\n\n\t #### Checking whether to consolidate $last_record_id and $record_id on $scaffold";
	print "\n\t #### Q Last:  $last_query_start\t $last_query_end";
	print "\n\t #### S Last:  $last_subject_start\t $last_subject_end";
	print "\n\t #### Q This:  $query_start\t $query_end";
	print "\n\t #### S This:  $subject_start\t $subject_end";
	print "\t #### Gap:   $subject_gap\n";

	# Calculate the gap between the query coordinates of the two matches
	# Note - this may be a negative number if the queries overlap 
	my $query_gap;	
	if ($orientation eq '+') {
		$query_gap = $query_start - $last_query_end;
	}
	elsif ($orientation eq '-') {
		$query_gap = $last_query_start - $query_end;
	}	
	else { print "\n\t orientation = $orientation\n\n";  die; }
	
	# Deal with contingencies that mean hits definitely should or should not be merged
	# 1. Hit is entirely within a previous hit it is redundant
	if ($last_subject_start <= $subject_start and $last_subject_end >= $subject_end) {
		return 1; # Effectively discarding current hit
	}

	### Deal with situations where hits are partially (but not completely) overlapping
	#   or where they are close enough to each other consider merging them
	
	if ($subject_gap < 1) {
		# Describe the merging event:
		$last_hit_ref->{subject_end} = $hit_ref->{subject_end};
		$last_hit_ref->{query_end}   = $hit_ref->{query_end};
		return 1;  # Hits have been merged
	}
	
	# For positive orientation hits
	# If the query_start of this hit is upstream of the query_end of the last, then its a distinct hit
	elsif ($orientation eq '+' and $query_start < $last_query_end) {
		return 0;
	}
	# For negative orientation hits
	# If the query_start of this hit < query_end of the last, then its a distinct hit
	elsif ($orientation eq '-' and ($last_query_start - $query_end) ) {
		return 0;
	}

	# If not, then check the intervening distance between the hits
	else {

		my $nt_query_gap = $query_gap * 3;
		if (($subject_gap - $nt_query_gap) < $defragment_range) {
			$last_hit_ref->{subject_end} = $hit_ref->{subject_end};
			$last_hit_ref->{query_end}   = $hit_ref->{query_end};
			return 1;  # Consolidated (joined) hits
		}
		else {
			return 0;  # Distinct hit
		}
	}
}

#***************************************************************************
# Subroutine:  update_db_loci  
# Description: update the BLAST results table with overlapping hit info
#***************************************************************************
sub update_db_loci {
	
	my ($self, $hits_ref, $retained_ref, $merged_ref) = @_;

	#$devtools->print_hash($merged_ref); exit;

	# Get relevant member variabless
	my $db_ref  = $self->{db};
	my $blast_results_table = $db_ref->{blast_results_table};
	my $extracted_table     = $db_ref->{extracted_table};

	# Update BLAST table with merged hits
	my $blast_updated = 0;
	my $blast_deleted = 0;
	my $extract_deleted = 0;
	my @ids = keys %$merged_ref;
	foreach my $record_id (@ids)  {

		my $where   = " WHERE record_id = $record_id ";
		my $hit_ref = $merged_ref->{$record_id};		
		delete $hit_ref->{record_id};
		#print "\n\t ## updating hit record ID $record_id in BLAST results";
		my $this_id = $blast_results_table->update($hit_ref, $where);
		$blast_updated++;
	
		# Update the 'BLAST_chains' table accordingly
		#$hit{this_query_id} = $this_id; 
	
		# Delete from Extracted_sequences (because we need to re-extract and assign)
		#print "\n\t ## deleting hit BLAST ID $record_id in Extracted_sequences";
		my $ex_where = " WHERE blast_id = $record_id ";
		$extracted_table->delete_rows($ex_where);
		$extract_deleted++;
	}

	# Delete redundant rows from BLAST results
	foreach my $hit_ref (@$hits_ref)  {
		my $record_id = $hit_ref->{record_id};
		unless ($retained_ref->{$record_id}) {
			# Delete rows
			#print "\n\t ## deleting hit record ID $record_id in BLAST results";
			my $where = " WHERE record_id = $record_id ";
			$blast_results_table->delete_rows($where);
			$blast_deleted++;

			#print "\n\t ## deleting hit BLAST ID $record_id in Extracted_sequences";
			my $ex_where = " WHERE blast_id = $record_id ";
			$extracted_table->delete_rows($ex_where);
			$extract_deleted++;
		}
		else { 
			#print "\n\t ## keeping hit $record_id";
		}
	}
}

#***************************************************************************
# Subroutine:  interactive_defragment 
# Description: console driven menu options for interactive defragment
#***************************************************************************
sub interactive_defragment {

	my ($self) = @_;

	my $db = $self->{db};

	# Display current settings
	my $redundancy_mode= $self->{redundancy_mode};
	my $defragment_range=$self->{defragment_range};
	print "\n\n\t\t Current settings (based on control file)";
	print "\n\t\t redundancy mode: $redundancy_mode";
	print "\n\t\t defragment_range: $defragment_range";

	# Get the ordered hits from the Extracted table
	
	my $extracted_table = $db->{extracted_table};
	
	# Set the fields to get values for
	my @fields = qw [ record_id organism assigned_name assigned_gene
	                  scaffold orientation
	                  subject_start subject_end
                      extract_start extract_end sequence_length ];

	# Order by ascending start coordinates within each scaffold in the target file
	my $where = "ORDER BY organism, target_name, scaffold, extract_start ";
	
	# Get the relevant loci
	my @hits;
	$extracted_table->select_rows(\@fields, \@hits, $where);
	my $numhits = scalar @hits;
	print "\n\n\t # There are currently $numhits distinct records in the Extracted_sequences table";


	my %defragmented;
	my $choice;
	do {

		my $max = 100000;
		my $question1 = "\n\n\t # Set the range for merging hits";
		my $t_range = $console->ask_int_with_bounds_question($question1, $defragment_range, $max);		
		my $t_mode = 2;		
		$self->preview_defragment(\%defragmented, \@hits, $t_range, $t_mode);

		my @cluster_ids = keys %defragmented;
		my $cluster_count;
		foreach my $id (@cluster_ids) {
		
			$cluster_count++;
			my $hits_ref = $defragmented{$id};
			my $cluster_size = scalar @$hits_ref;

			if ($cluster_size > 1) {
				print "\n";
				foreach my $hit_ref (@$hits_ref) {

					my $organism      = $hit_ref->{organism};			
					my $assigned_name = $hit_ref->{assigned_name};			
					my $assigned_gene = $hit_ref->{assigned_gene};			
					my $scaffold      = $hit_ref->{scaffold};			
					my $extract_start = $hit_ref->{extract_start};			
					my $extract_end   = $hit_ref->{extract_end};
					print "\n\t\t CLUSTER $cluster_count $organism: $assigned_name: $assigned_gene: $scaffold $extract_start-$extract_end";

				}			
			}
		}
		
		# Prompt for what to do next
		print "\n\n\t\t Option 1: preview new parameters";
		print "\n\t\t Option 2: apply these parameters";
		print "\n\t\t Option 3: exit";
		my $list_question = "\n\n\t # Choose an option:";
		$choice = $console->ask_list_question($list_question, 3);

	} until ($choice > 1);

	if ($choice eq 2) {
		$self->defragment(\%defragmented, 1);
	}
	elsif ($choice eq 3) {
		exit;
	}
}

#***************************************************************************
# Subroutine:  preview_defragment 
# Description:
#***************************************************************************
sub preview_defragment {

	my ($self, $defragmented_ref, $hits_ref, $range, $t_mode) = @_;
	
	# Iterate through consolidating as we go
	my $i = 0;
	my $j = 1;
	my %last_hit;
	my %name_counts;
	my $initialised = undef;
	my $extracted_count = 0;
	foreach my $hit_ref (@$hits_ref)  {

		# Get hit values
		my $record_id     = $hit_ref->{record_id};
		my $assigned_name = $hit_ref->{assigned_name};
		my $scaffold      = $hit_ref->{scaffold};
		my $orientation   = $hit_ref->{orientation};
		my $extract_start = $hit_ref->{extract_start};
		my $extract_end   = $hit_ref->{extract_end};
		my $query_start   = $hit_ref->{query_start};
		my $query_end     = $hit_ref->{query_end};
		
		# Get last hit values
		my $last_record_id     = $last_hit{record_id};
		my $last_scaffold      = $last_hit{scaffold};
		my $last_assigned_name = $last_hit{assigned_name};
		my $last_orientation   = $last_hit{orientation};
		my $last_extract_start = $last_hit{extract_start};
		my $last_extract_end   = $last_hit{extract_end};		
		my $gap;

		$extracted_count++;
		#print "\n\t Extracted_sequences count $extracted_count";
	    if ($initialised) {

            my $new = $self->compare_adjacent_hits($hit_ref, \%last_hit, $range);
            if ($new) {
                
                # Finish record
				#$devtools->print_hash(\%defragmented); die;
				#$devtools->print_array($array_ref); die;
				#unless ($array_ref) { die; }
				my $array_ref = $defragmented_ref->{$j};
				$self->finish_cluster($array_ref, \%name_counts, $j);

                # Increment the count
                $j++;
            
                # Initialise new Missillac record
                $self->initialise_cluster($defragmented_ref, $hit_ref, $j);
            }
            else {
                
                # Extend current Missillac record
                $self->extend_cluster($defragmented_ref, $hit_ref, $j);

            }
        }
		else {
			#print "\n\t\t First record ($j)";
            $initialised = 'true';
            
            # Initialise new Missillac record
            $self->initialise_cluster($defragmented_ref, $hit_ref, $j);
			#$devtools->print_hash($defragmented_ref); die;
		}

		# Update last hit data
		$last_hit{record_id}     = $record_id;
		$last_hit{assigned_name} = $assigned_name;
		$last_hit{scaffold}      = $scaffold;
		$last_hit{orientation}   = $orientation;
		$last_hit{extract_start} = $extract_start;
		$last_hit{extract_end}   = $extract_end;
	}
}

#***************************************************************************
# Subroutine:  defragment 
# Description: function that implements a defragmentation process
#***************************************************************************
sub defragment {

	my ($self, $defragmented_ref, $mode) = @_;

	my $db = $self->{db};
	my $extracted_table = $db->{extracted_table};

	my @cluster_ids = keys %$defragmented_ref;
	#$devtools->print_hash(\%defragmented);
	foreach my $id (@cluster_ids) {
		my $hits_ref = $defragmented_ref->{$id};
		#$devtools->print_array($hits_ref);
		#print "\n\t # CLUSTER $id";
		my @delete_ids;
		my $highest = undef;
		my $lowest  = undef;
		my $longest  = undef;
		my $longest_id  = undef;
		my %genes;
		foreach my $hit_ref (@$hits_ref) {
			
			my $record_id     = $hit_ref->{record_id};			
			my $assigned_name = $hit_ref->{assigned_name};			
			my $assigned_gene = $hit_ref->{assigned_gene};			
			my $scaffold      = $hit_ref->{scaffold};			
			my $extract_start = $hit_ref->{extract_start};			
			my $extract_end   = $hit_ref->{extract_end};
			my $orientation   = $hit_ref->{orientation};
			$genes{$assigned_gene} = $assigned_gene;
			
			#print "\n\t\t RECORD ID:\t $record_id, $assigned_name: $scaffold $extract_start-$extract_end";
			my $length = $extract_end - $extract_start;
			if ($longest) {		
				if ($length > $longest) {
					$longest_id = $record_id;
				}
			}
			else {
				$longest = $length;
				$longest_id = $record_id;
			}							
		}
		
		# DELETE all but longest hit in cluster
		foreach my $hit_ref (@$hits_ref) {
			
			my $record_id     = $hit_ref->{record_id};			
			unless ($record_id eq $longest_id) {
				#print "\n\t\t DELETING RECORD ID:\t $record_id";
				my $where = " WHERE Record_ID = $record_id ";
				$extracted_table->delete_rows($where);
		
			}	
		}
		
		# UPDATE gene information if we have hits from > 1 different gene
		my @genes = keys %genes;
		my $num_genes = scalar @genes;
		if ($num_genes > 1) {
			my $joined_genes = join ('-', @genes);
			my $set = " WHERE Record_ID = $longest_id ";
			my %set;
			$set{assigned_gene} = $joined_genes;
			$extracted_table->update(\%set, $set);
		}
	}
}

#***************************************************************************
# Subroutine:  compare_adjacent_hits
# Description: compare two hits 
#***************************************************************************
sub compare_adjacent_hits {

	my ($self, $hit_ref, $last_hit_ref, $range) = @_;

	# Get the current hit values
	my $name           = $hit_ref->{assigned_name};
	my $scaffold       = $hit_ref->{scaffold};	
	my $extract_start  = $hit_ref->{extract_start};
	my $extract_end    = $hit_ref->{extract_end};
	my $orientation    = $hit_ref->{orientation};			

	# Get the last hit values
	my $last_name        = $last_hit_ref->{assigned_name};
	my $last_scaffold    = $last_hit_ref->{scaffold};	
	my $last_start       = $last_hit_ref->{extract_start};
	my $last_end         = $last_hit_ref->{extract_end};
	my $last_orientation = $last_hit_ref->{orientation};			
	
	if ($scaffold ne $last_scaffold) {
		return 1;
	}
	
	if ($orientation ne $last_orientation) {
		return 1;
	}
	
	my $gap;
	$gap = $extract_start - $last_end;		
	#print "\n\t #\t CALC: '$scaffold': '$extract_start'-'$last_end' = $gap";

	# Test whether to combine this pair into a set
	if ($gap < $range) {  # Combine
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
	my $name        = $hit_ref->{assigned_name};
	my $scaffold    = $hit_ref->{scaffold};	
	#print "\n\t New record ($count) [$name]";
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
	my $name        = $hit_ref->{assigned_name};
	my $scaffold    = $hit_ref->{scaffold};	
    
	#print "\n\t Extending record ($count) [$name]";
    my $array_ref = $defragmented_ref->{$count};
    push (@$array_ref, $hit_ref);

}

#***************************************************************************
# Subroutine:  finish_cluster
# Description:  
#***************************************************************************
sub finish_cluster {

	my ($self, $array_ref, $name_counts_ref, $j) = @_;

	unless ($array_ref) { die; }
	
	# Create summary data for this annotation cluster
	my %cluster_data;
    foreach my $hit_ref (@$array_ref) {
   	
    	# Get the current hit values
		my $name        = $hit_ref->{assigned_name};
		my $scaffold    = $hit_ref->{scaffold};	
		my $start       = $hit_ref->{extract_start};
		my $end         = $hit_ref->{extract_end};
		my $orientation = $hit_ref->{orientation};			
		#print "\n\t\t  # $name: $scaffold:  $start $end";	

    }
}


############################################################################
# FUNCTIONS FOR RECORDING CROSS-MATCHING
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
# Subroutine:  write matrix 
# Description: write cross-matching results as a matrix
#***************************************************************************
sub write_matrix {

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
# INDEXING DB TABLES
############################################################################

#***************************************************************************
# Subroutine:  index BLAST results by record id
# Description: Index loci in BLAST_results table by the 'record_id' field
#***************************************************************************
sub index_blast_chains {
	
	my ($self, $data_ref, $where) = @_;

	# Get relevant variables and objects
	my $db = $self->{db};
	my $blast_chains_table = $db->{blast_chains_table}; 
	my @fields = qw [ blast_id extract_id ];
	my @blast_chain_rows;
	$blast_chains_table->select_rows(\@fields, \@blast_chain_rows, $where);	
	foreach my $hit_ref (@blast_chain_rows) {
		my $blast_id   = $hit_ref->{blast_id};
		my $extract_id = $hit_ref->{extract_id};
		if ($data_ref->{$blast_id}) { die; } # BLAST ID should be unique
		$data_ref->{$blast_id} = $extract_id;	
	}
}

#***************************************************************************
# Subroutine:  index_previously_executed_searches 
# Description: index BLAST searches that have previously been executed
#***************************************************************************
sub index_previously_executed_searches {
	
	my ($self) = @_;

	my $db = $self->{db};	
	my $searches_table = $db->{searches_table};
	unless ($searches_table) { die "\n\t Searches_performed table not loaded\n\n"; }
	my @data;
	my @fields = qw [ record_id organism data_type version target_name
                      probe_name probe_gene ];
	my $where = " ORDER BY Record_ID ";
	$searches_table->select_rows(\@fields, \@data, $where);
	
	# Index the executed searches
	my %done;
	foreach my $data_ref (@data) {
		
		# Get the query parameters
		my $organism    = $data_ref->{organism};
		my $data_type   = $data_ref->{data_type};
		my $version     = $data_ref->{version};
		my $target_name = $data_ref->{target_name};
		my $probe_name  = $data_ref->{probe_name};
		my $probe_gene  = $data_ref->{probe_gene};
	
		# Sanity checking
		unless ( $organism and $data_type and $version and $target_name 
             and $probe_name and $probe_gene) { 
			die;
		};
		
		# Create the unique key for this search
		my @genome = ( $organism , $data_type, $version );
		my $genome_id = join ('|', @genome);
		my $probe_id  = $probe_name . '_' .  $probe_gene;
		my @key = ( $genome_id, $target_name, $probe_id );
		my $key = join ('|', @key);

		# Store the query in a hash indexed by it's unique key
		$done{$key} = $data_ref;		
	}
	$self->{previously_executed_searches} = \%done;
}

############################################################################
# SHOWING PROGRAM TITLE INFORMATION & HELP MENU
###########################################################################

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
        $HELP  .= "\n\t -m=1  Create screening DB"; 
		$HELP  .= "\n\t -m=2  Screen"; 
		$HELP  .= "\n\t -m=3  Reassign"; 
		$HELP  .= "\n\t -m=4  Defragment"; 
		$HELP  .= "\n\t -m=5  Flush screening DB"; 
		$HELP  .= "\n\t -m=6  Drop screening DB"; 
		$HELP  .= "\n\t -m=7  Manage ancillary tables"; 
		$HELP  .= "\n\t -m=8  Format genome directory ($ENV{DIGS_GENOMES})\n"; 
        $HELP  .= "\n\t ### Utility functions"; 
		$HELP  .= "\n\t -u=1  Summarise genomes (short, by species)";
		$HELP  .= "\n\t -u=2  Summarise genomes (long, by target file)";
		$HELP  .= "\n\t -u=3  Extract sequences using track\n\n";
	print $HELP;
}


############################################################################
# EOF
############################################################################