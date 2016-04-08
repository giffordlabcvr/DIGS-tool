#!/usr/bin/perl -w
############################################################################
# Module:      DIGS.pm
# Description: Genome screening pipeline using reciprocal BLAST
# History:     December 2009: Created by Robert Gifford 
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

# Program components
use DIGS::ScreenBuilder; # Functions to set up screen

############################################################################
# Globals
############################################################################

# Base objects
my $fileio    = FileIO->new();
my $console   = Console->new();
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

############################################################################
# TOP LEVEL HANDLER
############################################################################

#***************************************************************************
# Subroutine:  run_digs_process
# Description: handler for main DIGS functions 
#***************************************************************************
sub run_digs_process {

	my ($self, $option, $ctl_file) = @_;

 	# Show title
	$self->show_title();  

	# Do initial set up and sanity checking for options that require it
	unless ($option > 6) { # Validate configurations that require a control file
		# An infile must be defined
		unless ($ctl_file) {  die "\n\t Option '$option' requires an infile\n\n"; }
		$self->initialise($ctl_file);
	}

	# If it exists, load the screening database specified in the control file
	if ( $option > 1 and $option < 7) {
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
	elsif ($option eq 2) { # Screen
		$self->screen();	
	}
	elsif ($option eq 3) { # Reassign data in Exracted table
		$self->reassign();	
	}
	elsif ($option eq 4) { # Flush screening DB
		my $db = $self->{db};
		$db->flush_screening_db();
	}
	elsif ($option eq 5) { # Drop screening DB 
		my $db = $self->{db};
		$db->drop_screening_db();    
	}
	elsif ($option eq 6) { # Add a table of data to the screening database
		$self->extend_screening_db();
	}
	elsif ($option eq 7) { # Add a table of data to the screening database
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

	my ($self, $option, $ctl_file) = @_;

 	# Show title
	$self->show_title();  

	# Do initial set up and sanity checking for options that require it
	unless ($option < 6) { # Validate configurations that require a control file
		# An infile must be defined
		unless ($ctl_file) {  die "\n\t Option '$option' requires an infile\n\n"; }
		$self->initialise($ctl_file);
	}

	# If it exists, load the screening database specified in the control file
	if ( $option > 2 and $option < 7) {
		my $db_name = $self->{db_name};
		unless ($db_name) { die "\n\t Error: no DB name defined \n\n\n"; }
		my $db_obj = ScreeningDB->new($self);
		$db_obj->load_screening_db($db_name);	
		$self->{db} = $db_obj; # Store the database object reference 
	}

	# Hand off to functions 
	if ($option eq 1)    { # Create a screening DB 
		my $target_db_obj = TargetDB->new($self);
		$target_db_obj->summarise_genomes_short();
	}
	elsif ($option eq 2)    { # Create a screening DB 
		my $target_db_obj = TargetDB->new($self);
		$target_db_obj->summarise_genomes_long();
	}
	else {
		print "\n\t  Unrecognized option '-m=$option'\n";

	}
}

############################################################################
# MAIN FUNCTIONS
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
	my $db_ref = $self->{db};
	unless ($db_ref) { die; }  # Sanity checking

	# Set up the screening queries
	my %queries;
	my $loader_obj = $self->{loader_obj};
	unless ($loader_obj) { die; }  # Sanity checking
	my $total_queries = $loader_obj->set_up_screen($self, \%queries);
	unless ($total_queries)  { 
		print "\n\t  Exiting without screening.\n\n";	
		exit;
	}
	$self->{total_queries} = $total_queries;

	# Iterate through and excute screens
	print "\n\t  Starting database-integrated genome screening\n";
	my @probes = keys %queries;

	my $num_queries = 0;
	my %crossmatching;
	$self->{crossmatching} = \%crossmatching;
	foreach my $probe_name (@probes) {
		
		# Get the array of queries for this target file
		my $probe_queries = $queries{$probe_name};
		foreach my $query_ref (@$probe_queries) {  
	
			# Increment query count
			$num_queries++;		
			$self->{num_queries} = $num_queries;

			# Do the 1st BLAST (probe vs target)
			$self->search($query_ref);	
			# Do the 2nd BLAST (hits from 1st BLAST vs reference library)
			$self->assign($query_ref);
		}	
	}
	
	# Cleanup
	my $output_dir = $loader_obj->{report_dir};
	my $command1 = "rm -rf $output_dir";
	system $command1;

	# Show cross matches at end
	$self->show_cross_matching();

	# Print finished message
	print "\n\n\t ### SCREEN COMPLETE ~ + ~ + ~";
}

#***************************************************************************
# Subroutine:  search
# Description: execute a search (i.e. run a BLAST query)
# Arguments:   $query_ref - data structure with the query details
#***************************************************************************
sub search {
	
	my ($self, $query_ref) = @_;

	# Get relevant member variables and objects
	my $db_ref          = $self->{db};
	my $blast_obj       = $self->{blast_obj};
	my $tmp_path        = $self->{tmp_path};
	my $min_length      = $self->{seq_length_minimum};
	my $redundancy_mode = $self->{redundancy_mode};
	my $total_queries   = $self->{total_queries};
	my $num_queries     = $self->{num_queries};	
	unless ($num_queries and $total_queries) { die; }
	my $percent_prog    = ($num_queries / $total_queries) * 100;
	my $f_percent_prog  = sprintf("%.2f", $percent_prog);

	# Sanity checking
	unless ($min_length)      { die "\n\t Default_min_seqlen undefined\n\n\n"; }
	unless ($blast_obj)       { die; } 
	unless ($tmp_path)        { die; } 
	unless ($redundancy_mode) { die; } 
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

	# Do the BLAST search
	print "\n\t  Screen: $num_queries (%$f_percent_prog done): '$organism' file '$target_name' with probe $probe_id";   
	$blast_obj->blast($blast_alg, $target_path, $probe_path, $result_file);
	
	# Parse out the alignment
	my @hits;
	$blast_obj->parse_tab_format_results($result_file, \@hits, $cutoff);

	# TODO: catch error from BLAST and don't update Status table	

	# Clean up - remove the result file
	my $rm_command = "rm $result_file";
	system $rm_command;

	# Rename coordinate fields to match DB
	foreach my $hit_ref (@hits) {
		$hit_ref->{subject_start} = $hit_ref->{aln_start};
		$hit_ref->{subject_end}   = $hit_ref->{aln_stop};
		$hit_ref->{query_end}     = $hit_ref->{query_stop};
	}

	# Store new new BLAST hits that meet conditions
	my $blast_results_table = $db_ref->{blast_results_table};
	my $i = 0;
	my $num_hits = scalar @hits;
	if ($num_hits > 0) {
		print "\n\t\t # $num_hits matches to probe: $probe_name, $probe_gene";
	}
	foreach my $hit_ref (@hits) {
	
		$i++;
		
		# Skip sequences that are too short
		my $start  = $hit_ref->{subject_start};
		my $end    = $hit_ref->{subject_end};
		if ($end - $start < $min_length) {  next; }
		
		# Record hit in BLAST results table
		$hit_ref->{organism}     = $organism;
		$hit_ref->{version}      = $version;
		$hit_ref->{data_type}    = $data_type;
		$hit_ref->{target_name}  = $target_name;
		$hit_ref->{probe_id}     = $probe_id;
		$hit_ref->{probe_name}   = $probe_name;
		$hit_ref->{probe_gene}   = $probe_gene;
		$hit_ref->{probe_type}   = $query_ref->{probe_type};
		$hit_ref->{hit_length}   = $hit_ref->{align_len};
		$blast_results_table->insert_row($hit_ref);
	} 

	# Consolidate the BLAST table based on results
	if ($redundancy_mode > 1) {
		$self->consolidate_hits($query_ref);
	}
	
	# Update the status table
	my $status_table = $db_ref->{status_table};
	$status_table->insert_row($query_ref);
}

#***************************************************************************
# Subroutine:  assign
# Description: assign sequences that matched probes in a BLAST search 
#***************************************************************************
sub assign {
	
	my ($self, $query_ref) = @_;
	
	# Get parameters from self
	my $db_ref = $self->{db};
	my $table  = $db_ref->{extracted_table}; 

	# Extract hits 
	my @extracted;
	my $new_sequences = $self->extract_unassigned_hits($query_ref, \@extracted);	
	
	# Iterate through the matches
	my $assigned_count = 0;
	my $crossmatch_count = 0;
	foreach my $hit_ref (@extracted) {

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
		if ($probe_key ne $assigned) {
			#print "\n\t\t #   cross match to '$assigned'";
			$crossmatch_count++;
			$self->update_cross_matching($probe_key, $assigned);
		}

		# Insert the data to the Extracted table
		my $extract_id = $table->insert_row(\%data);
	}

	print "\n\t\t # $new_sequences newly identified hits";
	if ($new_sequences > 0) {
		print "\n\t\t # $assigned_count extracted sequences matched to reference library";
		print "\n\t\t # $crossmatch_count cross-matched to something other than the probe";
	}
}

#***************************************************************************
# Subroutine:  reassign
# Description: reassign sequences in the extracted_table (for use after
#              the reference library has been updated)
#***************************************************************************
sub reassign {
	
	my ($self) = @_;

	# Set up to perform the reassign process
	my @assigned_seqs;
	$self->initialise_reassign(\@assigned_seqs);

	# Get data structures and variables from self
	my $blast_obj       = $self->{blast_obj};
	my $result_path     = $self->{report_dir};
	my $db              = $self->{db};
	my $extracted_table = $db->{extracted_table};
	unless ($extracted_table) { die; }
	
	# Iterate through the matches
	print "\n\n\t  Reassigning Extracted table\n";
	my $count = 0;
	foreach my $hit_ref (@assigned_seqs) {

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
			# Insert the data
			print "\n\t ##### Reassigned $blast_id from $previous_assign ($previous_gene)";
			print " to $assigned_name ($assigned_gene)";
			my $where = " WHERE Record_id = $blast_id ";
			$extracted_table->update($hit_ref, $where);
		}
	}
	# Cleanup
	my $output_dir = $self->{report_dir};
	my $command1 = "rm -rf $output_dir";
	system $command1;
}

############################################################################
# SECTION: Extend the screening database by adding ancillary tables
############################################################################

#***************************************************************************
# Subroutine:  extend_screening_db
# Description: console managemant of ancillary tables in the screening database
#***************************************************************************
sub extend_screening_db {

	my ($self) = @_;

	# Try to read the tab-delimited infile
	print "\n\n\t #### WARNING: This function expects a tab-delimited data table with column headers!";
	my $question1 = "\n\n\t Please enter the path to the file with the table data and column headings\n\n\t";
	my $infile = $console->ask_question($question1);
	unless ($infile) { die; }
	my @infile;
	$fileio->read_file($infile, \@infile);
	my @data;
	my $line_number = 0;
	foreach my $line (@infile) {
		$line_number++;
		if     ($line =~ /^\s*$/)  { next; } # discard blank line
		elsif  ($line =~ /^\s*#/)  { next; } # discard comment line 
		unless ($line =~ /\t/)     { die; }
		push (@data, $line);
	}
	my $header_row = shift @data;
	my @header_row = split ("\t", $header_row);
	print "\n\n\t The following column headers (i.e. table fields) were obtained\n";
	my $i;
	my %fields;
	my @fields;
	foreach my $element (@header_row) {
		chomp $element;
		$i++;
		$element =~ s/\s+/_/g;
		if ($element eq '') { $element = 'EMPTY_COLUMN_' . $i; } 
		print "\n\t\t Column $i: '$element'";
		push (@fields, $element);
		$fields{$element} = "varchar";
	}
	my $question3 = "\n\n\t Is this correct?";
	my $answer3 = $console->ask_yes_no_question($question3);
	if ($answer3 eq 'n') { # Exit if theres a problem with the infile
		print "\n\t\t Aborted!\n\n\n"; exit;
	}

	# Create new table or add to existing one
	my $db = $self->{db};
	unless ($db) { die; }
	my %extra_tables;
	my @extra_tables;
	my $table_to_use;

	my @choices = qw [ 1 2 3 4 ];
	print "\n\t\t 1. Create new ancillary table";
	print "\n\t\t 2. Append data to existing ancillary table";
	print "\n\t\t 3. Flush existing ancillary table and upload fresh data";
	print "\n\t\t 4. Drop an ancillary table\n";
	my $question4 = "\n\t Choose an option";
	my $answer4   = $console->ask_simple_choice_question($question4, \@choices);
	if ($answer4 == '1') {	# Create new table
		my $table_name_question = "\n\t What is the name of the new table?";
		my $table_name = $console->ask_question($table_name_question);
		$table_to_use = $db->create_ancillary_table($table_name, \@fields, \%fields);	
	}
	else {
		my @choices;
		$db->get_ancillary_table_names(\@extra_tables);
		my $table_num = 0;
		foreach my $table_name (@extra_tables) {
			$table_num++;
			$extra_tables{$table_num} = $table_name;
			print "\n\t\t Table $table_num: '$table_name'";
			push (@choices, $table_num);
		}
		my $question5 = "\n\n\t Apply to which of the above tables?";
		my $answer5   = $console->ask_simple_choice_question($question5, \@choices);
		$table_to_use = $extra_tables{$answer5};
	}
	unless ($table_to_use) { die; }
	my $dbh = $db->{dbh};
	unless ($dbh) { die "\n\t Couldn't retrieve database handle \n\n"; }
	my $anc_table = MySQLtable->new($table_to_use, $dbh, \%fields);
    $db->{$table_to_use} = $anc_table;

	if ($answer4 == '4') {	# Drop an ancillary table
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
			#if ($verbose) {
				print "\n\t Row count $row_count: uploading value '$value' to field '$field'";
			#}
			unless ($value) { 
				$value = 'NULL';
			}
			$insert{$field} = $value;
		}
		$anc_table->insert_row(\%insert);
	}
}

############################################################################
# ANCILLARY FUNCTIONS 
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

	# Store the ScreenBuilder object (used later)
	$self->{loader_obj} = $loader_obj; 
}

#***************************************************************************
# Subroutine:  initialise_reassign 
# Description: set up for reassigning the sequences in the Extracted table
#***************************************************************************
sub initialise_reassign {

	my ($self, $assigned_seqs_ref) = @_;

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
	$extracted_table->select_rows(\@fields, $assigned_seqs_ref);

	# Set up the reference library
	$loader_obj->setup_reference_library($self);

	# Transfer parameters from loader to this obj
	$self->{seq_length_minimum}    = $loader_obj->{seq_length_minimum};
	$self->{bit_score_min_tblastn} = $loader_obj->{bit_score_min_tblastn};
	$self->{bit_score_min_blastn}  = $loader_obj->{bit_score_min_blastn};
	$self->{blast_orf_lib_path}    = $loader_obj->{blast_orf_lib_path};
	$self->{blast_utr_lib_path}    = $loader_obj->{blast_utr_lib_path};
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
	my $assigned      = 1;

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

	return $assigned;
}

#***************************************************************************
# Subroutine:  extract_unassigned_hits
# Description: extract BLAST hit sequences that are not in Extracted table
#***************************************************************************
sub extract_unassigned_hits {
	
	my ($self, $query_ref, $extracted_ref) = @_;

	# Get paths, objects, data structures and variables from self
	my $blast_obj   = $self->{blast_obj};
	my $target_path = $query_ref->{target_path};
	my $target_name = $query_ref->{target_name};
	my $db_ref      = $self->{db};
	my $blast_results_table = $db_ref->{blast_results_table};

	# Index extracted BLAST results 
	my %extracted;
	my $where = " WHERE target_name = '$target_name' ";
	$db_ref->index_extracted_loci_by_blast_id(\%extracted, $where);

	# Get hits we are going to extract
	my @hits;
	$db_ref->get_blast_hits_to_extract($query_ref, \@hits);
	my $num_hits = scalar @hits;
	#print "\n\t ### There are $num_hits hits to extract";

	# Store all outstanding matches as sequential, target-ordered sets 
	my $new_hits = scalar @hits;
	foreach my $hit_ref (@hits) {
		
		# Skip previously extracted hits
		my $record_id   = $hit_ref->{record_id};
		if ($verbose) {
			print "\n\t Checking $record_id in extracted table";
		}
		if ($extracted{$record_id}) { 
			if ($verbose) {
				print "\t ALREADY EXTRACTED";
			}
			next;
		 }
		# Extract the sequence
		my $sequence = $blast_obj->extract_sequence($target_path, $hit_ref);
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
# Subroutine:  consolidate_hits
# Description: Consolidate the BLAST table based on results
#***************************************************************************
sub consolidate_hits {
	
	my ($self, $query_ref) = @_;

	# Get relevant member variables and objects
	my $db_ref  = $self->{db};
	my $blast_results_table = $db_ref->{blast_results_table};
	my $blast_chains_table  = $db_ref->{blast_chains_table};

	# Set the fields to get values for
	my @fields = qw [ record_id scaffold orientation
	                  subject_start subject_end
                      query_start query_end ];

	# Get the information for this query
	my $target_name = $query_ref->{target_name};
	my $probe_name  = $query_ref->{probe_name};
	my $probe_gene  = $query_ref->{probe_gene};

	my $where  = " WHERE target_name = '$target_name'";

	# Filter by gene name 
	if ($self->{redundancy_mode} eq 2) {
		$where .= " AND probe_gene = '$probe_gene' ";
	}
	elsif ($self->{redundancy_mode} eq 3) {
		$where .= " AND probe_name = '$probe_name'
                    AND probe_gene = '$probe_gene' ";
	}

	# Order by ascending start coordinates within each scaffold in the target file
	$where .= "ORDER BY scaffold, subject_start ";
	
	# Get the relevant loci
	my @hits;
	$blast_results_table->select_rows(\@fields, \@hits, $where);

	# Iterate through consolidating as we go
	my $i;
	my %last_hit;
	my %consolidated;
	my %retained;
	my $consolidated_count = 0;
	foreach my $hit_ref (@hits)  {

		# Get hit values
		$i++;
		my $record_id     = $hit_ref->{record_id};
		my $scaffold      = $hit_ref->{scaffold};
		my $orientation   = $hit_ref->{orientation};
		my $subject_start = $hit_ref->{subject_start};
		
		# Get last hit values
		my $last_record_id     = $last_hit{record_id};
		my $last_scaffold      = $last_hit{scaffold};
		my $last_orientation   = $last_hit{orientation};
		my $consolidated;

		if ($verbose) {
			print "\n\t $subject_start: RECORD ID $record_id";
		}

		# Keep if first in this loop process
		if ($i eq 1) {
			$retained{$record_id} = 1;
		}
		# ...or if first hit on scaffold
		elsif ($scaffold ne $last_scaffold) {
			$retained{$record_id} = 1;
		}
		# ...or if in opposite orientation to last hit
		elsif ($orientation ne $last_orientation) {
			$retained{$record_id} = 1;
		}
		else { # If we get this far we have two hits on the same scaffold
			
			# Check whether to consolidate hit on the same target scaffold
			$consolidated = $self->inspect_adjacent_hits(\%last_hit, $hit_ref);

			# Keep track of the outcome
			if ($consolidated) {
				$consolidated_count++;
				my %hit = %last_hit; # Make a copy
				$hit{hit_length} = ($hit{subject_end} - $hit{subject_start}) + 1;
				$consolidated{$record_id} = \%hit;
				$retained{$last_record_id} = 1;
			}
			else { $retained{$record_id} = 1; }
		}
		
		# Update the 'last hit' to the current one before exiting this iteration
		unless ($consolidated) {
			$last_hit{record_id}     = $record_id;
			$last_hit{scaffold}      = $scaffold;
			$last_hit{orientation}   = $orientation;
			$last_hit{subject_start} = $hit_ref->{subject_start};
			$last_hit{subject_end}   = $hit_ref->{subject_end};
			$last_hit{query_start}   = $hit_ref->{query_start};
			$last_hit{query_end}     = $hit_ref->{query_end};
		}
	}
	# Update BLAST table
	$self->update_db_loci(\@hits, \%retained, \%consolidated);
}

#***************************************************************************
# Subroutine:  inspect_adjacent_hits
# Description: determine whether two adjacent hits should be joined
#***************************************************************************
sub inspect_adjacent_hits {
	
	my ($self, $last_hit_ref, $hit_ref) = @_;

	# Get parameters for consolidating hits from self
	my $probe_buffer = $self->{threadhit_probe_buffer};
	my $gap_buffer   = $self->{threadhit_gap_buffer};
	my $max_gap      = $self->{threadhit_max_gap};
	unless ($probe_buffer and $gap_buffer and $max_gap) { die; } 

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

	#print "\n\n\t #### Checking whether to consolidate $last_record_id and $record_id on $scaffold";
	#print "\n\t #### Q Last:  $last_query_start\t $last_query_end";
	#print "\n\t #### S Last:  $last_subject_start\t $last_subject_end";
	#print "\n\t #### Q This:  $query_start\t $query_end";
	#print "\n\t #### S This:  $subject_start\t $subject_end";
	#print "\t #### Gap:   $subject_gap\n";

	# Calculate the gap between the query coordinates of the two matches
	# Note this may be a negative number if the queries overlap 
	my $query_gap;
	
	# TO REMOVE
	if ($orientation eq '+ve') { $orientation = '+'; }
	if ($orientation eq '-ve') { $orientation = '-'; }
	if ($orientation eq '+') {
		$query_gap = $query_start - $last_query_end;
	}
	elsif ($orientation eq '-') {
		$query_gap = $last_query_start - $query_end;
	}
	
	else { print "\n\t orientation = $orientation\n\n";  die; }
	
	### Deal with contingencies that mean hits definitely should or should not be consolidated

	#   1. Hit is entirely within a previous hit it is redundant
	if ($last_subject_start <= $subject_start and $last_subject_end >= $subject_end) {
		return 1; # Effectively discarding current hit
	}
	#   2. Hits are too far apart
	elsif ($max_gap) {
		if ($subject_gap > $max_gap) {
			return 0;  
		}
	}

	### Deal with situations where hits are partially (but not completely) overlapping
	#   or where they are close enough to each other consider consolidating them
	
	# Set the position based on threadhit_probe_biffer
	my $buffer_start = $query_start + $probe_buffer;

	if ($subject_gap < 1) {
		# Describe the consolidation:
		if ($verbose) {
			print "\t    SUBJECT GAP: Merged hit $record_id into hit $last_record_id (subject_gap: $subject_gap)";
		}
		$last_hit_ref->{subject_end} = $hit_ref->{subject_end};
		$last_hit_ref->{query_end}   = $hit_ref->{query_end};
		return 1;  # Consolidated (joined) hits
	}
	

	# For positive orientation hits
	# If the query_start of this hit is before the query_end of the last, then its a distinct hit
	elsif ($orientation eq '+' and $buffer_start < $last_query_end) {
		if ($verbose) {
			print "\n\t    DISTINCT +: start $query_start (buffer $buffer_start) is too much less than end of last query ($last_query_end) \tsubject gap: $subject_gap";
		}
		return 0;
	}
	# For negative orientation hits
	# If the query_start of this hit < query_end of the last, then its a distinct hit
	elsif ($orientation eq '-' and ($last_query_start - $query_end) > $probe_buffer ) {
		if ($verbose) {
			print "\n\t    DISTINCT -: ($last_query_start - $query_end) > $probe_buffer \tsubject gap: $subject_gap";
		}
		return 0;
	}

	# If not, then check the intervening distance between the hits
	else {

		my $nt_query_gap = $query_gap * 3;
		if (($subject_gap - $nt_query_gap) < $gap_buffer) {
			# Describe the consolidation:
			if ($verbose) {
				print "\n\t    Merged hit $record_id into hit $last_record_id: subject gap ($subject_gap) nt query gap ($nt_query_gap) ";
			}
			$last_hit_ref->{subject_end} = $hit_ref->{subject_end};
			$last_hit_ref->{query_end}   = $hit_ref->{query_end};
			return 1;  # Consolidated (joined) hits
		}
		else {
			if ($verbose) {
				print "\n\t    DISTINCT LAST: start $query_start (buffer $buffer_start) is too much less than end of last query ($last_query_end) \tsubject gap: $subject_gap";
			}
			return 0;  # Distinct hit
		}
	}
}

#***************************************************************************
# Subroutine:  update_db_loci  
# Description: update the BLAST results table with overlapping hit info
#***************************************************************************
sub update_db_loci {
	
	my ($self, $hits_ref, $retained_ref, $consolidated_ref) = @_;

	# Get relevant member variabless
	my $db_ref  = $self->{db};
	my $blast_results_table = $db_ref->{blast_results_table};
	my $extracted_table     = $db_ref->{extracted_table};

	# Update BLAST table with consolidated hits
	my $blast_updated = 0;
	my $blast_deleted = 0;
	my $extract_deleted = 0;
	my @ids = keys %$consolidated_ref;
	foreach my $record_id (@ids)  {

		my $where   = " WHERE record_id = $record_id ";
		my $hit_ref = $consolidated_ref->{$record_id};		
		delete $hit_ref->{record_id};
		#print "\n\t ## updating hit record ID $record_id in BLAST results";
		my $this_id = $blast_results_table->update($hit_ref, $where);
		$blast_updated++;
	
		# Update the 'BLAST_chains' table accordingly
		#$hit{this_query_id} = $this_id; 
	
		# Delete from Extracted (because we need to re-extract and assign)
		#print "\n\t ## deleting hit BLAST ID $record_id in Extracted";
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

			#print "\n\t ## deleting hit BLAST ID $record_id in Extracted";
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
		$HELP  .= "\n\t -m=4  Flush screening DB"; 
		$HELP  .= "\n\t -m=5  Drop screening DB"; 
		$HELP  .= "\n\t -m=6  Manage ancillary tables"; 
		$HELP  .= "\n\t -m=7  Format genome directory ($ENV{DIGS_GENOMES})\n"; 
        $HELP  .= "\n\t ### Utility functions"; 
		$HELP  .= "\n\t -u=1  Summarise the genome directory ($ENV{DIGS_GENOMES})\n\n"; 
	print $HELP;
}

############################################################################
# EOF
############################################################################
