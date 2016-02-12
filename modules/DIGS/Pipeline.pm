#!/usr/bin/perl -w
############################################################################
# Module:      Pipeline.pm
# Description: Genome screening pipeline using reciprocal BLAST
# History:     December 2009: Created by Robert Gifford 
############################################################################
package Pipeline;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;
############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::FileIO;
use Base::DevTools;
use Base::Console;
use Base::Sequence;    # For performing basic sequence manipulations

# Program components
use DIGS::ScreenBuild; # Functions to set up screen
use DIGS::Retrieve; 

############################################################################
# Globals
############################################################################

# Default minimum length of BLAST hit to extract 
my $default_min_seqlen = 100; # minimum sequence length of BLAST hit to extract
my $process_dir_warn   = 100; # folders in process directory to trigger warning 

# Base objects
my $fileio    = FileIO->new();
my $devtools  = DevTools->new();
my $console   = Console->new();
my $verbose   = undef;
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create new Pipeline 'object'
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {
		
		# Flags
		process_id             => $parameter_ref->{process_id},
		program_version        => $parameter_ref->{program_version},
		refresh_genomes        => $parameter_ref->{refresh_genomes},
		verbose                => $parameter_ref->{verbose},
		
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
# SECTION: DIGS top level handler subroutines
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
	unless ($option > 10) { # Validate configurations that require a control file
		# An infile must be defined
		unless ($ctl_file) { 
			die "\n\t Option '$option' requires an infile\n\n";
		}
		$self->initialise($ctl_file);
	}

	# If it exists, load the screening database specified in the control file
	if ( $option > 1) {
		my $db_name = $self->{db_name};
		unless ($db_name) { die "\n\t Error: no DB name defined \n\n\n"; }
		my $db_obj = ScreeningDB->new($self);
		$db_obj->load_screening_db($db_name);	
		$self->{db} = $db_obj; # Store the database object reference 
	}

	# Hand off to subroutines 
	if ($option eq 1)    { # Create a screening DB 
		$self->create_screening_db($ctl_file);
	}
	elsif ($option eq 2) { # Screen
		$self->screen();	
	}
	elsif ($option eq 3) { # Defragment
		$self->defragment();	
	}
	elsif ($option eq 4) { # Reassign data in Exracted table
		$self->reassign();	
	}
	elsif ($option eq 5) { # Flush screening DB
		my $db = $self->{db};
		$db->flush_screening_db();
	}
	elsif ($option eq 6) { # Drop screening DB 
		my $db = $self->{db};
		$db->drop_screening_db();    
	}
	elsif ($option eq 7) { # Summarise genomes 
		my $retrieve_obj = Retrieve->new($self);
		$retrieve_obj->run_data_retrieval_functions();
	}
	elsif ($option eq 8) { # Add a table of data to the screening database
		$self->extend_screening_db();
	}
}

############################################################################
# SECTION: Initialise
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
	print "\n\t ### Reading control file\n";
	my $loader_obj = ScreenBuild->new($self);
	$loader_obj->parse_control_file($ctl_file, $self);

	# Store the ScreenBuild.pm object (used later for some configurations)
	$self->{loader_obj} = $loader_obj; 
}

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

############################################################################
# SECTION: Running a round of paired BLAST screening
############################################################################

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
	print "\n\n\t ### Setting up a DIGS screen";
	my %queries;
	my $loader_obj = $self->{loader_obj};
	unless ($loader_obj) { die; }  # Sanity checking
	my $valid = $loader_obj->set_up_screen($self, \%queries);
	unless ($valid) {
		print "\n\n\t ### Exiting without screening\n\n\n";
		exit;
	}
	#$devtools->print_hash(\%queries); die; # DEBUG

	# Iterate through and excute screens
	print "\n\n\t ### Starting database-integrated genome screening\n";
	my @probes = keys %queries;
	foreach my $probe_name (@probes) {
		
		# Get the array of queries for this target file
		my $probe_queries = $queries{$probe_name};
		foreach my $query_ref (@$probe_queries) {  
			
			# Show status 
			my $probe_id    = $query_ref->{probe_id};
			my $organism = $query_ref->{organism}; 
			my $target_name = $query_ref->{target_name}; # name of target sequence file
			print "\n\n\t # Screening '$organism' file '$target_name' with probe $probe_id ";   
			
			# Do the 1st BLAST (probe vs target)
			$self->search($query_ref);	
			print "\n\t # 1st BLAST step completed, BLAST results table up-to-date...";

			# Do the 2nd BLAST (hits from 1st BLAST vs reference library)
			$self->assign($query_ref);
			print "\n\t # 2nd BLAST step completed, Extracted table up-to-date...";

			# Summarise extracted table
			#$db_ref->summarise_extracted_table();	
		}	
	}
	
	# Cleanup
	my $output_dir = $loader_obj->{report_dir};
	my $command1 = "rm -rf $output_dir";
	system $command1;

	# Print finished message
	print "\n\n\n\t ### SCREEN COMPLETE ~ + ~ + ~\n\n\n";
}

#***************************************************************************
# Subroutine:  search
# Description: execute a search (i.e. run a BLAST query)
# Arguments:   $query_ref - data structure with the query details
#***************************************************************************
sub search {
	
	my ($self, $query_ref) = @_;

	# Get relevant member variables and objects
	my $db_ref       = $self->{db};
	my $blast_obj    = $self->{blast_obj};
	my $tmp_path     = $self->{tmp_path};
	my $min_length   = $self->{seq_length_minimum};
	my $redundancy_mode = $self->{redundancy_mode};
	unless ($min_length) { $min_length = $default_min_seqlen; }
	
	# Sanity checking
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
	#$devtools->print_hash($query_ref); die; # DEBUG	

	# Do the BLAST search
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
	print "\n\t\t # $num_hits matches to probe: $probe_name, $probe_gene";
	print "\n\t\t # in  $organism, $data_type, $version, '$target_name'";
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
		if ($verbose) {
			print "\n\t # Match $i to probe: $probe_name, $probe_gene";
			print "\n\t # - in genome: $organism, $data_type, $version";
			print "\n\t # - target file: '$target_name'";
		}
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

############################################################################
# SECTION: Assign 
############################################################################

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
	$self->extract_unassigned_hits($query_ref, \@extracted);	
	
	# Iterate through the matches
	my $assigned_count = 0;
	foreach my $hit_ref (@extracted) {

		# Set the linking to the BLAST result table
		my $blast_id  = $hit_ref->{record_id};
		$hit_ref->{blast_id} = $blast_id;
		
		# Execute the 'reverse' BLAST (2nd BLAST in a round of paired BLAST)	
		my %data = %$hit_ref; # Make a copy
		$self->do_reverse_blast(\%data);
		$assigned_count++;

		# Insert the data
		my $extract_id = $table->insert_row(\%data);
	}

	# TODO: validate?

}

#***************************************************************************
# Subroutine:  reverse BLAST
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
		unless ($lib_path) { 
			print "\n\t NO UTR LIBRARY defined"; 
		}
	}
	elsif ($probe_type eq 'ORF') {
		$lib_path  = $self->{blast_orf_lib_path};
		unless ($lib_path) { 
			print "\n\t NO ORF LIBRARY defined"; 
		}
		$blast_alg = 'blastx';
	}
	else { die; }
	unless ($lib_path) { return; }

	# Execute the 'reverse' BLAST (2nd BLAST in a round of paired BLAST)	
	$blast_obj->blast($blast_alg, $lib_path, $query_file, $result_file);
	my @results;
	$blast_obj->parse_tab_format_results($result_file, \@results);
	#$devtools->print_array(\@results); die; # DEBUG 	

	# Get the best match from this file
	my $top_match = shift @results;
	my $query_start   = $top_match->{query_start};
	my $query_end     = $top_match->{query_stop};
	my $subject_start = $top_match->{aln_start};
	my $subject_end   = $top_match->{aln_stop};
	my $assigned_name = $top_match->{scaffold};	
	#$devtools->print_hash($top_match); # DEBUG
	
	unless ($assigned_name) {	
		
		print "\n\t ### No match found in reference library\n";
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
	}
	else {	

		# Split assigned to into (i) refseq match (ii) refseq description (e.g. gene)	
		my @assigned_name = split('_', $assigned_name);
		my $assigned_gene = pop @assigned_name;
		$assigned_name = join ('_', @assigned_name);
		#print " assigned to: $assigned_name: $assigned_gene!";
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
	}

	# Clean up
	my $command1 = "rm $query_file";
	my $command2 = "rm $result_file";
	system $command1;
	system $command2;
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
}

#***************************************************************************
# Subroutine:  reassign
# Description: reassign sequences in the extracted_table (for use after
#              the reference library has been updated)
#***************************************************************************
sub reassign {
	
	my ($self) = @_;

	# Set up to perform the reassign process
	print "\n\t # Reassigning Extracted table";
	my @assigned_seqs;
	$self->initialise_reassign(\@assigned_seqs);

	# Get data structures and variables from self
	my $blast_obj       = $self->{blast_obj};
	my $result_path     = $self->{report_dir};
	my $db              = $self->{db};
	my $extracted_table = $db->{extracted_table};
	unless ($extracted_table) { die; }
	
	# Iterate through the matches
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
		print "\n\t Redoing assign for record ID $blast_id assigned to $previous_assign";
		print "\n\t coordinates: $extract_start-$extract_end";
		$self->do_reverse_blast($hit_ref);
		
		my $assigned_name = $hit_ref->{assigned_name};
		my $assigned_gene = $hit_ref->{assigned_gene};
		if ($assigned_name ne $previous_assign 
		or  $assigned_gene ne $previous_gene) {
			
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
	if ($loader_obj->{reference_aa_fasta}) {
		$loader_obj->load_aa_fasta_reference_library();
	}
	if ($loader_obj->{reference_nt_fasta}) {
		$loader_obj->load_nt_fasta_reference_library();
	}

	# Transfer parameters from loader to this obj
	$self->{seq_length_minimum}    = $loader_obj->{seq_length_minimum};
	$self->{bit_score_min_tblastn} = $loader_obj->{bit_score_min_tblastn};
	$self->{bit_score_min_blastn}  = $loader_obj->{bit_score_min_blastn};
	$self->{blast_orf_lib_path}    = $loader_obj->{blast_orf_lib_path};
	$self->{blast_utr_lib_path}    = $loader_obj->{blast_utr_lib_path};
}

#***************************************************************************
# Subroutine:  reassign
# Description: reassign sequences in the extracted_table (for use after
#              the reference library has been updated)
#***************************************************************************
sub reassign2 {
	
	my ($self, $reextract, $reassign_set_ids) = @_;

	# Set up to perform the reassign process
	print "\n\t # Regenerating Extracted table";
	my @reassign_set;
	$self->initialise_reassign2($reextract, \@reassign_set, $reassign_set_ids);
	#$devtools->print_array(\@reassign_set); die;

	# Get data structures and variables from self
	my $blast_obj       = $self->{blast_obj};
	my $result_path     = $self->{report_dir};
	my $db              = $self->{db};
	my $extracted_table = $db->{extracted_table};
	unless ($extracted_table) { die; }
	
	# Iterate through the matches
	foreach my $hit_ref (@reassign_set) {

		# Set the linking to the BLAST result table
		my $blast_id        = $hit_ref->{blast_id};
		my $record_id       = $hit_ref->{record_id};

		# Execute the 'reverse' BLAST (2nd BLAST in a round of paired BLAST)	
		my $previous_assign = $hit_ref->{assigned_name};
		my $previous_gene   = $hit_ref->{assigned_gene};
		print "\n\t Redoing assign for BLAST record '$blast_id'";
		print ", previously assigned to $previous_assign ($previous_gene)";
		$self->do_reverse_blast($hit_ref);
		
		my $assigned_name = $hit_ref->{assigned_name};
		my $assigned_gene = $hit_ref->{assigned_gene};
		if ($assigned_name ne $previous_assign 
		or  $assigned_gene ne $previous_gene) {
			
			# Remove the fields we don't need to update 
			delete $hit_ref->{record_id};
			delete $hit_ref->{sequence};
		
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


#***************************************************************************
# Subroutine:  initialise_reassign2 
# Description: set up for reassigning the sequences in the Extracted table
#***************************************************************************
sub initialise_reassign2 {

	my ($self, $reextract, $assign_set, $assign_ids) = @_;

	# Create a unique ID and report directory for this run
	my $output_path = $self->{output_path};
	my $process_id  = $self->{process_id};
	my $db          = $self->{db};
	my $blast_obj   = $self->{blast_obj};
	my $db_name     = $db->{db_name};
	unless ($db and $db_name and $process_id and $output_path and $blast_obj) { die; }
	
	# Create report directory
	my $loader_obj = $self->{loader_obj};
	$loader_obj->create_output_directories($self);


	# Get the set of rows to reassign
	my @reassign_rows;
	if ($assign_ids) { 
		foreach my $data_ref (@$assign_ids) {
			my $record_id = $data_ref->{record_id};
			unless ($record_id) { $devtools->print_hash($data_ref); die; }
			my @rows;
			my $blast_table = $db->{blast_results_table};
			my @fields  = qw [ record_id probe_name probe_gene probe_type organism 
                               data_type version target_name scaffold 
                               subject_start subject_end orientation
                               query_start query_end hit_length  ];
			my $where = "WHERE Record_ID = $record_id";
			$blast_table->select_rows(\@fields, \@rows, $where);
			my $data_ref2 = shift @rows;
			my $genome_path = $self->{genome_use_path};
			my $organism    = $data_ref2->{organism};
			my $data_type   = $data_ref2->{data_type};
			my $version     = $data_ref2->{version};
			my $target_name = $data_ref2->{target_name};
			my $target_path = $genome_path . "Mammalia/$organism" . "/$data_type" . "/$version" . "/$target_name";

			# Extract the sequence
			my $sequence = $blast_obj->extract_sequence($target_path, $data_ref2);
			if ($sequence) {	
				my $seq_length = length $sequence; # Set sequence length
				$data_ref2->{sequence_length} = $seq_length;
				$data_ref2->{sequence} = $sequence;
				push (@$assign_set, $data_ref2);
			}
			else { die; }

			#$devtools->print_hash($data_ref2); die; # DEBUG
			push (@$assign_set, $data_ref2);
		}
	}
	else {   # If no set of rows to reassign are in the array, reassign the entire Extracted table
		# Get the minumum set of fields required for a reverse blast update of an Extracted table row
		my $extracted_table = $db->{extracted_table};
		my @fields  = qw [ record_id  blast_id probe_type
                           assigned_name assigned_gene sequence ];
		$extracted_table->select_rows(\@fields, $assign_set);
	}

	# Set up the reference library
	if ($loader_obj->{reference_aa_fasta}) {
		$loader_obj->load_aa_fasta_reference_library();
	}
	if ($loader_obj->{reference_nt_fasta}) {
		$loader_obj->load_nt_fasta_reference_library();
	}
	if ($loader_obj->{reference_glue}) {
		$loader_obj->load_glue_reference_library();
	}

	# Transfer parameters from loader to this obj
	$self->{seq_length_minimum}    = $loader_obj->{seq_length_minimum};
	$self->{bit_score_min_tblastn} = $loader_obj->{bit_score_min_tblastn};
	$self->{bit_score_min_blastn}  = $loader_obj->{bit_score_min_blastn};
	$self->{blast_orf_lib_path}    = $loader_obj->{blast_orf_lib_path};
	$self->{blast_utr_lib_path}    = $loader_obj->{blast_utr_lib_path};
}

############################################################################
# SECTION: Consolidate Fxns
############################################################################

#***************************************************************************
# Subroutine:  consolidate_hits
# Description: Consolidate the BLAST table based on results
#***************************************************************************
sub consolidate_hits {
	
	my ($self, $query_ref) = @_;

	# Get relevant member variables and objects
	my $db_ref  = $self->{db};
	my $blast_results_table = $db_ref->{blast_results_table};

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
			else {
				$retained{$record_id} = 1;
			}
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

	# Get relevant member variables and objects
	my $db_ref  = $self->{db};
	my $blast_results_table = $db_ref->{blast_results_table};
	my $extracted_table     = $db_ref->{extracted_table};

	# Delete redundant rows from BLAST results
	my $blast_deleted = 0;
	my $extract_deleted = 0;
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

	# Update BLAST table with consolidated hits
	my $blast_updated = 0;
	my @ids = keys %$consolidated_ref;
	foreach my $record_id (@ids)  {

		my $where   = " WHERE record_id = $record_id ";
		my $hit_ref = $consolidated_ref->{$record_id};		
		delete $hit_ref->{record_id};
		#print "\n\t ## updating hit record ID $record_id in BLAST results";
		$blast_results_table->update($hit_ref, $where);
		$blast_updated++;
		
		# Delete from Extracted (because we need to re-extract and assign)
		#print "\n\t ## deleting hit BLAST ID $record_id in Extracted";
		my $ex_where = " WHERE blast_id = $record_id ";
		$extracted_table->delete_rows($ex_where);
		$extract_deleted++;
	}
}

############################################################################
# SECTION: Defragment/Consolidation
############################################################################

#***************************************************************************
# Subroutine:  defragment
# Description:  
#***************************************************************************
sub defragment {

	my ($self) = @_;
	print "\n\t ### Consolidating hits into Loci\n";
    # Set up for consolidation
    my @scaffolds;
    my @assigned_names;
   	my %chunks;
    $self->set_up_consolidation(\@scaffolds, \@assigned_names, \%chunks);
    ####my $key = $version . '/' . $chunk_name . '|' . $scaffold;
    ####my $content = $organism . '/' . $data_type;
    ####$chunks_ref->{$key} = $content;
    ####push (@$scaffs_ref, $key2);
	#$devtools->print_hash(\%chunks); die;   
 
    $self->get_tax_group();
    
	#$devtools->print_hash($self); die;   
    my $cons_obj = Consolidation->new($self);

    # DO THE CONSOLIDATION FOR EACH Assigned and Scaffold
        foreach my $assigned (@assigned_names){

        my @CONSmain;  # This is where the results get stored
            foreach my $scaf (@scaffolds){

                    # Run the consolidate on plus strand
                    push(@CONSmain, $cons_obj->consolidate( $scaf, '+', $assigned));
                    print "\nDone Consolidating $assigned in $scaf + orientation\n";

                    # Run the consolidate on minus strand
                    push(@CONSmain, $cons_obj->consolidate( $scaf, '-', $assigned));
                    print "\nDone Consolidating $assigned in $scaf - orientation\n";
            }    
    
            # Insert the data
            my $db_name   = $self->{db_name};
            $cons_obj->insert_loci_data($db_name, \@CONSmain);
            print "\n##################################################################";
            #die;
    }
	print "\n\t ### Done Consolidating hits into Loci\n";
}

#***************************************************************************
# Subroutine:  set up consolidation
# Description: 
#***************************************************************************
sub set_up_consolidation {

        my ($self, $scaffs_ref, $assigned_ref, $chunks_ref) = @_;

        my $db_obj = $self->{db};
        unless ($db_obj) { die; }

        # GET A LIST OF UNIQUE SCAFFOLDS, associated with chunk info
        my $extracted_table = $db_obj->{extracted_table};
        my @fields = qw [ organism data_type version target_name scaffold ];
        my @data;
        $extracted_table->select_distinct(\@fields, \@data);
        foreach my $row (@data) {

                my $organism    = $row->{organism};
                my $data_type   = $row->{data_type};
                my $version     = $row->{version};
                my $chunk_name  = $row->{target_name};
                my $scaffold    = $row->{scaffold};
                ####### DBM CHANGE HERE!!!!!!!!!!
                chomp($organism);
                chomp($data_type);
                chomp($version);
                chomp($chunk_name);
                chomp($scaffold);

                my $key = $version . '/' . $chunk_name . '^' . $scaffold;
                my $content = $organism . '/' . $data_type . '/';
                my $key2 = $content . $key;
                $chunks_ref->{$key} = $content;

                push (@$scaffs_ref, $key2);
        }

    # GET A LIST OF ALL ASSIGNED_NAME FROM EXTRACTED
    @fields = ();
    @fields = qw [ assigned_name ];
    @data = ();
    $extracted_table->select_distinct(\@fields, \@data);
    foreach my $row (@data) {

        my $assigned_name = $row->{assigned_name};
        chomp ($assigned_name);
        push(@$assigned_ref, $assigned_name);
    }
   
}


#***************************************************************************
# Subroutine:  get_tax_group
# Description: get the taxonomic group of the screening targets
#***************************************************************************
sub get_tax_group{

    my ($self) = @_;

    #GET PATH TO CHUNK
    my $first_group;
    my $group;
    my $paths = $self->{target_paths};
	unless ($paths) {
		print "\npaths\t$paths\n";
		die;
	}
	my $first_path = shift(@$paths);
    my @path_bits = split(/\//,$first_path);
    if ($first_path =~ /^\//){
        $first_group = $path_bits[1];
    }else{
        $first_group = $path_bits[0];
    }
    foreach my $path (@$paths){
        @path_bits =();
        @path_bits = split(/\//,$path);
        if ($path =~ /^\//){
            $group = $path_bits[1];
        }else{
            $group = $path_bits[0];
        }
        if ($group ne $first_group){
            die "\n\tCannot consolidate from Screening using multiple taxonomic groups\n";
        }
    }
    if($group){
        $self->{group} = $group;
    }else{
        $self->{group} = $first_group;
    }
    #and then from extracted get Organism Data type and VErsion (including the target_name)
}

############################################################################
# SECTION: UCSC tracks subroutines
############################################################################

#***************************************************************************
# Subroutine:  UCSCtracks
# Description: handler for UCSC tracks
#***************************************************************************
sub UCSCtracks{

    my ($self,$mode) = @_;

    print "\n\t ### Formating results into BED file for UCSC genome browser\n";

	my %colors = (  LTR => '0,0,0',
        LEA => '96,96,96', 
        gag => '255,0,0',
        pro => '0,0,255',
        pol => '255,0,255',
        env => '0,255,0'
	);
	my @fields;
	my $track_name;
	my @data;
	my $db_obj = $self->{db};
    unless ($db_obj) { die; }

	if($mode == 1){
        my $extracted_table = $db_obj->{extracted_table};
		@fields = qw (Scaffold Assigned_name Assigned_gene Orientation Extract_start Extract_end);
		$extracted_table->select_distinct(\@fields, \@data);
		$track_name = $self->{ucsc_extracted};
	}
	else {
		my $loci_table = $db_obj->{loci_table};
		@fields = qw (Scaffold Assigned_name Genome_structure Orientation Extract_start Extract_end);
		$loci_table->select_distinct(\@fields, \@data);
		$track_name = $self->{ucsc_loci};
	}
	my $outfile = $self->{output_path} . '/' . $track_name . '.txt';

    my $out = "track name=\"$track_name\" description=\"$track_name\" visibility=3 itemRgb=On useScore=0";
    $out .= "\n";

	foreach my $row (@data) {
        my $chr = $row->{Scaffold};
        my $start = $row->{Extract_start};
        my $end = $row->{Extract_end};
		if($start == 0 or $end == 0){
			next;
		}
		my $gene;
		if($mode == 1){
			$gene=$row->{Assigned_gene};
		}else{
			$gene=$row->{Genome_structure};
		}
		my $feature = $row->{Assigned_name} . '_' . $gene;
        my $score = 100;
        my $strand;
        if($row->{Orientation} eq '+'){
           $strand = '+';
        }else{
           $strand = '-';
        }
        my $col = $colors{$gene};
        unless($col){
            $col = '255,128,0';
        }   
        my $len = $end - $start;
        $out .= "$chr\t$start\t$end\t$feature\t$score\t$strand\t$start\t$end\t$col\t1\t$len\t0\n";
    }   
    open(OUT, ">$outfile") || die "Cannot open out\n";
    print OUT "$out";
    close(OUT);
}

############################################################################
# SECTION: Extend the screening database
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
	#print "\n\t infile 'print_array(\@infile);
	#sleep 2;
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
		my $row_count++;
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
		#$devtools->print_hash(\%insert);
		$anc_table->insert_row(\%insert);
	}
}

############################################################################
# SECTION: Utility subroutines
############################################################################

#***************************************************************************
# Subroutine:  run_utility_function
# Description: handler for utility functions
#***************************************************************************
sub run_utility_function {

	my ($self, $option, $ctl_file) = @_;

 	# Show title
	$self->show_title();  

	# Do initial set up and sanity checking for options that require it
	my $db_obj;
	if  ( $option < 4) {
		unless ($ctl_file) { die "\n\t Utility option '$option' requires a control file\n\n"; }
		$self->initialise($ctl_file);
		my $db_name = $self->{db_name};
		unless ($db_name) { die "\n\t Error: no DB name defined \n\n\n"; }
		$db_obj = ScreeningDB->new($self);
		$db_obj->load_screening_db($db_name);	
		$self->{db} = $db_obj; # Store the database object reference 
	}	
	
	# Hand off to functions based on the 
	if ($option eq 1) {  # Run database integrity checks
		print "\t #### SCREENING DATABASE SUMMARY\n";
		$db_obj->summarise_db();     # Do summary
		print "\t #### SCREENING DATABASE VALIDATION REPORT\n";
		$self->validate_screening_db(); # Find orphaned rows
	}
	elsif ($option eq 2) { # Remove a specific search
		$db_obj->clean_up(); # Data removal function
	}
	elsif ($option eq 3) { 
		$self->repair(); # Run a repair screen
	}
	elsif ($option eq 4) { # Summarise target genome databases, file by file
		my $genome_obj = TargetDB->new($self);
		$genome_obj->summarise_genomes_long();    
	}
	elsif ($option eq 5) { # Summarise target genome databases, one line per unique database
		my $genome_obj = TargetDB->new($self);
		$genome_obj->summarise_genomes_short();    
	}
	elsif ($option eq 6) { # Retrieve data
		my $refseqlib_obj = RefSeqLibrary->new($self);
		$refseqlib_obj->{reference_glue} = $ctl_file;
		print " UNIMPLEMENTED \n\n\n"; exit;
		$refseqlib_obj->convert_reference_library();    
	}
}

#***************************************************************************
# Subroutine:  validate_screening_db
# Description: validate screening database 
#***************************************************************************
sub validate_screening_db {

	my ($self) = @_;

	# Index data in the BLAST results and Extracted tables	
	my $db_obj = $self->{db};
	my %blast_results;
	my %extracted;
	my %status;
	$db_obj->index_BLAST_results_by_record_id(\%blast_results);
	$db_obj->index_extracted_loci_by_blast_id(\%extracted);
	$db_obj->index_previously_executed_queries(\%status);

	# Find Extracted table rows without BLAST table rows (should never happen)
	my @extracted_orphans;
	$db_obj->find_extracted_orphans(\%blast_results, \%extracted, \@extracted_orphans);
	my $num_extracted_orphans = scalar @extracted_orphans;
	if ($num_extracted_orphans) {
		print "\n\t Identified $num_extracted_orphans orphaned Extracted rows ";
	}
	else {
		print "\n\t There were no orphaned rows in the Extracted table ";
	}
	
	# Find BLAST table rows without Status table rows (should never happen)
	my @blast_orphans;
	$db_obj->find_blast_orphans(\%blast_results, \%status, \@blast_orphans);
	my $num_blast_orphans = scalar @blast_orphans;
	if ($num_blast_orphans) {
		print "\n\t Identified $num_blast_orphans BLAST rows with no search entry in the Status table ";
	}
	else {
		print "\n\t There were no orphaned rows in the BLAST_results table ";
	}

	# Find BLAST table rows without Extracted rows (i.e. unfinished round of paired BLAST)
	my @unfinished_searches;
	$db_obj->find_unfinished_searches(\%blast_results, \%extracted, \@unfinished_searches);
	my $num_unfinished_searches = scalar @unfinished_searches;
	if ($num_unfinished_searches) {
		print "\n\t Identified $num_unfinished_searches BLAST rows with no entry in the Extracted table ";
	}
	else {
		print "\n\t There were no half-finished screens ";
	}

	# Find Status table rows without BLAST hits (expected to happen sometimes)
	my @empty_searches;
	my $total_searches =  $db_obj->find_empty_searches(\%blast_results, \%status, \@empty_searches);
	my $num_empty_searches = scalar @empty_searches;
	if ($num_empty_searches) {
		my $percentage = ( $num_empty_searches / $total_searches ) * 100; 
		my $f_percentage = sprintf("%.2f", $percentage);
		print "\n\t $num_empty_searches of $total_searches (%$f_percentage) of BLAST searches returned no hits ";
		#$devtools->print_array(\@empty_searches);
	}
	else {
		print "\n\t All searches produced at least one hit ";
	}

	#$devtools->print_array(\@blast_orphans);
	#$devtools->print_array(\@extracted_orphans); die;
	# TODO  Go to console dialogue offering option to repair and/or repeat 
}

#***************************************************************************
# Subroutine:  repair 
# Description: run screens selectively to 'repair' a database
#***************************************************************************
sub repair {

	my ($self) = @_;

	# Index data in the BLAST results and Extracted tables	
	my $db_obj = $self->{db};
	my %blast_results;
	my %status;
	$db_obj->index_BLAST_results_by_record_id(\%blast_results);
	$db_obj->index_previously_executed_queries(\%status);

	# Find Status table rows without BLAST hits 
	my @empty_searches;
	my @blast_orphans;
	$db_obj->check_search_status(\%blast_results, \%status, \@blast_orphans, \@empty_searches);
	my $num_blast_orphans = scalar @blast_orphans;
	if ($num_blast_orphans) {
		print "\n\t Identified $num_blast_orphans BLAST rows with no search entry in the Status table ";
	}
	else {
		print "\n\t There were no BLAST table rows orphaned from the Status table ";
	}
	my $num_empty_searches = scalar @empty_searches;
	if ($num_empty_searches) {
		print "\n\t A total of $num_empty_searches BLAST searches retrurned no hits ";
	}
	else {
		print "\n\t All searches produced at least one hit ";
	}
	#$devtools->print_array(\@empty_searches); die;
	#my $reextract_flag = 'true';
	#$self->reassign($reextract_flag, \@blast_orphans);
	#die;
}

############################################################################
# Command line console fxns
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
	my ($HELP)  = "\n\t  usage: $0 [option] -i=[control file]\n";
        $HELP  .= "\n\t -m=1  create a screening DB"; 
		$HELP  .= "\n\t -m=2  execute a round of paired BLAST screening"; 
		$HELP  .= "\n\t -m=3  defragment hits (create loci table)"; 
		$HELP  .= "\n\t -m=4  reassign sequences (e.g. after reference library update)"; 
		$HELP  .= "\n\t -m=5  flush a screening DB"; 
		$HELP  .= "\n\t -m=6  drop a screening DB"; 
		$HELP  .= "\n\t -m=7  retrieve data from a screening DB"; 
		$HELP  .= "\n\t -m=8  manage ancillary tables in screening database\n"; 
	
		$HELP  .= "\n\t -u=1  validate screening DB"; 
		$HELP  .= "\n\t -u=2  clean up screening DB"; 
		$HELP  .= "\n\t -u=3  repair screening DB"; 
		$HELP  .= "\n\t -u=4  summarise target genome directory (file by file - slow)"; 
		$HELP  .= "\n\t -u=5  summarise target genome directory (by unique database - fast)"; 
		$HELP  .= "\n\t -u=6  summarise a GLUE-formatted reference sequence library\n"; 

		$HELP  .= "\n\t -r    format genome folder before executing"; 
	print $HELP;
}

#***************************************************************************
# Subroutine:  show_utility_help_page
# Description: show help page information
#***************************************************************************
sub show_utility_help_page {

	my ($self) = @_;
	
	# Initialise usage statement to print if usage is incorrect
	my ($HELP)  = "\n\t  usage: $0 -u=[option] -i=[control file]\n";
		$HELP  .= "\n\t -u=1  validate screening DB"; 
		$HELP  .= "\n\t -u=2  clean up screening DB"; 
		$HELP  .= "\n\t -u=3  repair screening DB"; 
		$HELP  .= "\n\t -u=4  summarise target genome directory"; 
		$HELP  .= "\n\t -u=5  summarise a GLUE-formatted reference sequence library"; 
		$HELP  .= "\n\n";
	print $HELP;
}

############################################################################
# EOF
############################################################################
