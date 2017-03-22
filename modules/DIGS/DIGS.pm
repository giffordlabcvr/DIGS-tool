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

	# Declare empty data structures
	my %crossmatching;

	# Set member variables
	my $self = {
		
		# Global settings
		process_id             => $parameter_ref->{process_id},
		program_version        => $parameter_ref->{program_version},
		
		# Flags
		verbose                => $parameter_ref->{verbose},
		force                  => $parameter_ref->{force},
		
		# Member classes 
		blast_obj              => $parameter_ref->{blast_obj},

		# MySQL database connection parameters
		mysql_username         => $parameter_ref->{mysql_username}, 
		mysql_password         => $parameter_ref->{mysql_password},
		db_name                => '',   # Obtained from control file or user
		mysql_server           => '',   # Obtained from control file or user

		# Parameters for DIGS
		query_na_fasta         => '',   # Obtained from control file
		query_aa_fasta         => '',   # Obtained from control file
		aa_reference_library   => '',   # Obtained from control file
		na_reference_library   => '',   # Obtained from control file
		bitscore_min_tblastn   => '',   # Obtained from control file
		bitscore_min_blastn    => '',   # Obtained from control file
		seq_length_minimum     => '',   # Obtained from control file

		# Paths used in DIGS process
		genome_use_path        => $parameter_ref->{genome_use_path},
		output_path            => $parameter_ref->{output_path},
		reference_na_fasta     => '',   
		reference_aa_fasta     => '',   
		blast_threads          => '',   # Obtained from control file

		# Data structures
		crossmatching          => \%crossmatching,

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

	my ($self, $ctl_file, $option) = @_;

	# Initialise
	$self->show_title();  
	my $valid = $self->initialise($option, $ctl_file);

	if ($valid) {
	
		# Hand off to DIGS functions
		if ($option eq 1) { 
			
			# Check the target sequences are formatted for BLAST
			$self->prepare_target_files_for_blast();
		}
		elsif ($option eq 2) { 
	
			# Run a DIGS process
			$self->perform_digs();	
		}
		elsif ($option eq 3) { 
	
			# Reassign data in digs_results table
			$self->reassign();	
		}
		elsif ($option eq 4) {
		
			# Interactively defragment results 	
			$self->interactive_defragment();	
		}
		elsif ($option eq 5) { 
	
			# Combine digs_results into higher order locus structures
			$self->consolidate_loci();
		}
		else {
					
			# Show error
			print "\n\t  Unrecognized option '-m=$option'\n";
		}
	}

	# Show final summary and exit message
	$self->wrap_up($option);

}

############################################################################
# PRIMARY FUNCTIONS (TOP LEVEL)
############################################################################

#***************************************************************************
# Subroutine:  prepare_target_files_for_blast
# Description: create index files for all target databases
#***************************************************************************
sub prepare_target_files_for_blast {

	my ($self) = @_;

	# Format targets files for BLAST searching
	my $target_db_obj = TargetDB->new($self);
	$target_db_obj->format_targets_for_blast();
}

#***************************************************************************
# Subroutine:  perform_digs
# Description: do the core database-integrated genome screening processes
#***************************************************************************
sub perform_digs {

	my ($self, $mode) = @_;

	# Get handle for the 'searches_performed' table, updated in this loop
	my $db_ref              = $self->{db};
	my $searches_table      = $db_ref->{searches_table};

	# Iterate through the list of DIGS queries, dealing each in turn 
	# Each DIGS query constitutes a probe sequence and a target FASTA file
	my $queries_completed = 0;
	my $queries_ref = $self->{queries};
	unless ($queries_ref) { die; }   # Sanity checking
	my @probes = keys %$queries_ref; # Get the list of queries
	print "\n\t ### Starting database-integrated genome screening";
	foreach my $probe_name (@probes) {
		
		# Get the array of queries for this target file
		my $probe_queries = $queries_ref->{$probe_name};
		foreach my $query_ref (@$probe_queries) {  
	
			# Increment query count
			$queries_completed++;		
			$self->{queries_completed} = $queries_completed;

			# Do the 1st BLAST (probe vs target)
			$self->search_target_file_using_blast($query_ref);
		
			# For this target, create a non-redundant locus set
			my @new_hits;
			$self->compile_nonredundant_locus_set($query_ref, \@new_hits);
			
			# Extract newly identified or extended sequences
			my @extracted;
			my $target_path = $query_ref->{target_path};
			$self->extract_sequences_from_target_file($target_path, \@new_hits, \@extracted);	
			
			# Do the 2nd BLAST (hits from 1st BLAST vs reference library)
			$self->classify_sequences_using_blast(\@extracted, $query_ref);
			
			# Update tables in the screening database to reflect new information
			$self->update_db(\@extracted, 'digs_results_table', 1);
	
			# Update the searches_performed table, indicating search has completed
			$searches_table->insert_row($query_ref);
		
			# Show a status update in the console
			$self->show_digs_progress();			
		}	
	}
}

#***************************************************************************
# Subroutine:  reassign
# Description: classify sequences already in the digs_results_table 
#***************************************************************************
sub reassign {
	
	my ($self) = @_;

	# Get data structures, paths and flags from self
	my $blast_obj   = $self->{blast_obj};
	my $result_path = $self->{report_dir};
	my $verbose     = $self->{verbose};

	# Get the connection to the digs_results table (so we can update it)
	my $db          = $self->{db};
	my $digs_results_table = $db->{digs_results_table};
	unless ($digs_results_table) { die; }
	
	# Get the sequences to reassign 
	my $reassign_loci = $self->{reassign_loci};
	unless ($reassign_loci) { die; }
	my $num_to_reassign = scalar @$reassign_loci;
	print "\n\n\t  Reassigning $num_to_reassign hits in the digs_results table\n";

	# Iterate through the loci, doing the reassign process for each
	my $count = 0;
	foreach my $locus_ref (@$reassign_loci) {
		
		# Set the linking to the BLAST result table
		my $record_id       = $locus_ref->{record_id};	
		my $extract_start   = $locus_ref->{extract_start};
		my $extract_end     = $locus_ref->{extract_end};
		$locus_ref->{subject_start} = $extract_start;
		$locus_ref->{subject_end}   = $extract_end;
		delete $locus_ref->{extract_start};
		delete $locus_ref->{extract_end};
	
		# Execute the 'reverse' BLAST (2nd BLAST in a round of paired BLAST)	
		my $previous_assign = $locus_ref->{assigned_name};
		my $previous_gene   = $locus_ref->{assigned_gene};
		$self->classify_sequence_using_blast($locus_ref);
		
		$count++;
		if (($count % 100) eq 0) { print "\n\t  Checked $count rows"; }

		my $assigned_name = $locus_ref->{assigned_name};
		my $assigned_gene = $locus_ref->{assigned_gene};
		if ($assigned_name ne $previous_assign or  $assigned_gene ne $previous_gene) {
				
			if ($verbose) {  # Report the outcome
				print "\n\t\t      - reassigned: was previously '$previous_assign ($previous_gene)'";
			}
			
			# Update the matrix
			my $previous_key = $previous_assign . '_' . $previous_gene;
			my $assigned_key = $assigned_name . '_' . $assigned_gene;	
			$self->update_cross_matching($previous_key, $assigned_key);
				
			# Insert the data
			my $where = " WHERE record_id = $record_id ";
			delete $locus_ref->{record_id}; # Required to remove this
			delete $locus_ref->{organism};  # Update not required for this field
			$digs_results_table->update($locus_ref, $where);
		}
	}
	
	# Write out the cross-matching matrix
	$self->show_cross_matching();

	# Cleanup
	my $output_dir = $self->{report_dir};
	my $command1 = "rm -rf $output_dir";
	system $command1;
}

#***************************************************************************
# Subroutine:  interactive_defragment 
# Description: interactively defragment via the console
#***************************************************************************
sub interactive_defragment {

	my ($self) = @_;

	# Get a list of all the target files from the screening DB
	$self->{defragment_mode} = 'defragment';
	my $db = $self->{db};
	my $digs_results_table = $db->{digs_results_table};
	my @fields = qw [ organism target_datatype target_version target_name ];
	my @targets;
	$digs_results_table->select_distinct(\@fields, \@targets);

	# Settings for clustering
	my %settings;
	$settings{total_loci}     = '0';
	$settings{total_clusters} = '0';
	$settings{range}          = undef;
	$settings{reextract}      = 1;

	# Preview changes 
	my $choice = undef;
	do    { 
		$choice = $self->preview_defragment(\@targets, \%settings);
	}   until ($choice > 1);

    # Apply the changes	if option is chosen
	if    ($choice eq 2) { 
		$self->defragment_target_files(\@targets, \%settings);
	}
	elsif ($choice eq 3) { 
		print "\n"; exit;
	}
	else { die; } # Should never get here
}

#***************************************************************************
# Subroutine:  consolidate_loci
# Description: assemble digs_results rows into higher-order loci 
#***************************************************************************
sub consolidate_loci {

	my ($self) = @_;

	# Get the digs results sorted by scaffold and extract start
	my @sorted;
	$self->get_sorted_digs_results(\@sorted);
	my $total_loci = scalar @sorted;
	
	# Compose clusters of overlapping/adjacent BLAST hits and extracted loci
	print "\n\t  Consolidating assigned extracted sequences into loci";
	print "\n\t  $total_loci loci in the digs_results table prior to consolidation'";
	my $settings_ref = $self->{consolidate_settings};
	unless ($settings_ref) { die; }
	my %consolidated;
	$self->compose_clusters(\%consolidated, \@sorted, $settings_ref);
	
	# Check the output
	my @cluster_ids  = keys %consolidated;
	my $num_clusters = scalar @cluster_ids;
	if ($total_loci > $num_clusters) {
		my $range = $settings_ref->{range};
		print "\n\t  $num_clusters clusters of loci within '$range' bp of one another ";
	}
	
	# Update locus data based on consolidated results
	$self->derive_locus_table_from_clustered_digs_results(\%consolidated);
	
	# Return the number of clusters
	return $num_clusters;
}

############################################################################
# INTERNAL FUNCTIONS: MAIN DIGS LOOP
############################################################################

#***************************************************************************
# Subroutine:  search_target_file_using_blast
# Description: execute a similarity search and parse the results
#***************************************************************************
sub search_target_file_using_blast {
	
	my ($self, $query_ref) = @_;

	# Get relevant member variables and objects
	my $blast_obj    = $self->{blast_obj};
	my $tmp_path     = $self->{tmp_path};
	my $min_length   = $self->{seq_length_minimum};
	my $min_score    = $self->{bitscore_minimum};

	# Sanity checking
	unless ($min_length) { die; }
	unless ($min_score)  { die; }	
	unless ($blast_obj)  { die; } 
	unless ($tmp_path)   { die; } 

	# Get query details
	my $probe_id        = $query_ref->{probe_id};
	my $blast_alg       = $query_ref->{blast_alg};
	my $probe_name      = $query_ref->{probe_name};
	my $probe_gene      = $query_ref->{probe_gene};
	my $probe_type      = $query_ref->{probe_type};
	my $probe_path      = $query_ref->{probe_path};
	my $organism        = $query_ref->{organism};
	my $version         = $query_ref->{target_version};
	my $datatype        = $query_ref->{target_datatype};
	my $target_name     = $query_ref->{target_name};
	my $target_path     = $query_ref->{target_path};
	my $result_file     = $tmp_path . "/$probe_id" . "_$target_name.blast_result.tmp";
	unless ($probe_id and $blast_alg) { die; }

	# Do BLAST similarity search
	my $completed = $self->{queries_completed};	
	print "\n\n\t  $blast_alg: $completed: '$organism' ($version, $datatype)";
	print   "\n\t  target: '$target_name'";
	print   "\n\t  probe:  '$probe_id'";   
	$blast_obj->blast($blast_alg, $target_path, $probe_path, $result_file);
	# TODO: catch error from BLAST and don't update "Searches_performed" table	
	
	# Extract the results from tabular format BLAST output
	my @hits;
	$blast_obj->parse_tab_format_results($result_file, \@hits);
	my $rm_command = "rm $result_file";
	system $rm_command; # Remove the result file

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
		if ($min_score) { 
			# Skip sequences that have too low bit scores
			my $query_score = $hit_ref->{bitscore};
			if ($query_score < $min_score) { 
				unless ($skip) { # Don't count as a bit_score exclusion if already exclude via length
					$skip = 'true';
					$score_exclude_count++;
				}
			}
		}	
		unless ($skip) {		
			# Insert values into 'active_set' table
			$self->insert_row_in_active_set_table($query_ref, $hit_ref);
			$num_retained_hits++;			
		}
	} 

	# Show summary of BLAST results after filtering
	if ($score_exclude_count or $length_exclude_count) {
		print "\n\t\t # $num_retained_hits matches above threshold ";
		print "(excluded: $length_exclude_count < length; $score_exclude_count < bitscore)";
	}
	
	return $num_hits;
}

#***************************************************************************
# Subroutine:  insert_row_in_active_set_table
# Description: insert a BLAST result as a row into the active set table
#***************************************************************************
sub insert_row_in_active_set_table {

	my ($self, $query_ref, $hit_ref) = @_;

	# Get screening database table objects
	my $db_ref           = $self->{db};
	my $active_set_table = $db_ref->{active_set_table};
	unless ($db_ref)          { die; } 

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

	$hit_ref->{digs_result_id}  = 0;
	$hit_ref->{organism}        = $organism;
	$hit_ref->{target_version}  = $version;
	$hit_ref->{target_datatype} = $datatype;
	$hit_ref->{target_name}     = $target_name;
	$hit_ref->{probe_id}        = $probe_id;
	$hit_ref->{probe_name}      = $probe_name;
	$hit_ref->{probe_gene}      = $probe_gene;
	$hit_ref->{probe_type}      = $probe_type;
	$hit_ref->{subject_start}   = $hit_ref->{aln_start};  # Rename to match DB
	$hit_ref->{subject_end}     = $hit_ref->{aln_stop};   # Rename to match DB
	$hit_ref->{query_end}       = $hit_ref->{query_stop}; # Rename to match DB
	$active_set_table->insert_row($hit_ref);

}

#***************************************************************************
# Subroutine:  extract_sequences_from_target_file
# Description: extract sequences from target databases
#***************************************************************************
sub extract_sequences_from_target_file {

	my ($self, $target_path, $loci_ref, $extracted_ref) = @_;

	# Get paths, objects, data structures and variables from self
	my $blast_obj = $self->{blast_obj};
	my $verbose   = $self->{verbose};
	my $buffer    = $self->{extract_buffer};

	# Iterate through the list of sequences to extract
	my $new_loci = 0;
	foreach my $locus_ref (@$loci_ref) {
			
		# Add any buffer 
		if ($buffer) { 
			my $orientation = $locus_ref->{orientation};
			$self->add_buffer_to_sequence($locus_ref, $orientation); 
		}
	
		# Extract the sequence
		my $sequence   = $blast_obj->extract_sequence($target_path, $locus_ref);
		if ($sequence) {
			
			# If we extracted a sequence, update the data for this locus
			my $seq_length = length $sequence; # Set sequence length
			if ($verbose) { print "\n\t\t    - Extracted sequence: $seq_length nucleotides "; }
			$locus_ref->{extract_start}   = $locus_ref->{start};
			$locus_ref->{extract_end}     = $locus_ref->{end};
			$locus_ref->{sequence}        = $sequence;
			$locus_ref->{sequence_length} = $seq_length;
			push (@$extracted_ref, $locus_ref);
			$new_loci++;
		}
		elsif ($verbose) { 
			print "\n\t\t    # Sequence extraction failed ";
		}
	}	
	return $new_loci;
}

#***************************************************************************
# Subroutine:  add_buffer_to_sequence
# Description: eadd leading-and-trailing buffer to extract coordinates
#***************************************************************************
sub add_buffer_to_sequence {

	my ($self, $hit_ref, $orientation) = @_;

	my $buffer = $self->{extract_buffer};
		
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

#***************************************************************************
# Subroutine:  classify_sequences_using_blast
# Description: classify a set of sequences using blast
#***************************************************************************
sub classify_sequences_using_blast {

	my ($self, $extracted_ref, $query_ref) = @_;

	#$devtools->print_array($extracted_ref); die;
	my $verbose = $self->{verbose};
	my $assigned_count   = 0;
	my $crossmatch_count = 0;
	unless ($query_ref) { die; }
	foreach my $locus_ref (@$extracted_ref) { # Iterate through the matches

		# Execute the 'reverse' BLAST (2nd BLAST in a round of paired BLAST)
		my $blast_alg = $self->classify_sequence_using_blast($locus_ref);
		my $assigned  = $locus_ref->{assigned_name};
		unless ($assigned) { die; }
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
		print "\n\t\t # $assigned_count extracted sequences classified";
	}
	if ($verbose) {	
		print "\n\t\t # $crossmatch_count cross-matched to something other than the probe";
	}
}

#***************************************************************************
# Subroutine:  classify_sequence_using_blast
# Description: classify a nucleotide sequence using blast 
#***************************************************************************
sub classify_sequence_using_blast {

	my ($self, $locus_ref) = @_;

	# Get paths and objects from self
	my $result_path = $self->{tmp_path};
	my $blast_obj   = $self->{blast_obj};
	my $verbose     = $self->{verbose};
	unless ($blast_obj)   { die; } 
	unless ($result_path) { die; }
	
	# Get required data about the query sequence
	my $sequence   = $locus_ref->{sequence};
	my $probe_type = $locus_ref->{probe_type};
	unless ($probe_type) { die; } # Sanity checking
	unless ($sequence)   { die; } # Sanity checking

	# Make a FASTA query file
	$sequence =~ s/-//g;   # Remove any gaps that might happen to be there
	$sequence =~ s/~//g;   # Remove any gaps that might happen to be there
	$sequence =~ s/\s+//g; # Remove any gaps that might happen to be there
	my $fasta      = ">TEMPORARY\n$sequence";
	my $query_file = $result_path . '/TEMPORARY.fas';
	$fileio->write_text_to_file($query_file, $fasta);
	my $result_file = $result_path . '/TEMPORARY.blast_result';
	
	# Do the BLAST according to the type of sequence (AA or NA)
	my $blast_alg = $self->get_blast_algorithm($probe_type);
	my $lib_path  = $self->get_blast_library_path($probe_type);
	my $lib_file;
	if    ($probe_type eq 'ORF') {  $lib_file = $self->{aa_reference_library}; }
	elsif ($probe_type eq 'UTR') {  $lib_file = $self->{na_reference_library}; }
	else  { die; }
	unless ($lib_file)  { die; }

	# Execute the call to BLAST and parse the results
	$blast_obj->blast($blast_alg, $lib_path, $query_file, $result_file);
	my @results;
	$blast_obj->parse_tab_format_results($result_file, \@results);

	# Define some variables for capturing the result
	my $top_match = shift @results;
	my $query_start   = $top_match->{query_start};
	my $query_end     = $top_match->{query_stop};
	my $subject_start = $top_match->{aln_start};
	my $subject_end   = $top_match->{aln_stop};
	my $assigned_key  = $top_match->{scaffold};	
	my $assigned;

	# Deal with a query that matched nothing in the 2nd BLAST search
	unless ($assigned_key) {
		$self->set_default_values_for_unassigned_locus($locus_ref);	
		$assigned = undef;
	}
	else {	# Assign the extracted sequence based on matches from 2nd BLAST search

		# Split assigned to into (i) refseq match (ii) refseq description (e.g. gene)	
		my @assigned_key  = split('_', $assigned_key);
		my $assigned_gene = pop @assigned_key;
		my $assigned_name = shift @assigned_key;
		#$assigned_name = join ('_', @assigned_name);
		$locus_ref->{assigned_name}  = $assigned_name;
		$locus_ref->{assigned_gene}  = $assigned_gene;
		$locus_ref->{identity}       = $top_match->{identity};
		$locus_ref->{bitscore}       = $top_match->{bitscore};
		$locus_ref->{evalue_exp}     = $top_match->{evalue_exp};
		$locus_ref->{evalue_num}     = $top_match->{evalue_num};
		$locus_ref->{mismatches}     = $top_match->{mismatches};
		$locus_ref->{align_len}      = $top_match->{align_len};
		$locus_ref->{gap_openings}   = $top_match->{gap_openings};
		$locus_ref->{query_end}      = $query_end;
		$locus_ref->{query_start}    = $query_start;
		$locus_ref->{subject_end}    = $subject_end;
		$locus_ref->{subject_start}  = $subject_start;
		if ($verbose) { 
			my $id = $locus_ref->{record_id};
			print "\n\t\t    - Record '$id' assigned as '$assigned_name ($assigned_gene)'";
		 	print " via $blast_alg comparison to $lib_file";
		 }
		$assigned = $assigned_name . '_' . $assigned_gene;
	}

	# Clean up
	my $command1 = "rm $query_file";
	my $command2 = "rm $result_file";
	system $command1;
	system $command2;
	
	return $blast_alg;
}

#***************************************************************************
# Subroutine:  set_default_values_for_unassigned_locus
# Description: set default values for an unassigned extracted sequence
#***************************************************************************
sub set_default_values_for_unassigned_locus {

	my ($self, $hit_ref) = @_;

	$hit_ref->{assigned_name}    = 'Unassigned';
	$hit_ref->{assigned_gene}    = 'Unassigned';
	$hit_ref->{identity}         = 0;
	$hit_ref->{bitscore}         = 0;
	$hit_ref->{evalue_exp}       = 0;
	$hit_ref->{evalue_num}       = 0;
	$hit_ref->{mismatches}       = 0;
	$hit_ref->{align_len}        = 0;
	$hit_ref->{gap_openings}     = 0;
	$hit_ref->{query_end}        = 0;
	$hit_ref->{query_start}      = 0;
	$hit_ref->{subject_end}      = 0;
	$hit_ref->{subject_start}    = 0;
	
}

#***************************************************************************
# Subroutine:  get_blast_algorithm
# Description: determine which blast algorithm to use based on settings
#***************************************************************************
sub get_blast_algorithm {

	my ($self, $probe_type) = @_;
	
	my $blast_alg;
	if    ($probe_type eq 'UTR') { $blast_alg = 'blastn'; }
	elsif ($probe_type eq 'ORF') { $blast_alg = 'blastx'; }
	else { die "\n\t Unknown probe type '$probe_type '\n\n"; }
	
	return $blast_alg;
}

#***************************************************************************
# Subroutine:  get_blast_library_path
# Description: get path to a reference library, based on settings
#***************************************************************************
sub get_blast_library_path {

	my ($self, $probe_type) = @_;
	my $lib_path;
	
	if ($probe_type eq 'UTR') { 
		$lib_path = $self->{blast_utr_lib_path};
		unless ($lib_path) {
			$devtools->print_hash($self); 
			die "\n\t NO UTR LIBRARY defined";
		}
	}
	elsif ($probe_type eq 'ORF') { 
		$lib_path = $self->{blast_orf_lib_path};
		unless ($lib_path) {
			die "\n\t NO ORF LIBRARY defined";
		}
	}	
	return $lib_path;
}

#***************************************************************************
# Subroutine:  update_db
# Description: update the screening DB based on a completed round of DIGS
#***************************************************************************
sub update_db {

	my ($self, $extracted_ref, $table_name, $update) = @_;
		
	# Get parameters from self
	my $db_ref              = $self->{db};
	my $verbose             = $self->{verbose};
	my $digs_results_table  = $db_ref->{$table_name}; 
	my $active_set_table    = $db_ref->{active_set_table}; 
	my $blast_chains_table  = $db_ref->{blast_chains_table}; 

	# Iterate through the extracted sequences
	my $deleted = '0';
	foreach my $hit_ref (@$extracted_ref) {
		
		# Insert the data to the digs_results table
		my $digs_result_id;
		if ($update) {
			$digs_result_id = $digs_results_table->insert_row($hit_ref);
		}
		else {
			$digs_result_id = $hit_ref->{digs_result_id};
			my $where = " WHERE record_id = $digs_result_id ";
			my %update;
			$update{extract_start} = $hit_ref->{extract_start};
			$update{extract_end}   = $hit_ref->{extract_end};
			$digs_results_table->update(\%update, $where);		
		}
		
		# Insert the data to the BLAST_chains table
		my $blast_chains = $hit_ref->{blast_chains};
		if ($blast_chains) {		
			my @blast_ids = keys %$blast_chains;
			foreach my $blast_id (@blast_ids) {							
				my $data_ref = $blast_chains->{$blast_id};
				$data_ref->{digs_result_id} = $digs_result_id;	
				$blast_chains_table->insert_row($data_ref);
			}
		}
		unless ($digs_result_id) { die; }
		
		# Delete superfluous data from the digs_results table
		my $digs_result_ids_ref = $hit_ref->{digs_result_ids};
		foreach my $old_digs_result_id (@$digs_result_ids_ref) {			
			
			# Where we updated an existing record, keep that record
			unless ($old_digs_result_id eq $digs_result_id) {

				# Delete superfluous extract rows
				my $extracted_where = " WHERE record_id = $old_digs_result_id ";	
				if ($verbose) { print "\n\t\t    - Deleting redundant locus '$old_digs_result_id'"; }
				$digs_results_table->delete_rows($extracted_where);
				$deleted++;

				# Update extract IDs			
				my $chains_where = " WHERE record_id = $old_digs_result_id ";
				my %new_id;
				$new_id{digs_result_id} = $digs_result_id;	
				$blast_chains_table->update(\%new_id, $chains_where);
			}
		}
	}

	# Flush the active set table
	$active_set_table->flush();

	# Return the number
	return $deleted;
}

#***************************************************************************
# Subroutine:  show_digs_progress
# Description: show progress in DIGS screening
#***************************************************************************
sub show_digs_progress {

	my ($self) = @_;

	# Get the counts
	my $total_queries   = $self->{total_queries};
	my $completed       = $self->{queries_completed};	
	unless ($completed and $total_queries) { die; } # Sanity checking
	
	# Calculate percentage progress
	my $percent_prog    = ($completed / $total_queries) * 100;
	my $f_percent_prog  = sprintf("%.2f", $percent_prog);
	#print "\n\t\t  ";
	print "\n\t\t # done $completed of $total_queries queries (%$f_percent_prog)";
}

#***************************************************************************
# Subroutine:  wrap_up
# Description: clean-up functions etc prior to exiting program
#***************************************************************************
sub wrap_up {

	my ($self, $option) = @_;

	# Remove the output directory
	my $output_dir = $self->{report_dir};
	if ($output_dir) {
		my $command1 = "rm -rf $output_dir";
		system $command1;
	}
	
	# Show cross matching at end if verbose output setting is on
	my $verbose = $self->{verbose};
	if ($verbose and $option eq 2 and $option eq 3) { 
		$self->show_cross_matching();
	}

	# Print finished message
	print "\n\n\t ### Process completed ~ + ~ + ~";

}

############################################################################
# INTERNAL FUNCTIONS: CONSOLIDATE LOCI
############################################################################

#***************************************************************************
# Subroutine:  derive_locus_table_from_clustered_digs_results
# Description: compile locus information  and update the locus tables
#***************************************************************************
sub derive_locus_table_from_clustered_digs_results {

	my ($self, $consolidated_ref) = @_;

	unless ($consolidated_ref) { die; }
	
	# Get parameters and data structures
	my $verbose           = $self->{verbose};
	my $db_ref            = $self->{db};
	my $loci_table        = $db_ref->{loci_table};
	my $loci_chains_table = $db_ref->{loci_chains_table};
	
	# Flags for how to handle
	#my $reextract = undef;
	my $reextract = 'true';
	#my $annotate_ends = undef;
	my $annotate_ends = 'true';

	# Iterate through the clusters	
	my @cluster_ids  = keys %$consolidated_ref;
	foreach my $cluster_id (@cluster_ids) {

		# Get the loci in this cluster
		my $cluster_ref = $consolidated_ref->{$cluster_id};

		# Turn this cluster into an annotated locus
		my %locus;
		$self->derive_locus_structure(\%locus, $cluster_ref);
		#$devtools->print_hash(\%locus); die;

		# Extract the consolidate locus if the flag is set
		if ($reextract) {
			$self->extract_consolidated_locus(\%locus);
		}

		# Do the annotation for truncated versus non-truncated 	matches
		if ($annotate_ends) {
			#$self->annotate_consolidated_locus_flanks(\%locus);
		}

		# Insert the consolidated locus information
		my $locus_array = $locus{locus_array};
		my $locus_structure = join('-', @$locus_array);
		$locus{locus_structure} = $locus_structure;
	
		# Insert the data	
		my $locus_id  = $loci_table->insert_row(\%locus);
				
		# Create the links between the loci and digs_results tables
		foreach my $digs_result_ref (@$cluster_ref) {
			my $digs_result_id = $digs_result_ref->{record_id};
			my %chain_data;
			$chain_data{digs_result_id} = $digs_result_id;
			$chain_data{locus_id}       = $locus_id;
			$loci_chains_table->insert_row(\%chain_data);
		}		
	}
}

#***************************************************************************
# Subroutine:  derive_locus_structure
# Description: derive locus structure based on clustered digs results 
#***************************************************************************
sub derive_locus_structure {

	my ($self, $consolidated_ref, $cluster_ref) = @_;

	my $annotate_flanks = undef;
	my $initialised = undef;
	my $organism;
	my $version;
	my $target_name;
	my $datatype;
	my $assigned_name;
	my $feature;
	my $lowest;
	my $highest;
	my $scaffold;
	my $orientation;
	my $target_datatype;
	my $target_version;
	my @locus_structure;
	my $target_id;
	my $multiple_orientations = undef;
	my $last_element = undef;
	foreach my $element_ref (@$cluster_ref) {
		
		# Capture values from the previous iterations	
		my $last_feature     = $feature;
		my $last_orientation = $orientation;
		#my $last_scaffold    = $scaffold;
		
		# Get the data about this digs_results table row
		my $start       = $element_ref->{extract_start};
		my $end         = $element_ref->{extract_end};

		$feature        = $element_ref->{assigned_gene};
		$assigned_name  = $element_ref->{assigned_name};
		$version        = $element_ref->{target_version};
		$datatype       = $element_ref->{target_datatype};
		$orientation    = $element_ref->{orientation};
		$scaffold       = $element_ref->{scaffold};
		$organism       = $element_ref->{organism};
		$target_name    = $element_ref->{target_name};
		unless ($feature and $orientation) { die; } # Sanity checking
		my $record      = "$feature($orientation)";
		#print "\n\t RECORD $record";

		# Create a target key
		$organism        = $element_ref->{organism};
		$target_name     = $element_ref->{target_name};
		$target_datatype = $element_ref->{target_datatype};
		$target_version  = $element_ref->{target_version};
		my @genome = ( $organism , $target_datatype, $target_version );
		my $this_target_id = join ('|', @genome);
		if ($target_id) {
			unless ($this_target_id eq $target_id) { 
				print "\n\t WHY??? $target_id NE $this_target_id\n\n";
				#die; 
			} 
		}
		$target_id = $this_target_id;

		
		# Deal with first locus in a cluster
		unless ($initialised) {
			$highest      = $end;
			$lowest       = $start;
			$last_element = $element_ref;
			$initialised  = 'true';								
			push(@locus_structure, $record);
			next;		
		}


		# Capture information about coordinates 			
		if ($end > $highest) {
			$highest = $end;
		}
		if ($start < $lowest) {
			$lowest = $start;					
		}

		# Deal with loci that follow at least one previous locus
		if ($orientation eq $last_orientation
		and $feature eq $last_feature) {
			next;
		}
		if ($orientation ne $last_orientation) {
			$multiple_orientations = 'true';
		}

		# Add this element to the start or end of the locus array, based on orientation
		if ($multiple_orientations) { # If multiple orientations, base it on last element
			if ($start >= $last_element->{extract_start}) {
				push(@locus_structure, $record);
			}
			else {
				unshift(@locus_structure, $record);	
			}
		}
		elsif ($orientation eq '+') {
			#my $record = $feature;
			push(@locus_structure, $record);
		}
		elsif ($orientation eq '-') {
			unshift(@locus_structure, $record);			
		}
		$last_element = $element_ref;				
	}

	# Store the data
	$consolidated_ref->{organism}        = $organism;
	$consolidated_ref->{target_version}  = $version;
	$consolidated_ref->{target_name}     = $target_name;
	$consolidated_ref->{target_datatype} = $datatype;
	$consolidated_ref->{scaffold}        = $scaffold;
	$consolidated_ref->{target_id}       = $target_id;
	$consolidated_ref->{orientation}     = $orientation;
	$consolidated_ref->{start}           = $lowest;
	$consolidated_ref->{end}             = $highest;
	$consolidated_ref->{extract_start}   = $lowest;
	$consolidated_ref->{extract_end}     = $highest;
	$consolidated_ref->{assigned_name}   = $assigned_name;
	$consolidated_ref->{assigned_name}   = $assigned_name;
	$consolidated_ref->{locus_array}     = \@locus_structure;
	
}

#***************************************************************************
# Subroutine:  extract_consolidated_locus
# Description: 
#***************************************************************************
sub extract_consolidated_locus {

	my ($self, $consolidated_ref) = @_;

	my $db_ref    = $self->{db};
	my $verbose   = $self->{verbose};
	my $blast_obj = $self->{blast_obj};
	my $seq_len   = 0;
	
	my $genome_use_path  = $self->{genome_use_path};
	my $target_group_ref = $self->{target_groups};
	#my $target_path = $self->get_target_file_path($target_ref);

	my $organism        = $consolidated_ref->{organism};
	my $target_version  = $consolidated_ref->{target_version};
	my $target_datatype = $consolidated_ref->{target_datatype};
	my $target_name     = $consolidated_ref->{target_name};
	my $target_id       = $consolidated_ref->{target_id};
	my $lowest          = $consolidated_ref->{start};
	my $highest         = $consolidated_ref->{end};

	my $full_id = $target_id . '|' . $target_name;
	my $target_group = $target_group_ref->{$full_id};
	unless ($target_group) {
		print " \n\t No target group found for TARGET ID $full_id\n\n"; 
		#$devtools->print_hash($target_group_ref);
        sleep 1;
		return 0;
	}
	
	# Construct the path to this target file
	my @path;
	push (@path, $genome_use_path);
	push (@path, $target_group);
	push (@path, $organism);
	push (@path, $target_datatype);
	push (@path, $target_version);
	push (@path, $target_name);
	my $target_path = join ('/', @path);

	# Extract the sequence
	#print "\n\t\t    # TARGET: '$target_path'";
	my $sequence   = $blast_obj->extract_sequence($target_path, $consolidated_ref);
	my $seq_length = length $sequence; # Set sequence length
	if ($sequence) {
		
		# If we extracted a sequence, update the data for this locus
		if ($verbose) { print "\n\t\t    - Re-extracted sequence: $seq_length nucleotides "; }
		$consolidated_ref->{sequence}        = $sequence;
		$consolidated_ref->{sequence_length} = $seq_length;
	}
	elsif ($verbose) { 
		print "\n\t\t    # Sequence extraction failed ";
	}
}

#***************************************************************************
# Subroutine:  annotate_consolidated_locus_flanks
# Description: 
#***************************************************************************
sub annotate_consolidated_locus_flanks {

	my ($self, $consolidated_ref) = @_;

	my $db_ref = $self->{db};
	my $contigs_table = $db_ref->{contigs_table};
	my $lowest   = $consolidated_ref->{start};
	my $highest  = $consolidated_ref->{end};
	my $scaffold = $consolidated_ref->{scaffold};

	# Get the length of this contig
	my %data;
	my @fields = qw [ contig_id seq_length ];
	my $where = " WHERE contig_id = '$scaffold'";
	$contigs_table->select_row(\@fields, \%data, $where);
	my $contig_length = $data{seq_length};
	unless ($contig_length) { die; }

	# Check the start of the match
	my @locus_structure;
	if ($lowest eq 1) {  
		unshift(@locus_structure, 'T');
	}
	else {
		unshift(@locus_structure, 'X');
	}
	# Check the end of the match
	if ($highest eq $contig_length) { 
		push(@locus_structure, 'T');
	}
	else {
		push(@locus_structure, 'X');
	}
}

############################################################################
# INTERNAL FUNCTIONS: defragmenting results
############################################################################

#***************************************************************************
# Subroutine:  compile_nonredundant_locus_set
# Description: determine what to extract based on current results
# Note similar to defragment_target fxn
#***************************************************************************
sub compile_nonredundant_locus_set {
	
	my ($self, $query_ref, $to_extract_ref) = @_;

	# Get flag
	my $verbose         = $self->{verbose};

	# Compose SQL WHERE statement to retrieve relevant set of loci
	my $target_name     = $query_ref->{target_name};
	my $organism        = $query_ref->{organism};
	my $probe_name      = $query_ref->{probe_name};
	my $probe_gene      = $query_ref->{probe_gene};
	my $probe_type      = $query_ref->{probe_type};

	my $where  = " WHERE organism = '$organism' ";
	   $where .= " AND target_name = '$target_name' "; # Always limit by target
	   $where .= " AND probe_type  = '$probe_type' ";  # EITHER utrs OR orfs NOT BOTH 

	# Get the relevant set of DIGS results
	my @digs_results;
	$self->get_sorted_digs_results(\@digs_results, $where);
	my $num_loci = scalar @digs_results;
	if ($verbose) { print "\n\t\t # $num_loci previously extracted $probe_type loci"; }
		
	# Add the digs results to the BLAST hits in the active_set table
	$self->add_digs_results_to_active_set(\@digs_results);

	# Get sorted list of digs results and BLAST hits from active_set table
	my @combined;
	$self->get_sorted_active_set(\@combined);
	my $total_loci = scalar @combined;
	if ($verbose) {
		if ($total_loci > 0) {
			print "\n\t\t # $total_loci rows in active set (including $num_loci previously extracted) ";
		}
	}
	
	# Compose clusters of overlapping/adjacent loci
	my %settings;
	my %defragmented;
	$settings{range} = $self->{defragment_range};
	$settings{start} = 'subject_start';
	$settings{end}   = 'subject_end';
	$self->compose_clusters(\%defragmented, \@combined, \%settings);
	# DEBUG $self->show_clusters(\%defragmented);  # Show clusters
	my @cluster_ids  = keys %defragmented;
	my $num_clusters = scalar @combined;
	if ($verbose) {
		if ($total_loci > $num_clusters) {
		}
	}

	# Get a resolved list of non-redundant, non-overlapping loci to extract
	$self->merge_clustered_loci(\%defragmented, $to_extract_ref);
	my $num_new = scalar @$to_extract_ref;
	if ($num_new){
		print "\n\t\t # $num_new sequences to extract";
	}	
	else {
		print "\n\t\t # No new loci to extract";
	}	
}

#***************************************************************************
# Subroutine:  preview_defragment
# Description: preview a defragmentation process
#***************************************************************************
sub preview_defragment {

	my ($self, $targets_ref, $settings_ref) = @_;

	# Display current settings	
	my $defragment_range = $self->{defragment_range};
	unless ($defragment_range )  { die; } 
	print "\n\n\t\t Current settings (based on control file)";
	print "\n\t\t defragment range: $defragment_range";

	# Get the range	
	my $question1 = "\n\n\t # Set the range for merging hits";
	my $t_range = $console->ask_int_with_bounds_question($question1, $defragment_range, $maximum);		
	#my $t_range = 1000;

	# Preview this defragment
	$self->defragment_digs_results($targets_ref, $settings_ref, $t_range);

	# Summarise results
	my $total_loci     = $settings_ref->{total_loci};
	my $total_clusters = $settings_ref->{total_clusters};
	print "\n\t\t\t TOTAL LOCI:     $total_loci";
	print "\n\t\t\t TOTAL CLUSTERS: $total_clusters ";
	print "\n\n\t\t Option 1: preview new parameters";
	print "\n\t\t Option 2: apply these parameters";
	print "\n\t\t Option 3: exit";

	# Prompt for what to do next
	my $list_question = "\n\n\t # Choose an option:";
	my $choice = $console->ask_list_question($list_question, 3);
	$settings_ref->{range} = $t_range;
	#my $choice = 2;

	return $choice;
}

#***************************************************************************
# Subroutine:  defragment_digs_results
# Description: preview results of a defragment process (for interactive defragment)
#***************************************************************************
sub defragment_digs_results {

    my ($self, $targets_ref, $cluster_params, $t_range) = @_;
   
	# Apply the settings
	my $verbose = $self->{verbose};
	my $total_loci = '0';
	my $total_clusters = '0';
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
		
		# Get number of clusters
		my @cluster_ids  = keys %target_defragmented;
		my $num_clusters = scalar @cluster_ids;

		# Show clusters if verbose flag is set
		if ($verbose) { print "\n\n\t\t Interactive defrag: $num_hits hits in target $target_name"; }
		$total_loci     = $total_loci + $num_hits;
		$total_clusters = $total_clusters + $num_clusters;
	}

	$cluster_params->{total_loci}     = $total_loci;
	$cluster_params->{total_clusters} = $total_clusters;
		
}

#***************************************************************************
# Subroutine:  defragment_target_files
# Description: implement a defragmentation process for a set of target files
#***************************************************************************
sub defragment_target_files {

	my ($self, $targets_ref,  $settings_ref) = @_;

	my $genome_use_path  = $self->{genome_use_path};
	my $target_group_ref = $self->{target_groups};
	my $db               = $self->{db};
	my $t_range          = $settings_ref->{range};
	my $reextract        = $settings_ref->{reextract};
	if ($reextract) {
		unless ($genome_use_path and $target_group_ref)  { die; }
	}

	# Create a copy of the digs_results table (changes will be applied to copy)
	# TODO: fix this
	print "\n\t # Defragmenting using range '$t_range'\n";
	my $copy_name = $db->backup_digs_results_table();
	print "\n\t # Copied DIGS results to '$copy_name'\n";
	my $dbh = $db->{dbh};
	$db->load_digs_results_table($dbh, 'digs_results');	
	unless ($db->{digs_results_table}) { die; }
	
	# Iterate through the target files, applying the defragment process to each		
	foreach my $target_ref (@$targets_ref) {

		# Get the target details (and thus the target path)
		#my $target_path = $self->get_target_file_path($target_ref);
		my $organism        = $target_ref->{organism};
		my $target_datatype = $target_ref->{target_datatype};
		my $target_version  = $target_ref->{target_version};
		my $target_name     = $target_ref->{target_name};
		my @genome = ( $organism , $target_datatype, $target_version, $target_name );
		my $target_id       = join ('|', @genome);
		print "\n\t\t # Defragmenting hits in '$target_name'";

		my $target_path = 'NULL';
		if ($reextract) {
			my $target_group = $target_group_ref->{$target_id};
			unless ($target_group) { 
				$devtools->print_hash($target_group_ref);
				print "\n\t Didn't get target group name for target file with id '$target_id'\n\n";
				die; 
			}
			# Construct the path to this target file
			my @path;
			push (@path, $genome_use_path);
			push (@path, $target_group);
			push (@path, $organism);
			push (@path, $target_datatype);
			push (@path, $target_version);
			push (@path, $target_name);
			$target_path = join ('/', @path);
		}

		# Construct WHERE statement
		my $where  = " WHERE organism      = '$organism' ";
		$where    .= " AND target_datatype = '$target_datatype' ";
		$where    .= " AND target_version  = '$target_version' ";
		$where    .= " AND target_name     = '$target_name' "; 
		$settings_ref->{start}     = 'extract_start';
		$settings_ref->{end}       = 'extract_end';
		$settings_ref->{where_sql} = $where;
		$self->defragment_target($settings_ref, $target_path, 'digs_results');
	}
}

#***************************************************************************
# Subroutine:  defragment_target 
# Description: implement a defragmentation process for a single target file
# Note: similar to compile_nonredundant_locus_set fxn, diff details
#***************************************************************************
sub defragment_target {

	my ($self, $settings_ref, $target_path, $copy_name) = @_;
	
	# Create the relevant set of previously extracted loci
	my @combined;
	my %target_defragmented;
	my $where     = $settings_ref->{where_sql};
    my $copy_table_name = $copy_name . '_table';

	# Get digs results
	$self->get_sorted_digs_results(\@combined, $where);
	my $num_hits = scalar @combined;
	print "\n\t\t # $num_hits digs results to defragment ";
	#$devtools->print_array(\@combined); die;
		
	# Compose clusters of overlapping/adjacent BLAST hits and extracted loci
	$self->compose_clusters(\%target_defragmented, \@combined, $settings_ref);
	my @cluster_ids  = keys %target_defragmented;
	my $num_clusters = scalar @combined;
	if ($num_clusters < $num_hits) {
		#$self->show_clusters(\%target_defragmented);  # Show clusters
		my $range = $settings_ref->{range};
		print "...compressed to $num_clusters overlapping/contiguous clusters within '$range' bp of one another";
	}
	
	# Determine what to extract, and extract it
	my @loci;
	my $reextract = $settings_ref->{reextract};
	my $extended  = $self->merge_clustered_loci(\%target_defragmented, \@loci, $reextract);
	print "\n\t\t # $extended extensions to previously extracted sequences ";
	my $num_new   = scalar @loci;
	print "\n\t\t # $num_new loci to extract after defragment ";

	if ($reextract and $num_new) {

		# Extract newly identified or extended sequences
		my @extracted;
		$self->extract_sequences_from_target_file($target_path, \@loci, \@extracted);
		
		# Do the genotyping step for the newly extracted locus sequences
		my $assigned_count   = 0;
		my $crossmatch_count = 0;
		my $num_extracted = scalar @extracted;
		print "\n\t\t # Genotyping $num_extracted newly extracted sequences:";
		foreach my $hit_ref (@extracted) { # Iterate through loci		
			$self->classify_sequence_using_blast($hit_ref);
			$assigned_count++;
			my $remainder = $assigned_count % 100;
			if ($remainder eq 0) { print "\n\t\t\t # $assigned_count sequences classified "; }
		}
		print "\n\t\t\t # $assigned_count sequences classified";

		# Update DB
		my $num_deleted = $self->update_db(\@extracted, $copy_table_name, 1);
		print "\n\t\t\t # $num_deleted rows deleted from digs_results table\n";
	}
	else { # DEBUG
		# Update DB
		$self->prepare_locus_update(\@loci);
		$self->update_db(\@loci, $copy_table_name);
	}
	return $num_new;
}

############################################################################
# INTERNAL FUNCTIONS: clustering/merging overlapping/adjacent loci
############################################################################

#***************************************************************************
# Subroutine:  compose_clusters 
# Description: process a sorted list of loci and group into 'clusters' of
#              overlapping feature annotations
#***************************************************************************
sub compose_clusters {

	my ($self, $defragmented_ref, $loci_ref, $settings_ref) = @_;
	
	# Get settings
	my $start_token = $settings_ref->{start};
	my $end_token   = $settings_ref->{end};
	
	# Iterate through loci, grouping them into clusters when they are within range
	my $j = 1;
	my %last_locus;
	my %name_counts;
	my $initialised = undef;
	foreach my $locus_ref (@$loci_ref)  {

		# Get locus data
		my $record_id     = $locus_ref->{record_id};
		my $scaffold      = $locus_ref->{scaffold};
		my $target_name   = $locus_ref->{target_name};
		my $assigned_name = $locus_ref->{assigned_name};
		my $probe_gene    = $locus_ref->{probe_gene};
		my $assigned_gene = $locus_ref->{assigned_gene};
		my $orientation   = $locus_ref->{orientation};
		my $start         = $locus_ref->{$start_token};
		my $end           = $locus_ref->{$end_token};
		
		# Get last hit values
		my $last_scaffold      = $last_locus{scaffold};
		my $last_start         = $last_locus{$start_token};
	    if ($initialised) {

			# Sanity checking - are sequences in sorted order for this scaffold?
			if ( $scaffold eq $last_scaffold and $start < $last_start) { die; }
           
            # Work out wether to merge this hit with the last
            my $merge = $self->compare_adjacent_loci($locus_ref, \%last_locus, $settings_ref);
            
            unless ($merge) {
                 # Increment the count
                $j++;       
                # Initialise record
                $self->initialise_cluster($defragmented_ref, $locus_ref, $j);
            }
            else {             
                # Extend record
                $self->extend_cluster($defragmented_ref, $locus_ref, $j);
            }
        }
		else {
            $initialised = 'true'; # Set flag - we have seen at least one   
            $self->initialise_cluster($defragmented_ref, $locus_ref, $j);
		}

		# Update last hit data
		$last_locus{record_id}     = $record_id;
		$last_locus{assigned_name} = $assigned_name;
		$last_locus{assigned_gene} = $assigned_gene;
		$last_locus{probe_gene}    = $probe_gene;
		$last_locus{scaffold}      = $scaffold;
		$last_locus{orientation}   = $orientation;
		$last_locus{$start_token}  = $start;
		$last_locus{$end_token}    = $end;
	}	

}

#***************************************************************************
# Subroutine:  compare_adjacent_loci
# Description: compare two loci and determine whether to merge into one
# Note: this is a critical function in many respects - the logic for 
#       merging or not merging loci is implemented here
#***************************************************************************
sub compare_adjacent_loci {

	my ($self, $locus1_ref, $locus2_ref, $settings_ref) = @_;

	# Get settings
	my $range       = $settings_ref->{range};
	my $start_token = $settings_ref->{start};
	my $end_token   = $settings_ref->{end};
	my $verbose     = $self->{verbose};
	my $mode        = $self->{defragment_mode};


	# Get the current hit values
	my $name             = $locus1_ref->{assigned_name};
	my $gene             = $locus1_ref->{assigned_gene};					
	unless ($gene) { # If there is no assigned gene set, use probe gene	
		$gene = $locus1_ref->{probe_gene};
	}			
	my $scaffold         = $locus1_ref->{scaffold};	
	my $start            = $locus1_ref->{$start_token};
	my $end              = $locus1_ref->{$end_token};
	my $orientation      = $locus1_ref->{orientation};			


	# Get the last hit values
	my $last_name        = $locus2_ref->{assigned_name};
	my $last_gene        = $locus2_ref->{assigned_gene};		
	unless ($last_gene) {  # If there is no assigned gene set, use probe gene	
		$last_gene = $locus2_ref->{probe_gene};
	}	
	my $last_scaffold    = $locus2_ref->{scaffold};	
	my $last_start       = $locus2_ref->{$start_token};
	my $last_end         = $locus2_ref->{$end_token};
	my $last_orientation = $locus2_ref->{orientation};			

	
	# Exclude the obvious cases
	if ($scaffold ne $last_scaffold) { return 0; }  # different scaffolds

	# Check orientation
	if ($orientation ne $last_orientation) {
		unless ($mode eq 'consolidate') { 
			if ($verbose) {
				print "\n\t\t Identified pair of loci that are in range, but different orientations";
			}
			return 0;
		}
		else {
			unless ($last_gene and $gene) { die; } # Should never get here		
		}
	}

	# Take action depending on whether we are DEFRAGMENTING or CONSOLIDATING
	if ($gene and $last_gene) { 
		if ($mode eq 'defragment') {
			if ($gene ne $last_gene) { return 0; }  # different genes
		}
		elsif ($mode eq 'consolidate' or $mode eq 'consolidate2') { 
			# do nothing (these loci can be merged, even though different genes)
		}
		else { # Shouldn't get here
			die;
		}
	}
	else { # Shouldn't get here
		print "\n\t\t ERROR genes not found: ene $gene LAST $last_gene";;
		$devtools->print_hash($locus2_ref);
		die; 
	}
	
	# If on same scaffold in same orientation, determine how far apart 
	my $gap = $start - $last_end;		
	if ($verbose) {
		#print "\n\t\t    - Defragment calculation '$scaffold': '$start'-'$last_end' = $gap";
		print "\n\t\t    - Gap between loci = $gap";
	}

	# Test whether to combine this pair of loci into a single merged locus
	if ($gap <= $range) {  # Combine
		if ($verbose) {
			if ($last_name and $last_gene) {
				print "\n\t\t      - Added pair to cluster: $last_name";
				print "($last_gene [$last_orientation]), $name ($gene [$orientation])"; 		
			}
			else { print "\n\t\t      - Added pair to cluster"; }
		}
		return 1;
	}
	else { # Don't combine
		return 0;
	}
}

#***************************************************************************
# Subroutine:  initialise_cluster
# Description: create the first element in a cluster of associated loci
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
# Description: add a new element to an initialised cluster of associated loci
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
# Description: resolve each of multiple clustered loci to one locus
#***************************************************************************
sub merge_clustered_loci {
	
	my ($self, $defragmented_ref, $to_extract_ref, $reextract) = @_;

	#$devtools->print_hash($defragmented_ref); die;

	# Get screening database table objects
	my $db_ref           = $self->{db};
	my $searches_table   = $db_ref->{searches_table};
	my $active_set_table = $db_ref->{active_set_table};	
	my $extend_count     = '0';
	my @cluster_ids      = keys %$defragmented_ref;
	foreach my $id (@cluster_ids) {

		# Get data for this cluster
		my $cluster_ref = $defragmented_ref->{$id};
		unless ($cluster_ref) { die; }
		my $num_cluster_loci = scalar @$cluster_ref;
		my $extended = $self->merge_cluster($id, $cluster_ref, $to_extract_ref, $reextract);
		if   ($extended) { $extend_count = $extend_count + $extended; }
	}
	return $extend_count;
}

#***************************************************************************
# Subroutine:  merge_cluster
# Description: resolve a cluster of overlapping loci to one single locus
#***************************************************************************
sub merge_cluster {
	
	my ($self, $cluster_id, $cluster_ref, $to_extract_ref, $reextract) = @_;

	# Determine what to extract for this cluster
	my $verbose = $self->{verbose}; # Get 'verbose' flag setting
	my %new_blast_chains;
	my @merged_results;
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
	my $extended;
	my $extract = undef;
	my $num_cluster_loci = scalar @$cluster_ref;		
	if ($verbose) {
		if ($num_cluster_loci > 1) {
			print "\n\t\t    - Merging $num_cluster_loci loci in cluster";
			#$self->show_cluster($cluster_ref, $cluster_id);
		}
	}

	my $num_loci = '0';
	foreach my $locus_ref (@$cluster_ref) {
		
		$num_loci++;			
		my $record_id      = $locus_ref->{record_id};					
		my $digs_result_id = $locus_ref->{digs_result_id};					
		my $start          = $locus_ref->{extract_start};			
		my $end            = $locus_ref->{extract_end};
		$target_name       = $locus_ref->{target_name};					
		$target_datatype   = $locus_ref->{target_datatype};			
		$version           = $locus_ref->{target_version};
		$scaffold          = $locus_ref->{scaffold};			
		$orientation       = $locus_ref->{orientation};
		$organism          = $locus_ref->{organism};
		$probe_type        = $locus_ref->{probe_type};
		unless ($organism) { $organism = $locus_ref->{organism};   }
		unless ($start)    { $start = $locus_ref->{subject_start}; }
		unless ($end)      { $end   = $locus_ref->{subject_end};   }
		if ($verbose) {
			print "\n\t\t    - ID = '$digs_result_id' ($scaffold $orientation $start-$end)";
		}
		#$devtools->print_hash($hit_ref);
							
		# Check if this is a a previously extracted locus
		if ($digs_result_id) {
			$previous_digs_result_ids{$digs_result_id} = $locus_ref;
			$previous_digs_result_id = $digs_result_id;		
		}
		
		# Store this BLAST result as part of a chain if it is new
		else {
			my %data = %$locus_ref;
			$new_blast_chains{$record_id} = \%data; 				
		}
			
		# Record the start and stop parameters so we know whether or not to extend
		if ($lowest_start and $highest_end ) {		
			if ($start < $lowest_start) { $lowest_start = $start; }
			if ($end > $highest_end)    { $highest_end  = $end;   }
		}
		elsif ($start and $end) {
			$highest_end  = $end;
			$lowest_start = $start;
		}
		else { die; } # should never get here			
	}

	# Determine whether or not we need to extract sequences for this cluster
	# Extract if cluster is composed entirely of new loci (no previous result IDs)
	unless ($reextract) {
		unless ($previous_digs_result_id) { 
			if ($verbose) {
				print "\n\t\t # Cluster $cluster_id is comprised entirely of new loci ";
			}
			$extract = 'true';	
		}
	}
	
	# Is this a merge of multiple previously extracted loci?
	my @previous_digs_result_ids = keys %previous_digs_result_ids;
	my $num_previously_extracted_loci_in_cluster = scalar @previous_digs_result_ids;

	if ($num_previously_extracted_loci_in_cluster > 1) {
		my $combined = join (',', @previous_digs_result_ids);
		if ($verbose) {
			print "\n\t\t    - Merging previously extracted loci: ($combined)";
		}
		$extract = 'true';							
	}
	# If it includes a single extracted locus, does this locus need to be extended?	
	elsif ($num_previously_extracted_loci_in_cluster eq 1 and $num_loci > 1) {
			
		# get the indexed query 
		my $ref_id   = shift @previous_digs_result_ids;
		my $data_ref = $previous_digs_result_ids{$ref_id};
		my $start    = $data_ref->{extract_start};			
		my $end      = $data_ref->{extract_end};
		unless ($start)    { $start = $data_ref->{subject_start}; }
		unless ($end)      { $end   = $data_ref->{subject_end};   }
		unless ($lowest_start >= $start and $highest_end <= $end) {	
			$extended++;
			$extract = 'true';
			if ($verbose) {
				print "\n\t\t    - Extending locus: $start, $end: ($lowest_start-$highest_end) ";
			}
		}			
	}

	# If the locus needs to be extracted record the details
	if ($extract) {
		my %extract; # Set the extract params for this cluster
		$extract{target_name}     = $target_name;
		$extract{target_datatype} = $target_datatype;
		$extract{target_version}  = $version;
		$extract{organism}        = $organism;
		$extract{probe_type}      = $probe_type;
		$extract{digs_result_id}  = $previous_digs_result_id;
		$extract{target_name}     = $target_name;
		$extract{scaffold}        = $scaffold;
		$extract{start}           = $lowest_start;	
		$extract{end}             = $highest_end;
		$extract{orientation}     = $orientation;
		$extract{digs_result_ids} = \@previous_digs_result_ids;
		my $num_chains = scalar keys %new_blast_chains;
		if ($num_chains) { $extract{blast_chains} = \%new_blast_chains; }
		push (@$to_extract_ref, \%extract);	
	}
	
	return $extended;
}

#***************************************************************************
# Subroutine:  show_clusters
# Description: print information about clustered loci to the screen
#***************************************************************************
sub show_clusters {

	my ($self, $defragmented_ref) = @_;

	#$devtools->print_hash($defragmented_ref); die;

	my @cluster_ids = keys %$defragmented_ref;
	my $cluster_count;
	foreach my $id (@cluster_ids) {
		$cluster_count++;
		my $cluster_ref = $defragmented_ref->{$id};
		my $cluster_size = scalar @$cluster_ref;
		if ($cluster_size > 1) {
			$self->show_cluster($cluster_ref, $cluster_count);
		}	
	}
}

#***************************************************************************
# Subroutine:  show_cluster
# Description: print information about a cluster of loci to the screen
#***************************************************************************
sub show_cluster {

	my ($self, $cluster_ref, $cluster_id) = @_;

 	#$devtools->print_array($cluster_ref); die;	
	#print "\n";
	
	foreach my $locus_ref (@$cluster_ref) {
   		
   		#$devtools->print_hash($hit_ref); die;
		my $organism      = $locus_ref->{organism};
		my $assigned_name = $locus_ref->{probe_name};
		my $assigned_gene = $locus_ref->{probe_gene};
		my $orientation   = $locus_ref->{orientation};
		my $track_name    = $locus_ref->{track_name};
		
		unless ($assigned_name) {
			$assigned_name = $locus_ref->{assigned_name};
		}
		unless ($assigned_gene) {
			$assigned_gene = $locus_ref->{assigned_gene};
		}

		my $scaffold      = $locus_ref->{scaffold};
		my $start         = $locus_ref->{extract_start};
		my $end           = $locus_ref->{extract_end};
		unless ($start)   { $start = $locus_ref->{subject_start}; }
		unless ($end)     { $end   = $locus_ref->{subject_end};	  }

		my $digs_result_id    = $locus_ref->{digs_result_id};

		print "\n\t\t $cluster_id $organism: ";
		if ($track_name) {
			print "TRACK '$track_name' ";
		}
		print "$assigned_name: $assigned_gene: $scaffold $start-$end ($orientation)";
		if ($digs_result_id) {
			print " (extract ID: $digs_result_id)";
		}
	}			
}

############################################################################
# INTERNAL FUNCTIONS: recording cross-matching
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
# Subroutine:  add_digs_results_to_active_set 
# Description: get all extracted loci for this target file (& probe)
#***************************************************************************
sub add_digs_results_to_active_set {
	
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
	
	return $num_loci;
}

#***************************************************************************
# Subroutine:  get_sorted_digs_results
# Description: get digs_results table rows sorted by scaffold & coordinates
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
	my @fields = qw [ record_id organism 
	                  target_name target_version target_datatype
	                  assigned_name assigned_gene probe_type
	                  scaffold orientation
	                  bitscore gap_openings
	                  query_start query_end 
	                  mismatches align_len
                      evalue_num evalue_exp identity 
                      extract_start extract_end sequence_length ];
	$digs_results_table->select_rows(\@fields, $data_ref, $where);
	
	# Set record_ID as 'digs_result_id' in all results
	foreach my $row_ref (@$data_ref) {
		$row_ref->{digs_result_id} = $row_ref->{record_id};
	}
}

#***************************************************************************
# Subroutine:  get_digs_results_sequences
# Description: get sequences from the digs_results table
#***************************************************************************
sub get_digs_results_sequences {

	my ($self, $data_ref, $where) = @_;;
	
	# Get database tables
	my $db = $self->{db};
	my $digs_results_table  = $db->{digs_results_table};
		
	# Set the fields to get values for
	my @fields = qw [ record_id assigned_name assigned_gene 
	                  probe_type sequence ];
	$digs_results_table->select_rows(\@fields, $data_ref, $where);
}

#***************************************************************************
# Subroutine:  get sorted active set 
# Description: get active set rows, sorted by scaffold, in order of location
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
	my $title       = "DIGS (version: $version_num)";
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

	# Create help menu
	my $program_version = $self->{program_version};
	
    my $HELP   = "\n\t ### DIGS version $program_version";
       $HELP .= "\n\t ### usage: $0 m=[option] -i=[control file] -h=[help]\n";

       $HELP  .= "\n\t ### Main functions\n"; 
	   $HELP  .= "\n\t -m=1  Prepare target files (index files for BLAST)";		
	   $HELP  .= "\n\t -m=2  Do DIGS"; 
	   $HELP  .= "\n\t -m=3  Reassign loci"; 
	   $HELP  .= "\n\t -m=4  Defragment loci"; 
	   $HELP  .= "\n\t -m=5  Consolidate loci\n"; 
	   $HELP  .= "\n\t Genome path = '$ENV{DIGS_GENOMES}'";

	   $HELP  .= "\n\t Run  $0 -e to see information on utility functions\n\n"; 

	print $HELP;
}

############################################################################
# INTERNAL FUNCTIONS: initialisation
############################################################################

#***************************************************************************
# Subroutine:  initialise 
# Description: set up, depending on what option we are running
#***************************************************************************
sub initialise {

	my ($self, $option, $ctl_file) = @_;

	# Die with error message if control file required but not recieved as input
	unless ($ctl_file) {
		unless ($option eq 1) { die "\n\t Option '$option' requires an infile\n\n"; }
		if ($option eq 1) { return; }
	}	

	# Try opening control file
	my @ctl_file;
	my $valid = $fileio->read_file($ctl_file, \@ctl_file);
	unless ($valid) {  # Exit if we can't open the file
		die "\n\t ### Couldn't open control file '$ctl_file'\n\n\n ";
	}

	# If control file looks OK, store the path and parse the file
	$self->{ctl_file} = $ctl_file;
	my $loader_obj = ScreenBuilder->new($self);
	$loader_obj->parse_control_file($ctl_file, $self, $option);

	# Store the ScreenBuilder object (used later)
	$self->{loader_obj} = $loader_obj;

	# Create the output directories
	$loader_obj->create_output_directories($self);

	# Load/create the screening database
	my $db_name = $loader_obj->{db_name};
	unless ($db_name) { die "\n\t Error: no DB name defined \n\n\n"; }
	$self->initialise_screening_db($db_name);

	# SET-UP FOR DIGS SCREENING
	if ($option eq 2) { 
	
		# If we're doing a screen, set up for the screen 
		my $valid = $self->setup_for_digs();
		unless ($valid) { return 0; }
		$self->{defragment_mode} = 'defragment';
	}
	
	# SET-UP FOR REASSIGN
	if ($option eq 3) { 
	
		my $where = '';
		unless ($self->{force}) {
			# Option to enter a WHERE statement
			my $question = "\n\n\t  Enter a WHERE statement to limit reaasign (Optional)";
			$where = $console->ask_question($question);
		}

		# Get the assigned digs_results
		my @reassign_loci;
		$self->get_digs_results_sequences(\@reassign_loci, $where);
		$self->{reassign_loci} = \@reassign_loci;

		# Set up the reference library
		$loader_obj->setup_reference_libraries($self);
	}
	
	# DO SET-UP NEEDED FOR BOTH DEFRAGMENT & CONSOLIDATE
	if ($option eq 4 or $option eq 5) { 

		# Set target sequence files for screening
		my %targets;
		my $num_targets = $loader_obj->set_targets(\%targets);

		# Show error and exit if no targets found
		unless ($num_targets) {
			$loader_obj->show_no_targets_found_error();
		}

		# Set target sequence files for screening
		my %target_groups;
		$loader_obj->set_target_groups(\%targets, \%target_groups);
		$self->{target_groups} = \%target_groups; 
	}

	# DO SET-UP NEEDED FOR DEFRAGMENT ONLY
	if ($option eq 4) { 
		$self->{defragment_mode} = 'defragment';	
		# Set up the reference library
		$loader_obj->setup_reference_libraries($self);
	}
	# DO SET-UP NEEDED FOR CONSOLIDATE ONLY
	elsif ($option eq 5) { 
		$self->set_up_consolidate_tables();
		$self->{defragment_mode} = 'consolidate';

		# Get contig lengths and capture in a table
		$self->calculate_contig_lengths();	

		# Get the parameters for consolidation
		my $range = $self->{consolidate_range};
		my $d_range = $self->{defragment_range};
		unless ($range) { 
			my $question1 = "\n\n\t # Set the range for consolidating digs results";
			$range = $console->ask_int_with_bounds_question($question1, $d_range, $maximum);		
		}
	
		# Set the parameters for consolidation
		my %consolidate_settings;
		$consolidate_settings{range} = $range;
		$consolidate_settings{start} = 'extract_start';
		$consolidate_settings{end}   = 'extract_end';
		$self->{consolidate_settings} = \%consolidate_settings;
	}

	# Create log file
	my $report_dir = $self->{report_dir};
	my $process_id = $self->{process_id};
	my $log_file   = $report_dir . "/log.txt";
	$fileio->append_text_to_file($log_file, "DIGS process $process_id\n");
	$self->{log_file} = $log_file;
}

#***************************************************************************
# Subroutine:  setup_for_digs
# Description: prepare database and DIGS query list prior to screening
#***************************************************************************
sub setup_for_digs {

	my ($self) = @_;

	# Flush active set
	my $db         = $self->{db};
	my $loader_obj = $self->{loader_obj};
	unless ($loader_obj) { die; }  # Sanity checking
	unless ($db)         { die "\n\t Error: no DB defined \n\n\n"; }

	#print "\n\t  Flushing 'active_set' table\n";
	my $active_set_table = $db->{active_set_table};
	$active_set_table->flush();
	
	# Index previously executed searches
	my %done;
	$self->index_previously_executed_searches(\%done);
	
	# Get the list of queries that have been completed 
	my %queries;
	$loader_obj->{previously_executed_searches} = \%done;

	# Set up the DIGS screen
	my $total_queries = $loader_obj->setup_screen($self, \%queries);
	unless ($total_queries)  { 
		print "\n\t  Exiting DIGS setup";	
		return 0;;
	}
		
	# Record queries 
	$self->{queries}       = \%queries;
	$self->{total_queries} = $total_queries;

	return 1;
}

#***************************************************************************
# Subroutine:  initialise_screening_db
# Description: load a DIGS screening database (create if doesn't exist) 
#***************************************************************************
sub initialise_screening_db {

	my ($self, $db_name) = @_;

	# Create the screening DB object
	my $db_obj = ScreeningDB->new($self);

	# Check if this screening DB exists, if not then create it
	my $db_exists = $db_obj->does_db_exist($db_name);
	unless ($db_exists) {
		$db_obj->create_screening_db($db_name);	
	}
	
	# Load the table handles into screening database object
	print   "\n\n\t  Connecting to DB:  $db_name";
	$db_obj->load_screening_db($db_name);	
	$self->{db} = $db_obj; # Store the database object reference 
}

#***************************************************************************
# Subroutine:  set_up_consolidate_tables
# Description:
#***************************************************************************
sub set_up_consolidate_tables {

	my ($self) = @_;

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
	my $contigs_exists = $db_ref->does_table_exist('contigs');
	unless ($contigs_exists) {
		$db_ref->create_contigs_table($dbh);
	}

	# Load tables
	$db_ref->load_loci_table($dbh);
	$db_ref->load_loci_chains_table($dbh);
	$db_ref->load_contigs_table($dbh);

	# Get table references and set up for this consolidation process
	my $loci_table        = $db_ref->{loci_table};
	my $loci_chains_table = $db_ref->{loci_chains_table};
	my $contigs_table     = $db_ref->{contigs_table};
	$loci_table->flush();
	$loci_chains_table->flush();
	$contigs_table->flush();
}

#***************************************************************************
# Subroutine:  calculate_contig_lengths
# Description: 
#***************************************************************************
sub calculate_contig_lengths {

	my ($self, $contigs_ref) = @_;

	# Get database and tables
	my $db_ref = $self->{db};
	my $digs_results     = $db_ref->{digs_results_table};
	my $contigs_table    = $db_ref->{contigs_table};
	my $genome_use_path  = $self->{genome_use_path};
	my $target_group_ref = $self->{target_groups};
	
	my $question1 = "\n\n\t # Refresh contig length table";
	my $refresh = $console->ask_yes_no_question($question1);		
	#my $refresh = 'y';

	if ($refresh eq 'y') {
	
		my @targets;
		my @fields = qw [ organism target_datatype target_version target_name ];
		$digs_results->select_distinct(\@fields, \@targets);
		foreach my $target_ref(@targets) {
			
			# Read the file and get the lengths of each contig	
			# Get the target details (and thus the target path)	
			#my $target_path = $self->get_target_file_path($target_ref);
			my $organism        = $target_ref->{organism};
			my $target_name     = $target_ref->{target_name};
			my $target_datatype = $target_ref->{target_datatype};
			my $target_version  = $target_ref->{target_version};
			my @genome = ( $organism , $target_datatype, $target_version, $target_name );
			my $target_id       = join ('|', @genome);
			my $target_group = $target_group_ref->{$target_id};
			unless ($target_group) { die; 
				print " \n\t No target group found for TARGET ID $target_id\n\n"; 
        		sleep 1;
				next;
			}
		
			# Construct the path to this target file
			my @path;
			push (@path, $genome_use_path);
			push (@path, $target_group);
			push (@path, $organism);
			push (@path, $target_datatype);
			push (@path, $target_version);
			push (@path, $target_name);
			my $target_path = join ('/', @path);
			my @contigs;
			$fileio->read_fasta($target_path, \@contigs, 'true');
			
			foreach my $contig_ref (@contigs) {
			
				my $header = $contig_ref->{header};
				my $length = $contig_ref->{seq_length};
				unless ($length and $header) {
					$devtools->print_hash($contig_ref); die;
				}
				
				$contig_ref->{organism}        = $organism;
				$contig_ref->{target_datatype} = $target_datatype;
				$contig_ref->{target_version}  = $target_version;
				$contig_ref->{target_name}     = $target_name;
				$contig_ref->{scaffold}        = $header;
				print "\n\t # $header: $length";
				
				$contigs_table->insert_row($contig_ref);
			}			
		}
	}
}

############################################################################
# Development
############################################################################

#***************************************************************************
# Subroutine:  prepare_locus_update 
# Description: 
#***************************************************************************
sub prepare_locus_update {

	my ($self, $loci_ref) = @_;

	# Get parameters from self
	foreach my $hit_ref (@$loci_ref) {
	
		$hit_ref->{extract_start}   = $hit_ref->{start};
		$hit_ref->{extract_end}     = $hit_ref->{end};
		$hit_ref->{sequence}        = 'NULL';
		$hit_ref->{sequence_length} = 0;
		#$devtools->print_hash($hit_ref); die;
	}
}

############################################################################
# EOF
############################################################################
