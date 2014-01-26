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
		
		# Paths and member variables
		blast_bin_path         => $parameter_ref->{blast_bin_path},
		genome_use_path        => $parameter_ref->{genome_use_path},
		output_path            => $parameter_ref->{output_path},

		# Database variables
		db_name                => '',
		server                 => '',  
		username               => '',
		password               => '',
	
		# Member classes 
		blast_obj              => $parameter_ref->{blast_obj},
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# SECTION: Handler functions 
############################################################################

#***************************************************************************
# Subroutine:  run_screen_function
# Description: handler for various utility processes 
#***************************************************************************
sub run_screen_function {

	my ($self, $option, $ctl_file) = @_;

	# Initialise
	print "\n\t ### Initialising database-guided genome screening\n";
	$self->initialise($ctl_file);
	my $db = $self->{db};
	
	# Hand off to functions, based on the options received
	if ($option eq 1) { 
		$self->run_screen($ctl_file);	
	}
	elsif ($option eq 2) {  
		die "\n\t ### Unimplemented!\n\n\n";
		#$self->reassign($loader_obj);	
	}
	elsif ($option eq 3) {  
		$db->summarise_db();
	}
	elsif ($option eq 4) {  
		$self->retrieve();
	}
	elsif ($option eq 5) {  
		$db->flush_screening_db();
	}
	elsif ($option eq 6) {  
		$db->drop_screening_db();    
	}
}

############################################################################
# SECTION: Initialise
############################################################################

#***************************************************************************
# Subroutine:   initialise 
# Description:  do initial checks
#***************************************************************************
sub initialise {

	my ($self, $ctl_file) = @_;
	
	# Try opening control file first
	my @ctl_file;
	my $valid = $fileio->read_input_file($ctl_file, \@ctl_file);
	unless ($valid) {
		print "\n\t ### Couldn't open control file '$ctl_file'\n\n\n ";
		exit;
	}	

	# Check the size of the process directory
	$self->check_process_dir_status();
	
	# Get parameters
	#print "\n\t ### Reading control file";
	my $loader_obj = ScreenBuild->new($self);
	$loader_obj->parse_control_file($ctl_file);
		
	# Load screening database (includes some MacroLineage Tables
	$loader_obj->set_screening_db();
	$self->{db} = $loader_obj->{db};
	
}

############################################################################
# SECTION: Running a round of bidirectional BLAST screening
############################################################################

#***************************************************************************
# Subroutine:   run_screen
# Description:  handler function to set up and execute a round of screening
#***************************************************************************
sub run_screen {

	my ($self, $ctl_file) = @_;

	# Initialise target sequence library
	my $genome_obj = GenomeControl->new($self); 
	$genome_obj->refresh_genomes();
	
	# Set up the screening queries
	print "\n\n\t ### Setting up a DIGS screen from control file '$ctl_file'";
	my %queries;
	my $loader_obj = ScreenBuild->new($self);
	$loader_obj->set_up_screen($self, \%queries, $ctl_file);

	# Run reciprocal BLAST pipeline
	print "\n\n\t ### Starting database-guided genome screening\n";
	$self->do_screening_process(\%queries);
}

#***************************************************************************
# Subroutine:   do pipeline screen 
# Description:  run a round of bidrectional BLAST 
#***************************************************************************
sub do_screening_process {

	my ($self, $queries_ref) = @_;

	# Get relevant member variables and objects
	my $db_ref       = $self->{db};

	# Iterate through and excute screens
	my @probes = keys %$queries_ref;
	foreach my $probe_name (@probes) {
		
		# Get the array of queries for this target file
		my $probe_queries = $queries_ref->{$probe_name};
		foreach my $query_ref (@$probe_queries) {  
			
			# Show status 
			my $probe_id    = $query_ref->{probe_id};
			my $target_name = $query_ref->{target_name}; # $target refers to the target file
			print "\n\t # probing target $target_name with probe $probe_id ";   
			
			# Do the BLAST search
			$self->search($query_ref);	
	
			# Extract matches 
			my @extracted;
			$self->extract($query_ref, \@extracted);	
			
			# Classify extracted sequences by BLAST comparison to reference library
			$self->assign(\@extracted);	
		}	
	}
	
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
	my $db_ref     = $self->{db};
	my $blast_obj  = $self->{blast_obj};
	my $tmp_path   = $self->{tmp_path};
	my $min_length = $self->{seq_length_minimum};
	unless ($min_length) { $min_length = $default_min_seqlen; }
	unless ($db_ref and $blast_obj and $tmp_path) { die; } # Sanity checking

	# Get query details
	my $probe_id     = $query_ref->{probe_id};
	my $probe_name   = $query_ref->{probe_name};
	my $probe_gene   = $query_ref->{probe_gene};
	my $probe_path   = $query_ref->{probe_path};
	my $organism     = $query_ref->{organism};
	my $target_name  = $query_ref->{target_name};
	my $target_path  = $query_ref->{target_path};
	my $blast_alg    = $query_ref->{blast_alg};
	my $cutoff       = $query_ref->{bitscore_cutoff};
	my $result_file  = $tmp_path . "/$probe_id" . "_$target_name.blast_result.tmp";
	
	# Create path to BLAST 
	my $blast_path  = $self->{blast_bin_path};
	my $blast_bin_path;
	if ($blast_path) { $blast_bin_path = $blast_path . $blast_alg; }	
	else             { $blast_bin_path = $blast_alg;               }	

	# Do the BLAST search
	my $command  = "$blast_bin_path -query $probe_path -subject $target_path ";
	   $command .= " -out $result_file -outfmt 7"; 
	system $command;		
	#print $command; die; # DEBUG

	# Parse out the alignment
	my @hits;
	$blast_obj->parse_tab_format_results($result_file, \@hits, $cutoff);
	
	# Remove the result
	my $rm_command = "rm $result_file";
	system $rm_command;

	# Rename coordinate fields to match DB
	foreach my $hit_ref (@hits) {
		$hit_ref->{subject_start} = $hit_ref->{aln_start};
		$hit_ref->{subject_end}   = $hit_ref->{aln_stop};
		$hit_ref->{query_end}     = $hit_ref->{query_stop};
	}

	# Index loci previously extracted from this target file
	my %extracted;
	$self->index_extracted_loci($target_name, \%extracted);

	# Store new new matches that meet conditions
	my $blast_results_table = $db_ref->{blast_results_table};
	my $i = 0;
	foreach my $hit_ref (@hits) {
		$i++;
		
		# Check if redundant with Loci table
		my $skip = $self->check_if_locus_extracted($hit_ref, \%extracted);
		if ($skip) {  next; }

		# Skip sequences that are too short
		my $start  = $hit_ref->{subject_start};
		my $end    = $hit_ref->{subject_end};
		if ($end - $start < $min_length) {  next; }
		
		# Record hit in BLAST results table
		$hit_ref->{organism}    = $organism;
		$hit_ref->{chunk_name}  = $target_name;
		$hit_ref->{probe_id}    = $probe_id;
		$hit_ref->{probe_name}  = $probe_name;
		$hit_ref->{probe_gene}  = $probe_gene;
		$hit_ref->{probe_type}  = $query_ref->{probe_type};
		print "\n\t # Match number $i to $probe_name, $probe_gene in '$target_name'";
		$blast_results_table->insert_row($hit_ref);
	} 

	# Update the status table
	my $status_table = $db_ref->{status_table};
	$status_table->insert_row($query_ref);
		
}

############################################################################
# SECTION: Extract 
############################################################################

#***************************************************************************
# Subroutine:  extract
# Description: extract sequences that matched probes in a BLAST search 
#***************************************************************************
sub extract {
	
	my ($self, $query_ref, $extracted_ref) = @_;

	# Get paths
	my $blast_bin_path     = $self->{blast_bin_path};
	
	# Get data structures and variables from self
	my $db_ref = $self->{db};
	my $blast_results_table = $db_ref->{blast_results_table};
	my $extracted_table     = $db_ref->{extracted_table}; 

	# Index extracted BLAST results 
	my %extracted;
	my @fields = qw [ blast_id ];
	my @blast_ids;
	$extracted_table->select_rows(\@fields, \@blast_ids);	
	foreach my $hit_ref (@blast_ids) {
		my $blast_id = $hit_ref->{blast_id};
		$extracted{$blast_id} = $hit_ref;	
	}

	# Get query params we need
	my $probe_name  = $query_ref->{probe_name};
	my $probe_gene  = $query_ref->{probe_gene};
	my $target_name = $query_ref->{target_name}; # target file name
	my $target_path = $query_ref->{target_path};

	# Get all BLAST results from table (ordered by sequential targets)
	my $where = " WHERE Chunk_name = '$target_name'
	                AND probe_name = '$probe_name' 
                    AND probe_gene = '$probe_gene'
	              ORDER BY scaffold, subject_start";
	my @matches;
	@fields = qw [ record_id probe_name probe_gene probe_type
				   e_value_num e_value_exp bit_score align_len orientation 
	               organism scaffold chunk_name subject_start subject_end 
		           query_start query_end ];
	$blast_results_table->select_rows(\@fields, \@matches, $where); 
	#$devtools->print_array(\@matches); exit;

	# Store all outstanding matches as sequential, target-ordered sets 
	foreach my $hit_ref (@matches) {
		
		my $record_id   = $hit_ref->{record_id};
		if ($extracted{$record_id}) { next; }

		my $start       = $hit_ref->{subject_start};
		my $end         = $hit_ref->{subject_end};
		my $orientation = $hit_ref->{orientation};
		my $scaffold    = $hit_ref->{scaffold};
		
		# Parsing for blastdbcmd
		my @gi = split(/\|/,$scaffold);	
		if (scalar(@gi) > 1) {
			$scaffold = $gi[1];
		}

		# Create the command
		# Command example: 
		# /bin/blast/blastdbcmd -db hs_alt_HuRef_chrX.fa -entry 157734237 
		# -range 10-60 -strand minus
		my $command = $blast_bin_path . "blastdbcmd -db $target_path";
		$command .= " -entry $scaffold ";
		$command .= " -range $start-$end ";
		if ($orientation eq '-ve') { $command .= ' -strand minus '; }
		
		# Execute the command
		my @sequence = `$command`;
		shift @sequence;  # Remove header
		my $sequence = join ('', @sequence);
		$sequence =~ s/\n//g;
		
		# Check we got sequences
		unless ($sequence) {
			print "FAILED TO EXTRACT sequence using command '$command'\n";
		}
		else {	
			# Set sequence length
			my $seq_length = length $sequence;
			$hit_ref->{sequence_length} = $seq_length;
			$hit_ref->{sequence} = $sequence;
			push (@$extracted_ref, $hit_ref);
		}
	}
}

############################################################################
# SECTION: Assign 
############################################################################

#***************************************************************************
# Subroutine:  assign
# Description: assign sequences that matched probes in a BLAST search 
#***************************************************************************
sub assign {
	
	my ($self, $extracted_ref) = @_;
	
	# Get data structures and variables from self
	my $db_ref        = $self->{db};
	my $result_path   = $self->{tmp_path};
	my $blast_obj     = $self->{blast_obj};
	my $table = $db_ref->{extracted_table}; 
	unless ($result_path) { die; }

	# Iterate through the matches
	foreach my $hit_ref (@$extracted_ref) {
		
		# Copy the coordinates AS EXTRACT coordinates	
		my $extract_start = $hit_ref->{subject_start};
		my $extract_end   = $hit_ref->{subject_end};
		my $blast_id      = $hit_ref->{record_id};
		my $sequence      = $hit_ref->{sequence};
		my $genome        = $hit_ref->{organism};
		my $probe_type    = $hit_ref->{probe_type};
		
		# Set the linking to the BLAST result table
		$hit_ref->{blast_id} = $blast_id;
		unless ($sequence) { # Sanity checking
			print "\n\t ###### WARNING NO Sequence!"; 
			next;
		}	
		
		# Make a file for BLAST
		my $fasta      = ">$blast_id\n$sequence";
		my $query_file = $result_path . $blast_id . '.fas';
		$fileio->write_text_to_file($query_file, $fasta);
		my $result_file = $result_path . $blast_id . '.blast_result';
			
		# Do the BLAST accoridnmg to the type of sequence (AA or NA)
		print "\n\t ##### BLAST HIT $blast_id, $genome";
		my $blast_alg;
		my $lib_path;
		unless ($probe_type) { die; }
		if ($probe_type eq 'UTR') {
			$lib_path  = $self->{blast_utr_lib_path};
			$blast_alg = 'blastn';
		}
		elsif ($probe_type eq 'ORF') {
			$lib_path  = $self->{blast_orf_lib_path};
			$blast_alg = 'blastx';
		}
		else { die; }

		unless ($lib_path) { die "\n\n\t NO BLAST LIBRARY defined\n\n\n"; }
		$blast_obj->blast($blast_alg, $lib_path, $query_file, $result_file);
		my @results;
		$blast_obj->parse_tab_format_results($result_file, \@results);

		# Get the best match from this file
		my $top_match = shift @results;
		my $query_start   = $top_match->{query_start};
		my $query_end     = $top_match->{query_stop};
		my $subject_start = $top_match->{aln_start};
		my $subject_end   = $top_match->{aln_stop};
		my $assigned_to   = $top_match->{scaffold};	
		
		unless ($assigned_to) {	
			print "\n\t ### No match found in reference library\n";
			next; 
		}
		
		# Split assigned to into (i) refseq match (ii) refseq description (e.g. gene)	
		my @assigned_to = split('_', $assigned_to);
		my $assigned_to_gene = pop @assigned_to;
		$assigned_to = join ('_', @assigned_to);
		print " assigned to: $assigned_to: $assigned_to_gene!";
		#$devtools->print_hash($hit_ref);

		# Adjust to get extract coordinates using the top match
		my $real_st  = ($extract_start + $top_match->{query_start} - 1);	
		my $real_end = ($extract_start + $top_match->{query_stop} - 1);
		
		$hit_ref->{assigned_to}      = $assigned_to;
		$hit_ref->{assigned_to_gene} = $assigned_to_gene;
		$hit_ref->{realex_start}     = $real_st;
		$hit_ref->{realex_end}       = $real_end;
		$hit_ref->{extract_start}    = $extract_start;
		$hit_ref->{extract_end}      = $extract_end;
		$hit_ref->{identity}         = $top_match->{identity};
		$hit_ref->{mismatches}       = $top_match->{mismatches};
		$hit_ref->{gap_openings}     = $top_match->{gap_openings};
		$hit_ref->{query_end}        = $query_end;
		$hit_ref->{query_start}      = $query_start;
		$hit_ref->{subject_end}      = $subject_end;
		$hit_ref->{subject_start}    = $subject_start;
		
		# Insert the data
		my $extract_id = $table->insert_row($hit_ref);

		# Clean up
		my $command1 = "rm $query_file";
		my $command2 = "rm $result_file";
		system $command1;
		system $command2;
	}
}

############################################################################
# Utility functions
############################################################################

#***************************************************************************
# Subroutine:  index previously extracted loci
# Description: Index loci previously extracted from this target file
#***************************************************************************
sub index_extracted_loci {
	
	my ($self, $target_name, $previously_extracted_ref) = @_;

	# Get relevant variables and objects
	my $db_ref          = $self->{db};
	my $extracted_table = $db_ref->{extracted_table}; 
	my @locus_data;
	my @fields = qw [ record_id scaffold assigned_to
	                  extract_start extract_end ];
	my $where = " WHERE chunk_name = '$target_name' "; 
	$extracted_table->select_rows(\@fields, \@locus_data, $where);	
	foreach my $row_ref (@locus_data) {
		my $scaffold = $row_ref->{scaffold};
		if ($previously_extracted_ref->{$scaffold}) {
			my $array_ref = $previously_extracted_ref->{$scaffold};
			push (@$array_ref, $row_ref);
		}
		else {
			my @array;
			push (@array, $row_ref);
			$previously_extracted_ref->{$scaffold} = \@array;
		}
	}
}

#***************************************************************************
# Subroutine:  check if locus previously extracted
# Description: 
#***************************************************************************
sub check_if_locus_extracted {
	
	my ($self, $hit_ref, $extracted_ref) = @_;

	my $skip = undef;
	my $scaffold = $hit_ref->{scaffold};
	unless ($scaffold) { die; } # Sanity checking
	if ($extracted_ref->{$scaffold}) {
		# Get hit coordinates
		my $hit_start = $hit_ref->{subject_start};
		my $hit_end   = $hit_ref->{subject_end};
		my $loci_ref = $extracted_ref->{$scaffold};
		foreach my $locus_ref (@$loci_ref) {
			my $locus_start = $locus_ref->{extract_start};
			my $locus_end   = $locus_ref->{extract_end};
			if ($hit_start >= $locus_start and $hit_end <= $locus_end) {
				#$devtools->print_hash($locus_ref);
				my $locus_id      = $locus_ref->{record_id};
				my $locus_erv     = $locus_ref->{assigned_to};
				my $scaffold      = $hit_ref->{scaffold};
				my $subject_start = $hit_ref->{subject_start};
				my $subject_end   = $hit_ref->{subject_end};
				$skip = 'true';
				print "\n\t # Will skip hit in $scaffold ($subject_start-$subject_end)";
				print "\n\t # Because it is redundant with locus $locus_id ($locus_erv)";
			}
		}
	}
	return $skip;
}

#***************************************************************************
# Subroutine:  reassign
# Description: reassign sequences in the extracted_table (for use after
#              the reference library has been updated)
#***************************************************************************
sub reassign {
	
	my ($self, $loader_obj) = @_;

	# Set up to perform the reassign process
	print "\n\t # Reassigning Extracted table";
	my @assigned_seqs;
	$self->initialise_reassign($loader_obj, \@assigned_seqs);

	# Get data structures and variables from self
	my $blast_obj       = $self->{blast_obj};
	my $result_path     = $self->{report_dir};
	my $db              = $self->{db};
	my $extracted_table = $db->{extracted_table};
	unless ($result_path) { die; }
	
	# Iterate through the matches
	foreach my $row_ref (@assigned_seqs) {
		
		my $sequence        = $row_ref->{sequence};
		my $probe_type      = $row_ref->{probe_type};
		my $blast_id        = $row_ref->{record_id};
		my $previous_assign = $row_ref->{assigned_to};
		my $previous_gene   = $row_ref->{assigned_to_gene};
		my $extract_start   = $row_ref->{subject_start};
		my $extract_end     = $row_ref->{subject_end};
		my $genome          = $row_ref->{organism};
		unless ($sequence) { die "\n\t NO Sequence!"; }	# Sanity checking
		print "\n\t Redoing assign for record ID $blast_id assigned to $previous_assign";
		
		# Make a file for BLAST
		my $fasta      = ">$blast_id\n$sequence";
		my $query_file = $result_path . $blast_id . '.fas';
		$fileio->write_text_to_file($query_file, $fasta);
		my $result_file = $result_path . $blast_id . '.blast_result';
			
		# Do the BLAST according to the type of sequence (AA or NA)
		print "\n\t ##### BLAST BACK for $blast_id, $genome";
		my $blast_alg;
		my $lib_path;
		unless ($probe_type) { die; }
		if ($probe_type eq 'UTR') {
			$lib_path  = $self->{blast_utr_lib_path};
			$blast_alg = 'blastn';
		}
		elsif ($probe_type eq 'ORF') {
			$lib_path  = $self->{blast_orf_lib_path};
			$blast_alg = 'blastx';
		}
		else { die; }
		unless ($lib_path) { die "\n\n\t NO BLAST LIBRARY defined\n\n\n"; }
		$blast_obj->blast($blast_alg, $lib_path, $query_file, $result_file);
		my @results;
		$blast_obj->parse_tab_format_results($result_file, \@results);

		# Get the best match from this file
		my $top_match = shift @results;
		my $query_start   = $top_match->{query_start};
		my $query_end     = $top_match->{query_stop};
		my $subject_start = $top_match->{aln_start};
		my $subject_end   = $top_match->{aln_stop};
		my $assigned_to   = $top_match->{scaffold};	
	
		# Split assigned to into (i) refseq match (ii) refseq description (e.g. gene)	
		my @assigned_to = split('_', $assigned_to);
		my $assigned_to_gene = pop @assigned_to;
		$assigned_to = join ('_', @assigned_to);
		
		if ($assigned_to ne $previous_assign 
		or  $assigned_to_gene ne $previous_gene) {
			
			print "\n\t ##### Reassigned $blast_id from $previous_assign ($previous_gene)";
			print " to $assigned_to ($assigned_to_gene)";
			# Adjust to get extract coordinates using the top match
			my $real_st  = ($extract_start + $top_match->{query_start} - 1);	
			my $real_end = ($extract_start + $top_match->{query_stop} - 1);
			my %update_row;
			$update_row{assigned_to}   = $assigned_to;
			$update_row{assigned_to_gene} = $assigned_to_gene;
			$update_row{realex_start}  = $real_st;
			$update_row{realex_end}    = $real_end;
			$update_row{extract_start} = $extract_start;
			$update_row{extract_end}   = $extract_end;
			$update_row{identity}      = $top_match->{identity};
			$update_row{mismatches}    = $top_match->{mismatches};
			$update_row{gap_openings}  = $top_match->{gap_openings};
			$update_row{query_end}     = $query_end;
			$update_row{query_start}   = $query_start;
			$update_row{subject_end}   = $subject_end;
			$update_row{subject_start} = $subject_start;
				
			# Insert the data
			my $where = " WHERE Record_id = $blast_id ";
			$extracted_table->update(\%update_row, $where);
		}
		
		# CLEAN UP
		my $system1 = "rm $query_file";
		my $system2 = "rm $result_file";
		system $system1;
		system $system2;
	}
}

#***************************************************************************
# Subroutine:  initialise_reassign 
# Description: 
#***************************************************************************
sub initialise_reassign {

	my ($self, $loader_obj, $assigned_seqs_ref) = @_;

	$self->show_title();

	# Create a unique ID and report directory for this run
	my $output_path = $self->{output_path};
	my $process_id  = $self->{process_id};
	my $db          = $self->{db};
	my $db_name     = $db->{db_name};
	
	unless ($process_id and $output_path) { die; }
	my $report_dir  = $output_path . $process_id . '/';
	$fileio->create_unique_directory($report_dir);
	my $tmp_path = $report_dir . '/tmp';
	$fileio->create_unique_directory($tmp_path);
	$self->{tmp_path}   = $tmp_path . '/';
	$self->{report_dir} = $report_dir;
	$loader_obj->{report_dir} = $report_dir;
	
	# Get the assigned data
	my $extracted_table = $db->{extracted_table};
	my @fields  = qw [ record_id probe_type assigned_to assigned_to_gene 
	                       subject_start subject_end sequence organism ];
	$extracted_table->select_rows(\@fields, $assigned_seqs_ref);

	# Transfer parameters from loader to this obj

	# Set up the reference library
	if ($loader_obj->{reference_aa_fasta}) {
		$loader_obj->load_aa_fasta_reference_library();
	}
	if ($loader_obj->{reference_nt_fasta}) {
		$loader_obj->load_nt_fasta_reference_library();
	}

	$self->{seq_length_minimum}    = $loader_obj->{seq_length_minimum};
	$self->{bit_score_min_tblastn} = $loader_obj->{bit_score_min_tblastn};
	$self->{bit_score_min_blastn}  = $loader_obj->{bit_score_min_blastn};
	$self->{blast_orf_lib_path}    = $loader_obj->{blast_orf_lib_path};
	$self->{blast_utr_lib_path}    = $loader_obj->{blast_utr_lib_path};
}

############################################################################
# SECTION: Retrieve data from screening database
############################################################################

#***************************************************************************
# Subroutine:   retrieve 
# Description:  retrieve data from screening database
#***************************************************************************
sub retrieve {

	my ($self, $ctl_file) = @_;

	# Get params from self	
	my $select = $self->{select_list};
	my $where  = $self->{where_statement};
	unless ($select) { die; }

	my $db = $self->{db};
	my @select = split(/'/, $select);

	# Get params from self	
	my @sequences;
	$devtools->print_hash($self); die;
	$db->retrieve_sequences(\@sequences, \@select, $where);
	die;

}

############################################################################
# SECTION: Utility Fxns
############################################################################

#***************************************************************************
# Subroutine:  check_process_dir_status
# Description: check the size of the process directory and warn if it is
#              getting very large 
#***************************************************************************
sub check_process_dir_status {

	my ($self) = @_;

	# Read in the process directory
	my $output_path = $self->{output_path};
	unless ($output_path) { die; }

	my @process_dir;
	$fileio->read_directory_to_array($output_path, @process_dir);
	my $num_folders = scalar @process_dir;

	if ($num_folders > $process_dir_warn) {
		print "\n\t ### Process directory '$output_path' contains > $process_dir_warn folders";
		print "\n\t #     consider cleaning up the contents.";
		sleep 2;
	}
}

#***************************************************************************
# Subroutine:  show_title
# Description: show command line title blurb 
#***************************************************************************
sub show_title {

	$console->refresh();
	my $title       = 'Database-Integrated Genome Screening (DIGS)';
	my $version     = '1.0';
	my $description = 'Sequence database screening using BLAST and MySQL';
	my $author      = 'Robert J. Gifford';
	my $contact	    = '<robert.gifford@glasgow.ac.uk>';
	$console->show_about_box($title, $version, $description, $author, $contact);
}

############################################################################
# EOF
############################################################################
