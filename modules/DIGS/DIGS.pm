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
use DIGS::Initialise;    # Initialises the DIGS tool
use DIGS::ScreenBuilder; # To set up a DIGS run
use DIGS::Defragment;    # Defragment tools
use DIGS::Consolidate;   # Consolidate locus functions
use DIGS::Extract;       # Extracting sequences for FASTA files using BLAST
use DIGS::Classify;      # Classifying sequences using BLAST
use DIGS::CrossMatch;    # Recording cross-matching during DIGS

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
# Description: create new DIGS 'object'
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Create objects for screening
	my $crossmatch_obj = CrossMatch->new($parameter_ref);

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
		crossmatch_obj         => $crossmatch_obj,
		
		# MySQL database connection parameters
		mysql_username         => $parameter_ref->{mysql_username}, 
		mysql_password         => $parameter_ref->{mysql_password},
		db_name                => '',   # Obtained from control file or console
		mysql_server           => '',   # Obtained from control file or console

		# Parameters for DIGS
		query_na_fasta         => '',   # Obtained from control file
		query_aa_fasta         => '',   # Obtained from control file
		reference_na_fasta     => '',   # Obtained from control file
		reference_aa_fasta     => '',   # Obtained from control file
		

                aa_reference_library   => '',   # Obtained from control file
		na_reference_library   => '',   # Obtained from control file

		bitscore_min_tblastn   => '',   # Obtained from control file
		bitscore_min_blastn    => '',   # Obtained from control file
		seq_length_minimum     => '',   # Obtained from control file
		extract_buffer         => '',   # Obtained from control file
		defragment_mode        => '',   # Determined by input options
		
		# Paths used in DIGS process
		genome_use_path        => $parameter_ref->{genome_use_path},
		output_path            => $parameter_ref->{output_path},
		tmp_path               => '',   # Created during set up
		blast_threads          => '',   # Obtained from control file
		
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

	$self->show_title();  

    # Initialise the screen based on (i) options given and (ii) control file 
	my $valid = $self->initialise($ctl_file, $option);  
	
	if ($valid) {
	
		# Hand off to the appropriate function, depending on the option received
		$self->hand_off_to_digs_fxns($option);

		# Show final summary and exit message
		$self->wrap_up($option);
	}
}

############################################################################
# SETTING UP - Processing and validating input, handing off to subroutines
############################################################################

#***************************************************************************
# Subroutine:  initialise
# Description: check input file and options received are valid, & initialise
#***************************************************************************
sub initialise {

	my ($self, $ctl_file, $option) = @_;

	my $valid          = undef;
	my $initialise_obj = Initialise->new($self);	
	
	if ($ctl_file) { # Try to initialise using the input file
		$valid = $initialise_obj->initialise($self, $option, $ctl_file);
	}
	elsif ($option eq 1) { 
		$valid = 1; # An infile is optional for option -m=1
	}
	elsif ($option > 1 and $option <= 5) { # Show error if no infile
		print "\n\t  Option '-m=$option' requires an infile\n\n";	
	}


	if ($option > 5) { # # Unavailable option error
		print "\n\t  Unrecognized option '-m=$option'\n\n";
	}


	return $valid;
}

#***************************************************************************
# Subroutine:  hand_off_to_digs_fxns
# Description: hand off to the appropriate function, depending on options 
#***************************************************************************
sub hand_off_to_digs_fxns {

	my ($self, $option) = @_;

	# Main screening functions
	if ($option eq 1) {      # Check the target sequences are formatted for BLAST		
		$self->index_target_files_for_blast();	
	}
	elsif ($option eq 2) { 	 # Screen
		$self->do_digs();	
	}
	elsif ($option eq 3) {   # Reassign data in digs_results table
		$self->reassign();	
	}
	
	# Defragmenting screening results
	elsif ($option eq 4 or $option eq 5) {	

		if ($option eq 4) {
			# Create a defragmenter module
			my $defragment_obj = Defragment->new($self);
		    # Interactively defragment contiguous hits to the same gene 	
			$defragment_obj->interactive_defragment();
		}
		elsif ($option eq 5) {  
			# Create a consolidate module
			my $consolidate_obj = Consolidate->new($self);
			# Combine hits to different genes into higher order locus structures
            $consolidate_obj->consolidate_loci();
		}
	}
	else { 
		print  "\n\n\t # Unrecognised option m='$option'\n\n\n\n";
		exit;
	} # Shouldn't get here

}

############################################################################
# MAIN DIGS FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  index_target_files_for_blast
# Description: create index files for all target databases
#***************************************************************************
sub index_target_files_for_blast {

	my ($self) = @_;

	# Format targets files for BLAST searching
	my $target_db_obj = TargetDB->new($self);
	$target_db_obj->format_targets_for_blast();

}

#***************************************************************************
# Subroutine:  do_digs
# Description: do the core database-integrated genome screening processes
#***************************************************************************
sub do_digs {

	my ($self, $mode) = @_;

	# Iterate through the list of DIGS queries, dealing each in turn 
	# Each DIGS query constitutes a probe sequence and a target FASTA file
	my $current_query_num = 0;
	my $queries_ref = $self->{queries};
	unless ($queries_ref) { die; }   # Sanity checking
	my @probes = keys %$queries_ref; # Get the list of queries
	print "\n\t ### Starting database-integrated genome screening";
	foreach my $probe_name (@probes) {
		
		# Get the array of queries for this target file
		my $probe_queries = $queries_ref->{$probe_name};
		foreach my $query_ref (@$probe_queries) {  
	
			# Increment query count
			$current_query_num++;		
			$self->{current_query_num} = $current_query_num;

			# DIGS round one
			my @to_extract;
			my @to_delete;
			$self->run_digs_search_phase($query_ref, \@to_extract, \@to_delete);
			
			# DIGS round two
			$self->run_digs_classify_phase($query_ref, \@to_extract, \@to_delete);
		
			# Show a status update in the console
			$self->show_digs_progress();			
		}	
	}
}


#***************************************************************************
# Subroutine:  run_digs_search_phase
# Description: execute a BLAST query and get non-redundant results
#***************************************************************************
sub run_digs_search_phase {

	my ($self, $query_ref, $to_extract_ref, $to_delete_ref) = @_;

	# Do the 1st BLAST (probe vs target)
	$self->search_target_file_using_blast($query_ref);

	# Create a non-redundant result set for this query
	my @combined;
	$self->created_combined_locus_set($query_ref, \@combined);

	# Defragment the combined set of new and previously identified loci
	$self->create_nonredundant_hit_set(\@combined, $to_extract_ref, $to_delete_ref);

}

#***************************************************************************
# Subroutine:  run_digs_classify_phase
# Description: extract hits and classify 
#***************************************************************************
sub run_digs_classify_phase {

	my ($self, $query_ref, $to_extract_ref, $to_delete_ref) = @_;

	# Get/create the various objects we will use here
	my $db_ref         = $self->{db};
	my $classify_obj   = Classify->new($self);
	my $extract_obj    = Extract->new($self);
	
	# Extract newly identified or extended sequences
	my @extracted;
	foreach my $locus_ref (@$to_extract_ref) {			
		$locus_ref->{target_path} = $query_ref->{target_path};
		$extract_obj->extract_locus_sequence_using_blast($locus_ref);	
		my %copy_locus = %$locus_ref;
		push (@extracted, \%copy_locus);
	}	

	# Do BLAST-based classification (new/extended hits from 1st search vs reference library)
	$self->classify_sequences_using_blast(\@extracted, $query_ref);

	# Update the digs_results table
	$db_ref->update_db($to_delete_ref, \@extracted, 'digs_results_table');

	# Update the searches_performed table
	my $searches_table = $db_ref->{searches_table};
	$searches_table->insert_row($query_ref);

}

#***************************************************************************
# Subroutine:  reassign
# Description: re-classify sequences in the digs_results_table 
#***************************************************************************
sub reassign {
	
	my ($self) = @_;

	# Get data structures, paths and flags from self
	my $blast_obj      = $self->{blast_obj};
	my $crossmatch_obj = $self->{crossmatch_obj};
	my $result_path    = $self->{report_dir};
	my $verbose        = $self->{verbose};
	my $classifier     = Classify->new($self);
	
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
		$classifier->classify_sequence_using_blast($locus_ref);
		
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
			$crossmatch_obj->update_cross_matching($previous_key, $assigned_key);
				
			# Insert the data
			my $where = " WHERE record_id = $record_id ";
			delete $locus_ref->{record_id}; # Required to remove this
			delete $locus_ref->{organism};  # Update not required for this field
			$digs_results_table->update($locus_ref, $where);
		}
	}
	
	# Cleanup
	my $output_dir = $self->{report_dir};
	my $command1 = "rm -rf $output_dir";
	system $command1;
}

############################################################################
# INTERNALS - search, classify, defragment, and updated screening database
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
	my $db_ref       = $self->{db};

	# Sanity checking
	unless ($min_length) { die; }
	unless ($min_score)  { die; }	
	unless ($blast_obj)  { die; } 
	unless ($tmp_path)   { $devtools->print_hash($self); die; } 

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
	my $completed = $self->{current_query_num};	
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
			$db_ref->insert_row_in_active_set_table($query_ref, $hit_ref);
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
# Subroutine:  classify_sequences_using_blast
# Description: classify a set of sequences using blast
#***************************************************************************
sub classify_sequences_using_blast {

	my ($self, $extracted_ref, $query_ref) = @_;

	my $verbose = $self->{verbose};
	my $classifier = Classify->new($self);
	my $assigned_count   = 0;
	my $crossmatch_count = 0;
	unless ($query_ref) { die; }
	foreach my $locus_ref (@$extracted_ref) { # Iterate through the matches

		# Classify using BLAST
		my $blast_alg = $classifier->classify_sequence_using_blast($locus_ref);
		my $assigned  = $locus_ref->{assigned_name};
		unless ($assigned) { die; }
		if ($assigned) { $assigned_count++; }

		# Get the unique key for this probe
		my $probe_name  = $query_ref->{probe_name};
		my $probe_gene  = $query_ref->{probe_gene};
		my $probe_key   = $probe_name . '_' . $probe_gene; 		

		# Record cross-matching
		if ($probe_key ne $assigned) {
			$crossmatch_count++;
			my $crossmatch_obj = $self->{crossmatch_obj};
			$crossmatch_obj->update_cross_matching($probe_key, $assigned);
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
# Subroutine:  created_combined_locus_set
# Description: combine a new set of hits with existing hits for a given target file
#***************************************************************************
sub created_combined_locus_set {
	
	my ($self, $query_ref, $combined_ref) = @_;

	# Get flags and objects
	my $verbose        = $self->{verbose};
	my $db_ref         = $self->{db};

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
	$db_ref->get_sorted_digs_results(\@digs_results, $where);
	my $num_loci = scalar @digs_results;
	if ($verbose) { print "\n\t\t # $num_loci previously extracted $probe_type loci"; }
		
	# Add the digs results to the BLAST hits in the active_set table
	$db_ref->add_digs_results_to_active_set(\@digs_results);

	# Get sorted list of digs results and BLAST hits from active_set table
	$db_ref->get_sorted_active_set($combined_ref, $where);

	# Show output 
	my $total_loci = scalar @$combined_ref;
	if ($verbose) {
		if ($total_loci > 0) {
			print "\n\t\t # $total_loci rows in active set (including $num_loci previously extracted) ";
		}
	}
}

#***************************************************************************
# Subroutine:  create_nonredundant_hit_set
# Description: combine most recent set of hits with previous 
#              and derive a non-redundant set
#***************************************************************************
sub create_nonredundant_hit_set {

	my ($self, $combined_ref, $to_extract_ref, $to_delete_ref) = @_;

	# Create the defragment object
	my $defragment_obj = Defragment->new($self);
		
	# Compose clusters of overlapping/adjacent loci		
	my %clusters;
	$defragment_obj->compose_clusters(\%clusters, $combined_ref);
	# DEBUG $devtools->print_hash(\%clusters); die;

	# Derive a non-redundant locus set (i.e. one merged locus for each cluster)
	my %merged;
	my %singletons;
	$defragment_obj->merge_clustered_loci(\%clusters, \%merged, \%singletons, $to_delete_ref);

	# Work out what we need to extract and classify 
	$defragment_obj->get_loci_to_extract(\%merged, \%singletons, $to_extract_ref);

}

############################################################################
# CONSOLE OUTPUT 
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
	$console->refresh();
	my $program_version = $self->{program_version};
	
    my $HELP   = "\n\n\t ### DIGS version $program_version";
       $HELP .= "\n\t ### usage: $0 m=[option] -i=[control file] -h=[help]\n";

       $HELP  .= "\n\t ### Main functions\n"; 
	   $HELP  .= "\n\t -m=1  Prepare target files (index files for BLAST)";		
	   $HELP  .= "\n\t -m=2  Do DIGS"; 
	   $HELP  .= "\n\t -m=3  Reassign loci"; 
	   $HELP  .= "\n\t -m=4  Defragment loci"; 
	   $HELP  .= "\n\t -m=5  Consolidate loci"; 

       $HELP  .= "\n\n\t ### Summarising target databases\n"; 	   
	   $HELP  .= "\n\t -g=1  Summarise targets (brief summary, by species)";
	   $HELP  .= "\n\t -g=2  Summarise targets (long, by individual target file)\n";

       $HELP  .= "\n\t ### Managing DIGS screening DBs\n"; 
	   $HELP  .= "\n\t -d=1  Import tab-delimited data"; 
	   $HELP  .= "\n\t -d=2  Flush core tables"; 
	   $HELP  .= "\n\t -d=3  Drop tables";
	   $HELP  .= "\n\t -d=4  Drop a screening DB"; 
	   $HELP  .= "\n\t -d=5  Append data to 'digs_results' table"; 
	   $HELP  .= "\n\t -d=6  Extract sequences using tabular file"; 

	   $HELP  .= "\n\n\t Target path variable '\$DIGS_GENOMES' is set to '$ENV{DIGS_GENOMES}'";
	   $HELP  .= "\n\n"; 

	print $HELP;
}

#***************************************************************************
# Subroutine:  show_digs_progress
# Description: show progress in DIGS screening
#***************************************************************************
sub show_digs_progress {

	my ($self) = @_;

	# Get the counts
	my $total_queries   = $self->{total_queries};
	my $completed       = $self->{current_query_num};	
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
	my $crossmatch_obj = $self->{crossmatch_obj};
	my $verbose = $self->{verbose};
	if ($verbose and $option eq 2
	or  $verbose and $option eq 3) { 
		$crossmatch_obj->show_cross_matching();
	}

	# Print finished message
	print "\n\n\t ### Process completed ~ + ~ + ~";

}

############################################################################
# EOF
############################################################################
