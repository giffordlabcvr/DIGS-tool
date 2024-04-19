#!usr/bin/perl -w
############################################################################
# Module:      Classify.pm 
# Description: Capture information about cross-matching during DIGS
# History:     May 2017: Created by Robert Gifford 
############################################################################
package Classify;

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
# Description: create new Classify 'object'
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {
		
		# Flags
		verbose                => $parameter_ref->{verbose},

		# Parameters for reverse BLAST (hits versus reference sequence library)
		num_threads  => $parameter_ref->{rev_num_threads},
		word_size    => $parameter_ref->{rev_word_size},
		evalue       => $parameter_ref->{rev_evalue},
		penalty      => $parameter_ref->{rev_penalty},
		reward       => $parameter_ref->{rev_reward},
		gapopen      => $parameter_ref->{rev_gapopen},
		gapextend    => $parameter_ref->{rev_gapextend},
		dust         => $parameter_ref->{rev_dust},
		softmasking  => $parameter_ref->{rev_softmasking},
		seg          => $parameter_ref->{rev_seg},

		# Paths used in DIGS process
		tmp_path               => $parameter_ref->{tmp_path},
		blast_bin_path         => $parameter_ref->{blast_bin_path},
		# TODO: check why both these are neccessary
		aa_reference_library   => $parameter_ref->{aa_reference_library},
		na_reference_library   => $parameter_ref->{na_reference_library},
		blast_orf_lib_path     => $parameter_ref->{blast_orf_lib_path},
		blast_utr_lib_path     => $parameter_ref->{blast_utr_lib_path},

	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# INTERNAL FUNCTIONS: CLASSIFY
############################################################################

#***************************************************************************
# Subroutine:  classify_sequence_using_blast
# Description: classify a nucleotide sequence using blast 
#***************************************************************************
sub classify_sequence_using_blast {

	my ($self, $locus_ref) = @_;

	# Set up BLAST object with parameters for the reverse BLAST (from $self)
	my $blast_obj = BLAST->new($self);
	my $result_path = $self->{tmp_path};
	unless ($result_path) { die; } # Sanity checking
	
	# Get data about this probe sequence 
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

	# Set parameters for the reverse BLAST
	#$devtools->print_hash(\%blast_run_params);
	#$devtools->print_hash($self); die;
	# Do reverse BLAST search (hits versus reference sequence library)
	$blast_obj->blast($blast_alg, $lib_path, $query_file, $result_file, $self);
	
	# Parse the results
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
	my $success = 1;

	# Deal with a query that matched nothing in the 2nd BLAST search
	unless ($assigned_key) {
		$self->set_default_values_for_unassigned_locus($locus_ref);	
		$assigned = undef;
		$success = undef;
	}
	else {	# Assign the extracted sequence based on matches from 2nd BLAST search

		# Split assigned to into (i) refseq match (ii) refseq description (e.g. gene)	
		my @assigned_key  = split('_', $assigned_key);
		my $assigned_gene = pop @assigned_key;
		my $assigned_name = join ('_', @assigned_key);
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
		#$devtools->print_hash($locus_ref); die;
	
		my $id = $locus_ref->{record_id};
		print "\n\t\t# Assigned as '$assigned_name ($assigned_gene)'";
		print " via $blast_alg comparison to $lib_file";
		$assigned = $assigned_name . '_' . $assigned_gene;
	}

	# Clean up
	my $command1 = "rm $query_file";
	my $command2 = "rm $result_file";
	system $command1;
	system $command2;
	
	return $success;
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
			$devtools->print_hash($self); 
			die "\n\t NO ORF LIBRARY defined";
		}
	}	
	return $lib_path;
}

############################################################################
# EOF
############################################################################
