#!usr/bin/perl -w
############################################################################
# Module:      Extract.pm 
# Description: Functions for extracting sequences from FASTA files 
# History:     May 2017: Created by Robert Gifford 
############################################################################
package Extract;

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

my $maximum_extract_length = 30000;

1;

############################################################################
# MAIN LOOP
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create new Extract 'object'
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {
		
		# Global settings
		process_id             => $parameter_ref->{process_id},
		program_version        => $parameter_ref->{program_version},
		
		# Flags
		verbose                => $parameter_ref->{verbose},
		force                  => $parameter_ref->{force},
		
		# Paths 
		blast_bin_path         => $parameter_ref->{blast_bin_path},

		# Parameters for DIGS
		extract_buffer         => '',   # Obtained from control file

	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# INTERNAL FUNCTIONS: EXTRACT
############################################################################

#***************************************************************************
# Subroutine:  extract_locus_sequence_using_blast
# Description: extract sequence from a  BLAST-indexed FASTA file using BLAST 
#***************************************************************************
sub extract_locus_sequence_using_blast {

	my ($self, $locus_ref) = @_;

	# Interface to BLAST
	my %blast_params;
	$blast_params{blast_bin_path} = $self->{blast_bin_path};
	my $blast_obj = BLAST->new(\%blast_params);
	
	# Get paths, obhects and variables from self
	my $verbose   = $self->{verbose};
	my $buffer    = $self->{extract_buffer};

	# Add any buffer 
	if ($buffer) { 
		my $orientation = $locus_ref->{orientation};
		$self->add_buffer_to_sequence($locus_ref, $orientation); 
	}

	my $seq_length = ($locus_ref->{end} - $locus_ref->{start}) + 1;
	#$seq_length = length $sequence; # Set sequence length
	if ($seq_length > $maximum_extract_length) { 
	
		print "\n\t # Sequence too long - skipping!!"; 
		sleep 1;
		return;
	}

	# Extract the sequence
	my $target_path = $locus_ref->{target_path};
	my $sequence  = $blast_obj->extract_sequence($target_path, $locus_ref);
	if ($sequence) {
		
		# If we extracted a sequence, update the data for this locus
		if ($verbose) { print "\n\t\t\t # Extracted sequence: $seq_length nucleotides "; }
		$locus_ref->{extract_start}   = $locus_ref->{start};
		$locus_ref->{extract_end}     = $locus_ref->{end};
		$locus_ref->{sequence}        = $sequence;
		$locus_ref->{sequence_length} = $seq_length;
	}
}

#***************************************************************************
# Subroutine:  add_buffer_to_sequence
# Description: eadd leading-and-trailing buffer to extract coordinates
#***************************************************************************
sub add_buffer_to_sequence {

	my ($self, $locus_ref, $orientation) = @_;

	my $buffer = $self->{extract_buffer};
		
	if ($orientation eq '-') {
		$locus_ref->{start} = $locus_ref->{start} + $buffer;
		$locus_ref->{end}   = $locus_ref->{end} - $buffer;
		if ($locus_ref->{end} < 1) { # Don't allow negative coordinates
			$locus_ref->{end} = 1;
		}	
	}
	else {
		$locus_ref->{start} = $locus_ref->{start} - $buffer;
		if ($locus_ref->{start} < 1) { # Don't allow negative coordinates
			$locus_ref->{start} = 1;
		}	
		$locus_ref->{end} = $locus_ref->{end} + $buffer;
	}
}

############################################################################
# EOF
############################################################################
