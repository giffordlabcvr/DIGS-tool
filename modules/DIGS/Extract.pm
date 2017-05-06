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

# Program components
use DIGS::Initialise;    # Initialises the DIGS tool
use DIGS::ScreenBuilder; # To set up a DIGS run
use DIGS::Defragment;    # Cluster/defragment/consolidate tools

############################################################################
# Globals
############################################################################

# Base objects
my $fileio    = FileIO->new();
my $console   = Console->new();
my $devtools  = DevTools->new();
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
# Subroutine:  extract_sequences_using_blast
# Description: extract sequences from a  BLAST-indexed FASTA file, using BLAST 
#***************************************************************************
sub extract_sequences_using_blast {

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

############################################################################
# EOF
############################################################################