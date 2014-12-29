#!/usr/bin/perl -w
############################################################################
# Script:       Mutation 
# Description:  Object for representing a reference sequence, primarily 
#               for use when aligning a set of nucleic acid sequences to an
#               amino acid reference sequence using the LAP program. 
# History:      Rob Gifford, Novemeber 2006: Creation
############################################################################
package Mutation;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base Classes
use Base::FileIO;
use Base::DevTools;

############################################################################
# Globals
############################################################################
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create a new reference sequence object 
# Arguments:   $input_ref: reference to an array containing input file data
#              assumed to be in standard reference sequence format
#***************************************************************************
sub new {

	my ($invocant, $parameters) = @_;
	my $class = ref($invocant) || $invocant;
	
	# set optional parameters to NULL or collect values if given
	my $mutation_class;
	my $sequence_id;
	if ($parameters->{class}) { 
		$mutation_class = $parameters->{class}; 
	}
	else { 
		$mutation_class = 'NULL';
	}
	if ($parameters->{sequence_id}) { 
		$sequence_id = $parameters->{sequence_id};
	}
	else { 
		$sequence_id = 'NULL';
	}
	
	# Member variables
	my $self = {
		sequence_id       => $sequence_id,
		class             => $mutation_class,
		gene              => $parameters->{gene},
		refseq_aa         => $parameters->{refseq_aa},
		position          => $parameters->{position},
		aa                => $parameters->{aa},
		codon             => $parameters->{codon},
		ref_codon         => $parameters->{ref_codon},
		num_transitions   => $parameters->{num_transitions},
		num_transitions   => $parameters->{num_transversions},
		num_changes       => $parameters->{num_changes},
		num_syn           => $parameters->{num_syn},
		num_nonsyn        => $parameters->{num_nonsyn},
		'A->T'            => $parameters->{'A->T'},
		'A->C'            => $parameters->{'A->C'},
		'A->G'            => $parameters->{'A->G'},
		'T->A'            => $parameters->{'T->A'},
		'T->C'            => $parameters->{'T->C'},
		'T->G'            => $parameters->{'T->G'},
		'C->A'            => $parameters->{'C->A'},
		'C->T'            => $parameters->{'C->T'},
		'C->G'            => $parameters->{'C->G'},
		'G->A'            => $parameters->{'G->A'},
		'G->T'            => $parameters->{'G->T'},
		'G->C'            => $parameters->{'G->C'},
	};

	bless ($self, $class);
	return $self;
}

############################################################################
# Public Member Functions
############################################################################

#***************************************************************************
# Subroutine:   create mutation key 
# Description: 
# Arguments:   
#***************************************************************************
sub create_mutation_key {

	my ($self) = @_;	

	my @key_elements;
	push(@key_elements, $self->{gene});
	push(@key_elements, $self->{refseq_aa});
	push(@key_elements, $self->{position});
	push(@key_elements, $self->{aa});

	my $mutation_key = join(':', @key_elements); 

	return $mutation_key;
}

#***************************************************************************
# Subroutine:   create position key 
# Description: 
# Arguments:   
#***************************************************************************
sub create_position_key {

	my ($self) = @_;	

	my @key_elements;
	push(@key_elements, $self->{gene});
	push(@key_elements, $self->{position});

	my $mutation_key = join(':', @key_elements); 

	return $mutation_key;
}


#***************************************************************************
# Subroutine:   increment_frequency
# Description: 
# Arguments:   
#***************************************************************************
sub increment_frequency {

	my ($self) = @_;	

	if ($self->{frequency}) {
		$self->{frequency} = $self->{frequency} + 1;
	}
	else {
		$self->{frequency} = 1;
	}
}


#***************************************************************************
# Subroutine:   describe
# Description:  
#***************************************************************************
sub describe {

	my ($self) = @_;	

	# get member variables
	my $gene      = $self->{gene};	
	my $refseq_aa = $self->{refseq_aa};	
	my $position  = $self->{position};	
	my $aa        = $self->{aa};	
	my $class     = $self->{class};	
	my $codon     = $self->{codon};	
	
	print "\n\n";	
	print "\n\t ### Describing Mutation object!!";
	print "\n\t ### Gene:      $gene";
	print "\n\t ### RefSeq AA: $refseq_aa";
	print "\n\t ### Position:  $position";
	print "\n\t ### AA:        $aa";
	print "\n\t ### Class:     $class";
	print "\n\t ### Codon:     $codon";

	if ($self->{seqset}) {
		print "\n\t ### Seqset:    $self->{seqset}";
	}
	if ($self->{naive}) {
		print "\n\t ### Naive:     $self->{naive}";
	}
	if ($self->{treated}) {
		print "\n\t ### Treated:   $self->{treated}";
	}

}



