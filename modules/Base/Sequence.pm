#!/usr/bin/perl -w
############################################################################
# Script:       Sequence 
# Description:  Object for representing a sequence with common sequence
#               processing functions
# History:      Rob Gifford, Decemeber 2006: Creation
############################################################################
package Sequence;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

############################################################################
# Globals
############################################################################
1;

# Create hash with IUPAC codes for nucleic acids
my %iupac;
$iupac{A} = 'adenine';
$iupac{T} = 'thymine';
$iupac{G} = 'guanine';
$iupac{C} = 'cytosine';

# Twofold degenerate nucleotides
my @R = qw [ G A ];
$iupac{R} = \@R;
my @Y = qw [ T C ];
$iupac{Y} = \@Y;
my @K = qw [ G T ];
$iupac{K} = \@K;
my @M = qw [ A C ];
$iupac{M} = \@M;
my @S = qw [ G C ];
$iupac{S} = \@S;
my @W = qw [ A T ];
$iupac{W} = \@W;

# Threefold degenerate nucleotides
my @B = qw [ G T C ];
$iupac{B} = \@B;
my @D = qw [ G A T ];
$iupac{D} = \@D;
my @H = qw [ A C T ];
$iupac{H} = \@H;
my @V = qw [ G C A ];
$iupac{V} = 1;

# Fourfold degenerate nucleotides
my @N = qw [ A G C T ];
$iupac{N} = \@N;

my @dinucs = qw [ AA AG AC AT GA GG GC GT CA CG CC CT TA TG TC TT ];

# Create base objects
my $fileio    = FileIO->new();
my $devtools  = DevTools->new();
my $console   = Console->new();

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Parameters
#    Member variables:
#       sequence_table: reference to a table object for the sequence table
#***************************************************************************
sub new {

	my ($invocant, $sequence, $header, $id) = @_;
	
	my $class = ref($invocant) || $invocant;

	# Member variables
	my $self = {
		iupac         => \%iupac,
		dinucs        => \@dinucs,
		sequence      => $sequence,
		sequence_id   => $id,
		header        => $header,
	};

	bless ($self, $class);
	return $self;
}

############################################################################
# Manipulations of sequence string
############################################################################

#***************************************************************************
# Subroutine:  translate
# Description: Translate na sequence to aa
# Arguments:   $sequence: sequence to translate
# Arguments:   $frame to translate in: (possible values 1, 2, 3)
#***************************************************************************
sub translate {

	my ($self, $sequence, $frame) = @_;

	my @sequence = split ('', $sequence);
	
	# do frame adjustment
	if ($frame) {
		if    ($frame eq 1)  {  }
		elsif ($frame eq 2)  { shift (@sequence); }
		elsif ($frame eq 3)  { shift (@sequence);	
			                   shift (@sequence); }
		else {
			die "\n\t translation in frame '$frame' not implemented\n\n";
		}
	}

	# Translate the sequence
	my $codon;
	my $translation;
	my $codon_index;
	foreach my $nt (@sequence) {
		$codon .= $nt;
		$codon_index++;
		if ($codon_index eq 3) {
			if ($codon =~ /-/) {
				$translation .= '-';
			}
			# can't do degenerate codons
			elsif ($codon =~ 'N') {
				$translation .= '?';
			}
			else {
				my @codon_list;
				$self->get_codon_list($codon, \@codon_list);
				# if the codon is degenerate arbitrarily pick first one in the list
				my $list_codon = $codon_list[0];
				my $aa = $self->codon2aa($list_codon);
				$translation .= $aa;
			}
			$codon = '';
			$codon_index = 0;
		}
	}
	return $translation;
}

#***************************************************************************
# Subroutine:  reverse and complement
# Description: reverse and complement sequence
# Arguments:   $sequence: the sequence string to reverese and complement 
#***************************************************************************
sub reverse_and_complement {

	my ($self, $sequence) = @_;
    $sequence = uc $sequence;
	
	my @sequence = split ('', $sequence);

	@sequence = reverse(@sequence);

	my $rc_sequence;
	foreach my $nt (@sequence) {
		
		if    ($nt eq 'A') { $nt = 'T'; }
		elsif ($nt eq 'T') { $nt = 'A'; }
		elsif ($nt eq 'G') { $nt = 'C'; }
		elsif ($nt eq 'C') { $nt = 'G'; }
		elsif ($nt eq 'N') { $nt = 'N'  }
	
		else {
			print "\n\t ERROR: input was $sequence";
			die "\n\t This routine cannot yet deal with NT $nt\n"; 
		}
	
		$rc_sequence .= $nt;
	}
	return $rc_sequence;
}


#***************************************************************************
# Subroutine:  extract_subsequence
# Description: extract subseqeunce region from a given sequence 
# Arguments:   $sequence: the sequence string to extract from 
#              $start: region start position
#              $stop:  region start position
#***************************************************************************
sub extract_subsequence {
	
	my ($self, $sequence, $start, $stop) = @_;

	# Sanity checking
	unless ($start) { $start = 1; }
	unless ($sequence and $start and $stop) { 
		die;
	}
	my @chars = split ('', $sequence);
	my $subsequence = '';
	#print "\n\t coordinates $start and $stop";
	my $i = 0;
	foreach my $char (@chars) {
		$i++;
		if ($i >= $start and $i <= $stop) {
			$subsequence .= $char;
		}
	}
	return $subsequence;
}

#***************************************************************************
# Subroutine:  get translations 
# Description: Get all possible translations of a given codon
# Arguments:   $data_ref: reference to hash to store all distinct AAs that 
#              can be obtained by translating the (possibly degenerate) 
#              codon
#***************************************************************************
sub get_translations {

	my ($self, $codon) = @_;

	# translate the codon (to all possible if degenerate)
	my @codon_list;
	$self->get_codon_list($codon, \@codon_list);
	my $translations = undef;
	my %seen_translations;
	foreach my $possible_codon (@codon_list) {
	
		my $translation .= $self->codon2aa($possible_codon);
		unless ($seen_translations{$translation}) {
			$translations .= $translation;
			$seen_translations{$translation} = 1;
		}
	}

	return $translations;
}

#***************************************************************************
# Subroutine:  get codon list
# Description: get a list with all the possible codons from a degenerate 
#              codon
# Arguments:   $codon: the degenerate codon to deal with
#              $array_ref: reference to an array to store the codon list
#***************************************************************************
sub get_codon_list {

	my($self, $codon, $array_ref) = @_;

	# populate a hash of arrays with degeneracy codes
    my (%degeneracy_code) = (
    
		'A' => ['A'],
		'T' => ['T'],
		'C' => ['C'],
		'G' => ['G'],
		'-' => ['-'],
		'M' => ['A' , 'C'],
		'R' => ['A' , 'G'],
		'W' => ['A' , 'T'],
		'S' => ['G' , 'C'],
		'Y' => ['C' , 'T'],
		'K' => ['G' , 'T'],
		'H' => ['A' , 'C', 'T'],
		'V' => ['A' , 'G', 'C'],
		'D' => ['A' , 'G', 'T'],
		'B' => ['G' , 'T', 'C'],
		'N' => ['A' , 'G', 'C' , 'T'],
	);
	
	my @bases = split('', $codon);
	my $first  = $degeneracy_code{$bases[0]};
	my $second = $degeneracy_code{$bases[1]};
	my $third  = $degeneracy_code{$bases[2]};

	foreach my $first_pos (@$first) {

		foreach my $second_pos (@$second) {
		
			foreach my $third_pos (@$third) {
		
				my $codon_option = $first_pos . $second_pos . $third_pos;
				push (@$array_ref, $codon_option);
			}
		}
	}
}

#***************************************************************************
# Subroutine:  codon2aa
# Description: does translation accorsing to standard genetic code
#              (taken from the o'reilly book 'beginning perl for 
#              bioinformatics).
#***************************************************************************
sub codon2aa {
   
	my ($self, $codon) = @_;
    
	$codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if (exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }
	else {
		if ($codon =~ /[-~MRWSYKHVDBN]/) {  
			#print "\n\t Can't translate codon $codon";
			return '-';
		}	
		else {
			print "Bad codon '$codon'!!\n";
			return '-';
		}
	}
}

############################################################################
# END OF FILE
############################################################################
