#!/usr/bin/perl -w
############################################################################
# Script:       RefSeq 
# Description:  Class for representing a reference sequence
# History:      Rob Gifford, Novemeber 2010: Creation
############################################################################
package RefSeq;

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
use Base::BioIO;
use Base::Sequence;

# Component Classes
use GLUE::Mutation;
use GLUE::RefSeqParser;

############################################################################
# Globals
############################################################################
my $bioio    = BioIO->new();
my $fileio   = FileIO->new();
my $devtools = DevTools->new();
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

	my ($invocant, $params) = @_;
	my $class = ref($invocant) || $invocant;

	# DEBUG
	#$devtools->print_hash($params); die;

	# Member variables
	my $max_mixtures   = 4;
	my $self = {
		
		# Basic member variables
		sequence       => $params->{sequence},
		seq_len        => $params->{seq_len},
		start          => $params->{start},
		stop           => $params->{stop},
		genes          => $params->{genes},
		hash_genes     => $params->{hash_genes},
		features       => $params->{features},
		indexed_seq    => $params->{indexed_seq},
		hash_lists     => $params->{hash_lists},
		listnames      => $params->{listnames},
		typical        => $params->{typical},
		
		# Flags/settings
		max_mixtures   => $max_mixtures,

	};

	# Store the metadata
	my $metadata = $params->{metadata};
	my @keys = keys %$metadata;
	foreach my $key (@keys) {
		my $value = $metadata->{$key};
		$self->{$key} = $value;
	}

	bless ($self, $class);
	return $self;
}

############################################################################
# Public Member Functions
############################################################################

#***************************************************************************
# Subroutine:  compare to aligned sequence
# Description: determine what genes and mutations are present in a 
#              sequence that is profile aligned to this reference sequence
# Arguments:   $aln_seq_ref: aligned sequence in standard VGLUE object
#***************************************************************************
sub compare_to_aligned_sequence {

	my ($self, $aln_seq_ref) = @_;
	
	# Create variables for tracking progress and storing data
	my @syn_mutations;  
	my @nonsyn_mutations;  
	my %syn_mutations;  
	my %nonsyn_mutations;  
	my %genes;

	# Set up the data structures with the aligned sequence
	$aln_seq_ref->{syn_mutations}         = \@syn_mutations;
	$aln_seq_ref->{nonsyn_mutations}      = \@nonsyn_mutations;
	$aln_seq_ref->{hash_syn_mutations}    = \%syn_mutations;
	$aln_seq_ref->{hash_nonsyn_mutations} = \%nonsyn_mutations;
	$aln_seq_ref->{gene_data}             = \%genes;
	$aln_seq_ref->{num_changes}           = 0;
	$aln_seq_ref->{num_syn}               = 0;
	$aln_seq_ref->{num_nonsyn}            = 0;
	$aln_seq_ref->{num_transitions}       = 0;
	$aln_seq_ref->{num_transversions}     = 0;

	# Create indexed sequence
	#my $position = $aln_seq_ref->{aln_start};
	my $position = $self->{start};
	my $sequence = $aln_seq_ref->{sequence};
	unless ($position and $sequence) {
		print "\n\t '$position' and '$sequence'";
		die;
	}
	unless ($position) { $position=1; }
	my %indexed_seq;
	my @sequence = split ('', $sequence);
	foreach my $nt (@sequence) {
		$indexed_seq{$position} = $nt;
		$position++;
	}
	$aln_seq_ref->{indexed_seq} = \%indexed_seq;
	
	# Get mutations in each gene relative to the reference sequence
	$self->get_mutations($aln_seq_ref);
	
	# Work out the proportion of each gene present
	$self->get_genes($aln_seq_ref);
	#$devtools->print_hash($aln_seq_ref); die;
}

#***************************************************************************
# Subroutine:  get mutations 
# Description: get mutations relative to the reference sequence
#***************************************************************************
sub get_mutations {

	my ($self, $aln_seq_ref) = @_;

	# Create objects/get data from self
	my $seq_obj     = Sequence->new(); # Create sequence utility class
	#my $sequence    = $aln_seq_ref->{sequence};
	my $header = $aln_seq_ref->{header};

	my $sequence    = $aln_seq_ref->{sequence};
	#my $aln_start   = $aln_seq_ref->{aln_start};
	#my $aln_stop    = $aln_seq_ref->{aln_stop};
	my $aln_start   = $self->{start};
	my $aln_stop    = $self->{stop};
	unless ($aln_start and $aln_stop) { die; }

	unless ($aln_start and $aln_stop)  { die; }
	my $indexed_seq = $aln_seq_ref->{indexed_seq};
	unless ($sequence)  { die; }
	unless ($aln_start) { $aln_start=1; }
	#unless ($sequence and $aln_start) { die; }
	
	# Iterate through the nt positions
	my $position = $aln_start - 1;
	
	my @sequence = split ('',$sequence);

	#if ($header=~ /JX106336/) {
	#	print "<BR><BR> SEQUENCE  $sequence";
	#	print "<BR><BR> ALN START $aln_start";
	#	print "<BR><BR> ALN START $aln_stop <BR><BR>";
	#	my $length = scalar @sequence;
	#	print "<BR> SEQ LENGTH ($length)";
	#}
	foreach my $nt (@sequence) {
		
		$position++;
		#if ($header=~ /JX106336/) {
		#	print "<BR> in nt processing loop at postion ($position)";
		#}
		unless ($position >= 3) { next; } 
		
		# Get alignment data at this nucleotide position
		my $gene_codon  = $self->check_for_inframe_codon($position);
		unless ($gene_codon)   { next; }
		#if ($header=~ /JX106336/) {
		#	print "<BR> got inframe codon ($position)";
		#}
		# Get the codon for this amino acid position
		my $nt1 = $indexed_seq->{$position - 2};
		my $nt2 = $indexed_seq->{$position - 1};
		my $nt3 = $indexed_seq->{$position};
		my $codon;

		#if ($header=~ /JX106336/) {
		#	print "<BR> nucs $nt1 + $nt2 + $nt3";
		#}
		if ($nt1 and $nt2 and $nt3) {
			$codon = $nt1 . $nt2 . $nt3;
		}
		else { next; }
		#$devtools->print_hash($indexed_seq); die;	

		my $valid_codon = $seq_obj->validate_codon($codon);
		unless ($valid_codon) { next; }
		
		# Record a mutation
		#if ($header=~ /JX106336/) {
		#	print "<BR> codon valid: entering record mutation loop at postion ($position)";
		#}
		$self->record_mutation($aln_seq_ref, $codon, $position);
		
		# Record gene codon
		$self->record_gene_codon($aln_seq_ref, $codon, $gene_codon);
	}
	#$devtools->print_hash($aln_seq_ref);
	#die;
}
	
#***************************************************************************
# Subroutine:  check_for_inframe_codon 
# Description: work out if a nucleotide position corresponds to the last
#              nucleotide of a codon that is in-frame for a gene in this
#              reference sequence
#***************************************************************************
sub check_for_inframe_codon {

	my ($self, $codon_end_position) = @_;

	# Create the sequence object
	my $seq_obj = Sequence->new();
	
	# Iterate through the genes in the reference sequence
	my $genes = $self->{genes};
	my $gene_codon = undef;
	foreach my $gene (@$genes) {
		
		# Get gene data
		my $name  = $gene->{name};
		my $start = $gene->{start};
		my $stop  = $gene->{stop};
		
		# Check if we're in range for for this particular gene
		if ($codon_end_position >=$start and $codon_end_position <= $stop) {
		
			# Check if we're in frame
			my $relative_position = $codon_end_position - ($start);
			my $remainder = ($relative_position % 3);
			if ($remainder eq 0) { # In frame - return codon coordinates
				my $aa_pos = $relative_position / 3;
				unless ($aa_pos eq 0) {	
					$gene_codon = $name . ':' . $aa_pos; 
				}
			}
		}
		
		# DEBUG
		#print "\n\t inframe for $name";	
		#print "\n\t align seq position $codon_end_position";
		#print "\n\t relative position $relative_position";
	}
	return $gene_codon;
}

#***************************************************************************
# Subroutine:  record_mutation
# Description: compare codon at given position to reference sequence,
#              and record any amino acid mutation(s). There may be more 
#              than one mutation at a single position due to nucleotide
#              degeneracy,
#***************************************************************************
sub record_mutation {

	my ($self, $align_seq_ref, $codon, $position) = @_;

	# Get the residue at in the reference sequence at the given position
	my $seq_obj            = Sequence->new();
	my $genes              = $self->{genes};
	my $indexed_na_seq     = $self->{indexed_seq};
	my $max_mixtures       = $self->{max_mixtures};
	my $exclude            = $self->{exclude};
	my $mutations_ref      = $align_seq_ref->{nonsyn_mutations};
	my $hash_mutations     = $align_seq_ref->{hash_nonsyn_mutations};
	my $syn_mutations_ref  = $align_seq_ref->{syn_mutations};
	my $hash_syn_mutations = $align_seq_ref->{hash_syn_mutations};
	
	# DEBUG
	#$devtools->print_hash($indexed_na_seq);
	#print "<BR> ###########  RECORD MUTATION routine position $position\n\n\n";
	#my $indexed_aln_seq = $align_seq_ref->{indexed_seq};
	#$devtools->print_hash($indexed_aln_seq);
	#die;

	# Get refseq codon and AA
	my $ref_codon .= $indexed_na_seq->{$position - 2};
	   $ref_codon .= $indexed_na_seq->{$position - 1};
	   $ref_codon .= $indexed_na_seq->{$position};
	my $ref_aa     = $seq_obj->get_translations($ref_codon);

	# Get the amino acid(s) for the query seq
	my $translations = $seq_obj->get_translations($codon);
	unless ($translations) {  die; }
	my $num_translations = length $translations;

	# Create data structure for the mutation 
	my %mutation;
	$mutation{refseq_aa}       = $ref_aa;
	$mutation{ref_codon}       = $ref_codon;
	$mutation{codon}           = $codon;
	$mutation{translations}    = $num_translations;
	$self->set_position_coordinates($position, \%mutation);
	$seq_obj->compare_codons($ref_codon, $codon, \%mutation); # Compare codons

		
	#print "<BR> MUTATION OUTPUT";
	#$devtools->print_hash(\%mutation);
	#my $pos      = $mutation{position};
	#my $gene     = $mutation{gene};
	#print "<BR> GENE $gene; AA POSITION $pos";

	if ($translations eq $ref_aa) { 
		
		# If there are no coding differences look for synonymous diffs
		if ($codon ne $ref_codon) {
			
			$mutation{aa} = $translations;
			my $mutation_obj = Mutation->new(\%mutation);
			$mutation_obj->{ref_codon} = $ref_codon;
			$mutation_obj->{codon}     = $codon;	
			my $mutation_key = $mutation_obj->create_mutation_key();
			push(@$syn_mutations_ref, $mutation_obj);
			$hash_syn_mutations->{$mutation_key} = $mutation_obj;
			#return; 
		}
	}
	else {

		# Apply degeneracy rules	
		my @translations;
		if ($num_translations > $max_mixtures) {
			push (@translations, 'X');
		}
		else {
			@translations = split('', $translations);
		}
				

		# Record each possible translation as a separate mutation
		foreach my $aa (@translations) {
			if ($aa eq $ref_aa) { next; } # skip AA that matches reference sequence
			my $mutation_obj = Mutation->new(\%mutation);
			$mutation_obj->{aa} = $aa;
			
			# For excluding regions that are known to be problematic
			my $position_key = $mutation_obj->create_position_key();
			if ($exclude->{$position_key}) { next; }
			
			# Store the mutation
			my $mutation_key = $mutation_obj->create_mutation_key();
			#print "<BR> ###########  RECORDING MUTATION $mutation_key\n\n\n";
			push(@$mutations_ref, $mutation_obj);
			$hash_mutations->{$mutation_key} = $mutation_obj;
		}
	}
	#print "<BR><BR> ############# <BR>";
	
	# Add to the totals for the sequence
	$align_seq_ref->{num_syn}           = $align_seq_ref->{num_syn} + $mutation{num_syn};
	$align_seq_ref->{num_nonsyn}        = $align_seq_ref->{num_nonsyn} + $mutation{num_nonsyn};
	$align_seq_ref->{num_changes}       = $align_seq_ref->{num_changes} + $mutation{num_changes};
	$align_seq_ref->{num_transitions}   = $align_seq_ref->{num_transitions} + $mutation{num_transitions};
	$align_seq_ref->{num_transversions} = $align_seq_ref->{num_transversions} + $mutation{num_transversions};
	

	#$devtools->print_hash($align_seq_ref); die;
	#return 1;
}


#***************************************************************************
# Subroutine:  set_position_coordinates 
# Description: use reference sequence object to get relative coordinates
#              from position in entire aligned sequence.
#***************************************************************************
sub set_position_coordinates {

	my ($self, $position, $data_ref) = @_;

	# Get gene details from self
	my $genes = $self->{genes};
	
	my $located;
	my $highest_stop = 0;
	my $gene_position;
	foreach my $gene (@$genes) {
		
		my $name  = $gene->{name};
		my $start = $gene->{start};
		my $stop  = $gene->{stop};
		if ($stop > $highest_stop) { $highest_stop = $stop; }	
		
		# Check if we're in range
		if ($position >=$start and $position <= $stop) {
		
			# Check if we're in frame
			my $relative_position = $position - ($start);
			#print "\n\t relative: '$relative_position' = pos '$position' - (start '$start' - 1)";
			my $remainder = ($relative_position % 3);
			if ($remainder eq 0) {
				$gene_position = ($relative_position / 3); # this is how you adjust position
				if ($gene_position eq 0) {
					#print "\n\t gene '$name', start '$start', stop '$stop' pos '$position'";
					#print "relative $relative_position' remainder '$remainder' ";
					#print "gene pos '$gene_position'";
				}
				else {
					$data_ref->{gene}      = $name; 
					$data_ref->{position}  = $gene_position; 
					$located = 'true';
				}
			}
		}
	}
	unless ($located) {
		#if ($position > $highest_stop) {	
		#	die "position $position out of range for reference sequence"; 
		#}
		#else {
		#	die "position $position out of frame for reference sequence"; 
		#}
	}
	
	# DEBUG
	#print "\n\t position $position";
	#$devtools->print_hash($data_ref);

}

#***************************************************************************
# Subroutine:  record_gene_codon 
# Description: record a codon in a gene
#***************************************************************************
sub record_gene_codon {

	my ($self, $aln_seq_ref, $codon, $gene_codon) = @_;

	# Work out gene proportions in sequence
	my $aln_seq_genes_ref  = $aln_seq_ref->{gene_data};
	
	my @gene_codon = split(':', $gene_codon);
	my $gene_name       = $gene_codon[0];
	my $gene_position   = $gene_codon[1];
	
	
	# Record valid codon
	if ($aln_seq_genes_ref->{$gene_name}) {
		my $gene_data_ref   = $aln_seq_genes_ref->{$gene_name};
		if ($gene_data_ref->{codons}) {
			my $alnseq_gene_data_ref = $gene_data_ref->{codons};
			$alnseq_gene_data_ref->{$gene_codon} = $codon; 
		}
		else {
			my %alnseq_gene_data;
			$alnseq_gene_data{$gene_codon} = $codon; 
			$gene_data_ref->{codons} = \%alnseq_gene_data;
		}
	}
	else {
		my %alnseq_gene_data;
		$alnseq_gene_data{$gene_codon} = $codon; 
		my %gene_data;
		$gene_data{codons} = \%alnseq_gene_data;
		$aln_seq_genes_ref->{$gene_name} = \%gene_data;
	}
}

#***************************************************************************
# Subroutine:  get_genes
# Description: get genes
#***************************************************************************
sub get_genes {

	my ($self, $aln_seq_ref) = @_;

	# Work out gene proportions in sequence
	my %ref_genes;
	my $ref_genes = $self->{genes};
	foreach my $ref_gene (@$ref_genes) {
		my $gene_name = $ref_gene->{name};
		$ref_genes{$gene_name} = $ref_gene;
	}
	
	# Work out how much of each gene present
	my $aln_genes_ref = $aln_seq_ref->{gene_data}; 
	my $gene_name;
	my $gene_data_ref;
	while (( $gene_name, $gene_data_ref ) = each %$aln_genes_ref) {
		
		# Get reference gene length
		my $ref_gene  = $ref_genes{$gene_name};
		my $ref_start = $ref_gene->{coding_start};
		my $ref_stop  = $ref_gene->{coding_stop};
		unless ($ref_start and $ref_stop) {
			die "gene data not properly initialised";
		}
		my $ref_gene_length = $ref_stop - ($ref_start - 1);

		# Work out what proportion of the gene we have
		my $codon_ref = $gene_data_ref->{codons};
		my @codons = keys %$codon_ref;
		my $num_codons = scalar @codons;
		#print "<BR> Number of codons $num_codons";
		
		my $gene_prop   = $num_codons / $ref_gene_length;
		my $f_gene_prop = sprintf("%.1f", $gene_prop);
		my $minimum    = 0.1;
		if ($f_gene_prop >= $minimum) {
			$gene_data_ref->{minimum_present} = 1;
		}
		$gene_data_ref->{gene_length} = $num_codons;
		
		#print "\n\t $f_gene_prop = coverage of $gene_length / $ref_gene_length";
		#my $alnseq_gene_data  = $gene_data_ref->{codons};
		#my @num_valid  = keys %$alnseq_gene_data;
		#my $num_valid = scalar @num_valid;
		#my $gene_prop = $num_valid / $aa_gene_length;
		#print "\n\t $gene_prop = coverage of $num_valid / gene length  $aa_gene_length";
		## Record the formatted results
		#my $coverage  = sprintf("%.1f", $gene_prop);
		#$gene_data_ref->{coverage} = $coverage;
		
		# Record the gene as present in this sequence if > minimum present
		#my $minimum    = $gene->{minimum};
	}

}

############################################################################
# RETRIEVE FEATURES
############################################################################

#***************************************************************************
# Subroutine:  get_gene_names 
# Description: get gene names in an array
#***************************************************************************
sub get_gene_names {

	my ($self, $array_ref) = @_;

	my $genes_ref  = $self->{genes};
	foreach my $gene_ref (@$genes_ref) {
		my $gene_name = $gene_ref->{name};
		push(@$array_ref, $gene_name);
	}
}

#***************************************************************************
# Subroutine:  get_gene_coordinate_string
# Description:
#***************************************************************************
sub get_gene_coordinate_string {

	my ($self, $target_gene_name) = @_;

	my @string;
	my $genes_ref = $self->{genes};
	foreach my $gene_ref (@$genes_ref) {

		my $gene_name  = $gene_ref->{name};
		unless ($gene_name eq $target_gene_name) {
			next;
		} 

		# Iterate through the exons
		my $exon = 0;
		my $starts = $gene_ref->{starts};
		my $exons  = $gene_ref->{exons};
		
		foreach my $exon_start (@$starts) {
			$exon++;
			my $exon_stop = $exons->{$exon_start};
			push (@string, $exon_start);
			push (@string, $exon_stop);
		}
	}
	
	my $string = join(',', @string);
	return $string;
}

#***************************************************************************
# Subroutine:  get_orfs 
# Description: get nucleic acid sequences of ORFs
#***************************************************************************
sub get_orfs {

	my ($self, $orf_ref) = @_;

	my $seq_obj  = Sequence->new();
	my $sequence  = $self->{sequence};
	my $genes_ref = $self->{genes};
	foreach my $gene_ref (@$genes_ref) {

		my $gene_name  = $gene_ref->{name};
		my $gene_start = $gene_ref->{start};
		my $gene_stop  = $gene_ref->{stop};
		my $gene_coding_start = $gene_ref->{coding_start};
		my $gene_coding_stop  = $gene_ref->{coding_stop};
		my $gene_seq = $seq_obj->extract_subsequence($sequence, $gene_start, $gene_stop);

		# Iterate through the exons
		my $exon = 0;
		my $starts = $gene_ref->{starts};
		my $exons  = $gene_ref->{exons};
		my $orf = '';
		foreach my $exon_start (@$starts) {
			$exon++;
			my $exon_stop = $exons->{$exon_start};
			
			#print "\n\t###\t\t exon $exon: spans $exon_start-$exon_stop"; 
			if ($exon_start < $gene_coding_start) {
				$exon_start = $gene_coding_start;
			}
			if ($exon_stop > $gene_coding_stop) {
				$exon_stop = $gene_coding_stop;
			}
			#print "\n\t###\t\t exon coding $exon_start-$exon_stop"; 
			
			my $exon_seq = $seq_obj->extract_subsequence($sequence, $exon_start, $exon_stop);
			$orf .= $exon_seq;
		}
		my $f_gene_name = $gene_name;
		$orf_ref->{$gene_name} = $orf;
	}
}

#***************************************************************************
# Subroutine:  get_translated_orfs   
# Description:
#***************************************************************************
sub get_translated_orfs {

	my ($self, $orf_ref) = @_;

	my $seq_obj  = Sequence->new();
	my $sequence  = $self->{sequence};
	my $genes_ref = $self->{genes};
	foreach my $gene_ref (@$genes_ref) {

		my $gene_name  = $gene_ref->{name};
		my $gene_start = $gene_ref->{start};
		my $gene_stop  = $gene_ref->{stop};
		my $gene_coding_start = $gene_ref->{coding_start};
		my $gene_coding_stop  = $gene_ref->{coding_stop};
		my $gene_seq = $seq_obj->extract_subsequence($sequence, $gene_start, $gene_stop);

		# Iterate through the exons
		my $exon = 0;
		my $starts = $gene_ref->{starts};
		my $exons  = $gene_ref->{exons};
		my $orf = '';
		foreach my $exon_start (@$starts) {
			$exon++;
			my $exon_stop = $exons->{$exon_start};
			
			#print "\n\t###\t\t exon $exon: spans $exon_start-$exon_stop"; 
			if ($exon_start < $gene_coding_start) {
				$exon_start = $gene_coding_start;
			}
			if ($exon_stop > $gene_coding_stop) {
				$exon_stop = $gene_coding_stop;
			}
			#print "\n\t###\t\t exon coding $exon_start-$exon_stop"; 
			
			my $exon_seq = $seq_obj->extract_subsequence($sequence, $exon_start, $exon_stop);
			my $exon_orf = $seq_obj->translate($exon_seq);
			$orf .= $exon_orf;
		}
		# format gene name (consistent)
		#my $f_gene_name;
		#my @name = split ($gene_name);
		#my $i = undef;
		#foreach my $char (@name) {
		#	if ($i) { $f_gene_name  = uc $char; }
		#	else    { $f_gene_name .= lc $char; }
		#	$i++;
		#}

		#$orf_ref->{$f_gene_name} = $orf;
		$orf_ref->{$gene_name} = $orf;
	}
}

#***************************************************************************
# Subroutine:  get_utrs
# Description: get the untranslated regions for this refseq
#***************************************************************************
sub get_utrs {

	my ($self, $utr_ref) = @_;

	my $seq_obj  = Sequence->new();
	my $sequence     = $self->{sequence};
	my $features_ref = $self->{features};
	foreach my $feature_ref (@$features_ref) {
		
		#$devtools->print_hash($feature_ref); die;
		my $gene_name  = $feature_ref->{name};
		my $gene_start = $feature_ref->{start};
		my $gene_stop  = $feature_ref->{stop};
	
		# Iterate through the exons
		my $exon = 0;
		my $starts = $feature_ref->{starts};
		my $exons  = $feature_ref->{exons};
		my $utr = '';
		foreach my $exon_start (@$starts) {
			$exon++;
			my $exon_stop = $exons->{$exon_start};
			#print "\n\t###\t\t exon coding $exon_start-$exon_stop"; 
			my $exon_seq = $seq_obj->extract_subsequence($sequence, $exon_start, $exon_stop);
			$utr .= $exon_seq;
		}
		#print "$utr";
		my $f_gene_name = uc $gene_name;
		$utr_ref->{$f_gene_name} = $utr;
	}
	#$devtools->print_array($features_ref); die;
	#$devtools->print_hash($utr_ref); die;
}

#***************************************************************************
# Subroutine:  append_new_data
# Description: 
#***************************************************************************
sub append_new_data {

	my ($self, $data_ref) = @_;

	# Store the metadata
	#my $metadata = $data_ref->{metadata};
	#my @keys = keys %$metadata;
	#foreach my $key (@keys) {
	#	my $value = $metadata->{$key};
	#	print "<BR> Adding $key : $value";
	#	$self->{$key} = $value;
	#}
	#$devtools->print_web_hash($data_ref); die;

	# Store the metadata
	my $genes_ref = $self->{genes};
	my $new_genes_ref = $data_ref->{genes};
	foreach my $gene_ref (@$new_genes_ref) {
		#my $value = $new_genes_ref->{$gene};
		#$devtools->print_web_hash($gene_ref); die;
		push (@$genes_ref, $gene_ref);
	}
	#$devtools->print_web_hash($data_ref); die;
	#print "<BR> Adding $gene";
	#$self->{$key} = $value;

	# Basic member variables
	my $coordinates_ref = $data_ref->{coordinates};
	my $group           = $data_ref->{virus_group};
	my $subgroup        = $data_ref->{virus_subgroup};
	my $family          = $data_ref->{family};
	my $genome_type     = $data_ref->{genome_type};
	my $accession       = $data_ref->{accession};
	my $refseq_name     = $data_ref->{name};
	my $genes           = $data_ref->{genes};
	my $hash_genes      = $data_ref->{hash_genes};
	my $features        = $data_ref->{features};
	my $indexed_seq     = $data_ref->{indexed_seq};
	my $hash_lists      = $data_ref->{hash_lists};
	my $listnames       = $data_ref->{listnames};
	my $typical         = $data_ref->{typical};
	#print "<BR> WAGA CONTROL EXIT FOR DEVELOPMENT"; die;
}

############################################################################
# WRITING & OUTPUT FXNS 
############################################################################

#***************************************************************************
# DEVELOPMENT
# Subroutine:  write_self_to_text
# Description: create the standard VGLUE text file
#***************************************************************************
sub write_self_to_text {

	my ($self, $directory_path) = @_;

	my @file;  # array to store the outfile
	
	# Do metadata section
	my $refseq_name       = $self->{name};
	my $full_name         = $self->{full_name};
	my $virus_group       = $self->{virus_group};
	my $virus_genus       = $self->{virus_genus};
	my $virus_subgroup    = $self->{virus_subgroup};
	my $virus_class       = $self->{virus_class};
	my $virus_supertribe  = $self->{virus_supertribe};
	my $virus_tribe       = $self->{virus_tribe};
	my $virus_state       = $self->{virus_state};
	my $genome_coverage   = $self->{genome_coverage};
	my $genome_state      = $self->{genome_state};
	my $virus_subfamily   = $self->{virus_subfamily};
	my $sequence_type     = $self->{sequence_type};
	unless ($virus_subfamily) { die; }

	push (@file, "Begin Metadata;");
	push (@file, "\nname             = $refseq_name;");
	push (@file, "\nfull_name        = $full_name;");
	push (@file, "\nvirus_family     = Retroviridae;");
	push (@file, "\nvirus_subfamily  = $virus_subfamily;");
	push (@file, "\nvirus_supertribe = $virus_supertribe;");
	push (@file, "\nvirus_tribe      = $virus_tribe;");
	push (@file, "\nvirus_genus      = $virus_genus;");
	push (@file, "\nvirus_subgroup   = $virus_subgroup;");
	if ($self->{genome_coverage}) {
		my $genome_coverage    = $self->{genome_coverage};
		push (@file, "\ngenome_coverage  = $genome_coverage;");
	}
	if ($self->{genome_state}) {
		my $genome_state    = $self->{genome_state};
		push (@file, "\ngenome_state     = $genome_state;");
	}
	unless ($self->{sequence_type}) { die; }
	push (@file, "\nsequence_type    = $sequence_type;");

	if ($self->{host_sci_name}) {
		my $host_sci_name    = $self->{host_sci_name};
		push (@file, "\nhost_sci_name    = $host_sci_name;");
	}
	if ($self->{host_common_name}) {
		my $host_common_name    = $self->{host_common_name};
		push (@file, "\nhost_common_name = $host_common_name;");
	}
	if ($self->{accession}) {
		my $accession    = $self->{accession};
		push (@file, "\naccession        = $accession;");
	}

	push (@file, "\nEndblock;");

	{
		# Do features section
		my $features_ref   = $self->{features};	
		push (@file, "\n\nBegin Features;");
		my $feature_ref = shift @$features_ref;	
		my $ftype      = $feature_ref->{type};
		my $fname      = $feature_ref->{name};
		my $ffull_name = $feature_ref->{full_name};
		if ($ffull_name eq 'leader') { $ffull_name = 'LEA'; } 

		my $fstart = $feature_ref->{start};
		my $fstop  = $feature_ref->{stop};
		my @feature_string;	
		push (@feature_string, $ftype);
		push (@feature_string, $ffull_name);
		push (@feature_string, $fname);
		push (@feature_string, $fstart);
		push (@feature_string, $fstop);
		if ($ftype eq 'UTR' and $fname eq 'LTR') {
			my $feature_string = join ("\t", @feature_string);
			push (@file, "\n$feature_string");
		}
		else {
			#print "\n\t No 5' LTR in $refseq_name";
			sleep 1;
		}

	}


	{
		# Do features section
		my $features_ref   = $self->{features};	
		my $feature_ref = shift @$features_ref;	
		my $ftype      = $feature_ref->{type};
		my $fname      = $feature_ref->{name};
		my $ffull_name = $feature_ref->{full_name};
		my $fstart = $feature_ref->{start};
		my $fstop  = $feature_ref->{stop};
		my @feature_string;	
		push (@feature_string, $ftype);
		push (@feature_string, $ffull_name);
		push (@feature_string, $fname);
		push (@feature_string, $fstart);
		push (@feature_string, $fstop);
		if ($ftype eq 'UTR' and $fname eq 'LEA') {
			my $feature_string = join ("\t", @feature_string);
			push (@file, "\n$feature_string");
		}
		else {
			#print "\n\t No LEA in $refseq_name";
			sleep 1;
		}

	}
	
	# Do features section
	my $features_ref   = $self->{genes};	
	foreach my $feature_ref (@$features_ref) {
		#print "<BR> WRITING $name";
		my @feature_string;
		my $type      = $feature_ref->{type};
		my $name      = $feature_ref->{name};
		my $full_name = $feature_ref->{full_name};
		push (@feature_string, $type);
		push (@feature_string, $full_name);
		push (@feature_string, $name);
		my $starts    = $feature_ref->{starts};
		my $exons_ref = $feature_ref->{exons};
		foreach my $start (@$starts) {
			my $stop = $exons_ref->{$start};
			push (@feature_string, $start);
			push (@feature_string, $stop);
		}
		my $coordinate_string = join("\t", @feature_string);
		push (@file, "\n$coordinate_string");
	}

	{
		# Do features section
		my $features_ref   = $self->{features};	
		my $feature_ref = shift @$features_ref;	
		my $ftype      = $feature_ref->{type};
		my $fname      = $feature_ref->{name};
		my $ffull_name = $feature_ref->{full_name};
		my $fstart = $feature_ref->{start};
		my $fstop  = $feature_ref->{stop};
		my @feature_string;	
		push (@feature_string, $ftype);
		push (@feature_string, $ffull_name);
		push (@feature_string, $fname);
		push (@feature_string, $fstart);
		push (@feature_string, $fstop);
		if ($ftype eq 'UTR' and $fname eq 'LTR') {
			my $feature_string = join ("\t", @feature_string);
			push (@file, "\n$feature_string");
		}
		else {
			#print "\n\t No 3' LTR in $refseq_name";
			sleep 1;
		}

	}

	push (@file, "\nEndblock;");
	
	
	push (@file, "\n\nORIGIN");
	my $sequence = $self->{sequence};	
	my @sequence = split('', $sequence);
	my $i = 0;
	my $seqline = '';
	foreach my $nt (@sequence) {
		$i++;
		$seqline .= $nt;
		if ($i eq 80) {
			push(@file, "\n$seqline");
			$i = 0;
			$seqline = '';
		}
	}

	#$devtools->print_web_hash($self);
	chomp $seqline;
	push(@file, "\n$seqline");
	push (@file, "\n//");
	push (@file, "\nEnd;");
	my $path = $directory_path . $refseq_name;
	$fileio->write_file($path, \@file);
}

#***************************************************************************
# DEVELOPMENT
# Subroutine:  create_linear_formatted 
# Description: create the representation of a reference sequence as a 
#              linear genome formatted for text printing
# (this is more or less finished but needs to be updated using the exon
# model for reference sequences
#***************************************************************************
sub create_linear_formatted {

	my ($self, $linear_ref) = @_;
	
	# Get base and component objects
	my $sequence  = $self->{sequence};
	my $genes_ref = $self->{genes};
	#$devtools->print_array($genes_ref);

	# Write full sequence part
	my $layer = 0;
	$linear_ref->{$layer} = $sequence;
	
	# Get ORFs in each frame
	my $seq_obj = Sequence->new();
	my %orfs;
	$self->get_translated_orfs(\%orfs);

	my %layers;
	foreach my $feature (@$genes_ref) {
		
		#$devtools->print_hash($feature); die;
		my $domain_type  = $feature->{type};
		my $name         = $feature->{name};
		my $start        = $feature->{start};
		my $stop         = $feature->{stop};
		my $frame        = $seq_obj->get_frame($start);
		my $orf          = $orfs{$name};
		print "\n\n\t $name is in frame $frame: ORF $orf";
		if ($layers{$frame}) {
			my $orfs_ref = $layers{$frame};
			$orfs_ref->{$start} = $orf;
		}
		else {
			my %orfs;
			$orfs{$start} = $orf;
			$layers{$frame} = \%orfs;
		}
	}

	my @frames = qw [ 1 2 3 ];
	foreach my $frame (@frames) {
		
		my $linear = '';
		my $current_len = length $linear;
		my $orfs_ref = $layers{$frame};
		my @starts   = sort by_number keys %$orfs_ref;
		#$devtools->print_array(\@starts);
		
		foreach my $start (@starts) {
			my $gap_len = $start - $current_len;
			my $leader = ' ' x $gap_len;
			$linear .= $leader;
			my $prot = $orfs_ref->{$start};
			
			print "\n\n\t Frame $frame, start $start: $prot\n\n";
			
			# Create amino spacing
			my @prot = split('', $prot);
			my $i = 0;
			foreach my $aa (@prot) {
				$i++;
				if ($i eq 1) { $linear .= "$aa ";  }
				else         { $linear .= " $aa "; }
			}
			$current_len = length $linear;
		}
		
		if ($linear_ref->{$frame}) {
			my $tmp = $linear_ref->{$frame};
			#print "\n\tbig $tmp";
			$tmp .= $linear;
			$linear_ref->{$frame} = $tmp;
			#exit;	
		}
		else {
			$linear_ref->{$frame} = $linear;
		}
	}
	#$devtools->print_hash($linear_ref);
}

#***************************************************************************
# Subroutine:  write_linear_formatted_seq 
# Description: 
#***************************************************************************
sub write_linear_formatted_seq {

	my ($self, $linear_ref, $outfile) = @_;
	
	my @layers   = qw [ 3 2 1 0 ];
	my %formatted;
	my $line_total = 0;
	foreach my $layer (@layers) {

		my $sequence = $linear_ref->{$layer};
		my %page_formatted;
		my $char_number = 0;
		my $line_number = 0;
		my $seq_line    = '';
		my @sequence = split('', $sequence);
		foreach my $char (@sequence) {
			
			$char_number++;
			$seq_line .= $char;
			if ($char_number % 80 eq 0) {
				# Write numbers for seqlines
				$line_number++;
				my $start = ($char_number - 80) + 1;
				my $num_len = length $start;
				my $leader;
				my $f_seq_line;
				if  ($layer eq 0) { 
					my $space_len =  7 - $num_len;
					$leader = $start;
					$leader .= ' ' x $space_len;
					$seq_line = $seq_line . "  $char_number";
				}
				else {
					$leader = '       ';
				}
				
				$f_seq_line .= $leader . $seq_line;
				$page_formatted{$line_number} = $f_seq_line;
				$seq_line = '';
			}
		}
		
		# Get the trailing line
		$line_number++;
		$line_total = $line_number;
		my $start = ($char_number - 80) + 1;
		my $num_len = length $start;
		my $leader;
		my $f_seq_line;
		if  ($layer eq 0) { 
			my $space_len =  7 - $num_len;
			$leader = $start;
			$leader .= ' ' x $space_len;
			$seq_line = $seq_line . "  $char_number";
		}
		else {
			$leader = '       ';
		}
			
		$f_seq_line .= $leader . $seq_line;
		$page_formatted{$line_number} = $f_seq_line;
		$formatted{$layer} = \%page_formatted; 
	}
	
	# Write out page formatted sequnce
	#$devtools->print_hash(\%formatted);	
	
	my @output;
	my $line_number = 0;
	do {
		$line_number++;
		foreach my $layer (@layers) {
			my $layer_ref  = $formatted{$layer};
			my $layer_line = $layer_ref->{$line_number};
			$layer_line .= "\n";
			if ($layer eq 0) {
				$layer_line .= "\n";
			}
			push (@output, "$layer_line");
		}
	} until ($line_number eq $line_total);

	# Write out data
	$fileio->write_file($outfile, \@output);
}

############################################################################
# Writing Functions
############################################################################

#***************************************************************************
# Subroutine:  create_new_refseq
# Description: 
#***************************************************************************
sub create_new_refseq {

	my ($self, $data_ref, $refseqfile_text) = @_;

	#$devtools->print_hash($data_ref); die;
	
	my $feature_type    = $data_ref->{feature_type};
	my $feature_name    = $data_ref->{feature_name};
	my $csv_coordinates = $data_ref->{csv_coordinates};
	my $coordinates_ref = $data_ref->{coordinates};
	my $sequence        = $data_ref->{sequence};
	my $genus           = $data_ref->{genus};
	my $family          = $data_ref->{family};
	my $subgroup        = $data_ref->{subgroup};
	my $seq_type        = $data_ref->{seq_type};
	my $eve_type        = $data_ref->{eve_type};
	my $genome_type     = $data_ref->{genome_type};
	my $accession       = $data_ref->{accession};
	my $refseq_name     = $data_ref->{refseq_name};

	# If an append then add feature to existing
	push (@$refseqfile_text, "Begin Metadata;\n");
	push (@$refseqfile_text, "name=$refseq_name;\n");
	push (@$refseqfile_text, "family=$family;\n");
	push (@$refseqfile_text, "group=$genus;\n");
	if ($accession) {
		push (@$refseqfile_text, "accession=$accession;\n");
	}
	push (@$refseqfile_text, "Endblock;\n\n\n");
	
	
	# Add the features if new
	my @features;
	if ($feature_type and $feature_name and $csv_coordinates) {
		push (@features, "Begin Features;\n");
		
		my @feature_string;
		push (@feature_string, $feature_type);
		push (@feature_string, $feature_name);
		push (@feature_string, $feature_name);
		my @coordinates = split(',', $csv_coordinates);
		my $start = undef;
		my %coordinates;
		foreach my $position (@coordinates) {
			if ($start) { 
				my $stop = $position;
				$coordinates{$stop} = $start;
				$start = undef;
			}
			else { 
				$start = $position;
			}	
		}
		if ($start) { die; } # Not an even number of coordinates
		push (@feature_string, @coordinates);
		my $feature_string = join("\t", @feature_string);
		push (@features, "$feature_string\n");
		push (@features, "Endblock;\n\n\n");
		push (@$refseqfile_text, @features);

	}


	# Add the sequence
	push (@$refseqfile_text, "ORIGIN;\n");
	my $count = 0;
	my @sequence = split ('', $sequence);
	my $seqstring = '';
	foreach my $char (@sequence) {
		$count++;
		$seqstring .= $char;
		my $remainder = $count % 80;
		unless ($remainder) {
			$seqstring .= "\n";
		}
	}
	push (@$refseqfile_text, $seqstring);
	push (@$refseqfile_text, "//\n");

	# To do: Add any new mutation
}

############################################################################
# DESCRIBE 
############################################################################

#***************************************************************************
# Subroutine:  describe 
# Description: describe 
#***************************************************************************
sub describe {

	my ($self) = @_;

	my $seq_obj  = Sequence->new();
	print "\n\t### Refseq description";

	my $sequence  = $self->{sequence};
	my $name      = $self->{name};
	my $seq_len   = $self->{seq_len};
	my $start     = $self->{start};
	my $stop      = $self->{stop};
	my $genes_ref = $self->{genes};

	print "\n\t### Name:     $name";
	print "\n\t### Start:    $start";
	print "\n\t### Stop:     $stop";
	print "\n\t### Length:   $seq_len";
	print "\n\t### Sequence: $sequence";

	print "\n\t### Locations of genes within this refseq:";
	foreach my $gene_ref (@$genes_ref) {

		my $gene_name  = $gene_ref->{name};
		my $gene_start = $gene_ref->{start};
		my $gene_stop  = $gene_ref->{stop};
		my $gene_coding_start = $gene_ref->{coding_start};
		my $gene_coding_stop  = $gene_ref->{coding_stop};
		print "\n\t###\t $gene_name spans $gene_start-$gene_stop"; 
		print "\n\t###\t coding spans $gene_coding_start-$gene_coding_stop"; 
		#$devtools->print_hash($gene_ref);
		my $gene_seq = $seq_obj->extract_subsequence($sequence, $gene_start, $gene_stop);
		

		# Iterate through the exons
		my $exon = 0;
		my $starts = $gene_ref->{starts};
		my $exons  = $gene_ref->{exons};
		my $orf = '';
		foreach my $exon_start (@$starts) {
			$exon++;
			my $exon_stop = $exons->{$exon_start};
			
			print "\n\t###\t\t exon $exon: spans $exon_start-$exon_stop"; 
			if ($exon_start < $gene_coding_start) {
				$exon_start = $gene_coding_start;
			}
			if ($exon_stop > $gene_coding_stop) {
				$exon_stop = $gene_coding_stop;
			}
			print "\n\t###\t\t exon coding $exon_start-$exon_stop"; 
			
			my $exon_seq = $seq_obj->extract_subsequence($sequence, $exon_start, $exon_stop);
			my $exon_orf = $seq_obj->translate($exon_seq);
			$orf .= $exon_orf;
		}
		
		print "\n\t$gene_name\n$orf \n\n"; 
	
	}
	print "\n\t### END\n\n\n";

}

#***************************************************************************
# Subroutine:  by number
# Description: by number - for use with perl 'sort'  (cryptic but works) 
#***************************************************************************
sub by_number { $a <=> $b }	

############################################################################
# EOF 
############################################################################
