#!/usr/bin/perl -w
############################################################################
# Script:       RefSeqAlignment.pm 
# Description:  Functions for working with GLUE alignments 
# History:      Rob Gifford, April 2011: Creation
############################################################################
package RefSeqAlignment;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::IO;
use Base::FileIO;
use Base::SeqIO;
use Base::BioIO;
use Base::DevTools;

############################################################################
# Globals
############################################################################
my $io       = IO->new();
my $fileio   = FileIO->new();
my $seqio    = SeqIO->new();
my $bioio    = BioIO->new();
my $devtools = DevTools->new();
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Create a new GLUE MSA object
# Arguments:   $refseq:   A reference to a RefSeq object 
#              $seqs_ref: A reference to an array of sequences
#            -  -  -     in FASTA OR VGLUE stitchback format  -  -  -  
# optional:    $data_ref: A reference to an array containing a table of values
# Values are linked to sequences by sequence ids (first column in table)
# The first row of the table contains column headings
#***************************************************************************
sub new {

	my ($invocant, $params_ref) = @_;
	
	my $class = ref($invocant) || $invocant;
		
	my %mutation_frequencies;
	my %mutation_counts;
	my %codon_frequencies;
	my %codon_counts;
	my %gene_counts;

	my $self = {
		
		# Member variable set at instantiation
		refseq               => $params_ref->{refseq},
		refseq_name          => $params_ref->{refseq}->{name},
		sequences            => $params_ref->{sequences},
		hash_sequences       => $params_ref->{hash_sequences},
		
		# Linking sequences to data
		sequence_ids         => $params_ref->{sequence_ids},
		header_id_link       => $params_ref->{header_id_link},
	
		# Alignment Statistics
		mutation_frequencies => \%mutation_frequencies,
		mutation_counts      => \%mutation_counts,
		codon_frequencies    => \%codon_frequencies,
		codon_counts         => \%codon_counts,
		gene_counts          => \%gene_counts,
		
		# Coordinates describing coverage relative to reference sequence
		lowest_start         => $params_ref->{lowest_start},
		highest_stop         => $params_ref->{highest_stop},
		msa_length           => $params_ref->{msa_length},
		
		# Data structures for stratifiying datasets
		fields_array         => $params_ref->{fields_array},
		field_tracker        => $params_ref->{field_tracker},

	};

	bless ($self, $class);
	return $self;
}

############################################################################
# Get alignment statistics
############################################################################

#***************************************************************************
# Subroutine:  profile sequences
# Description: determine what genes and mutations are present in each
#              sequence in the alignment, by comparison to the reference
#***************************************************************************
sub profile_sequences {

	my ($self) = @_;

	# Get data from self
	my $seq_ids_ref      = $self->{sequence_ids};
	my $aligned_seqs_ref = $self->{hash_sequences};
	my $refseq           = $self->{refseq};
	my $refgenes         = $refseq->{genes};
	#$devtools->print_hash($aligned_seqs_ref); die;
	
	# Iterate through, profiling each sequence in the alignment
	foreach my $seq_id (@$seq_ids_ref) {
		# Get this sequence and profile against the reference
		my $seq_ref = $aligned_seqs_ref->{$seq_id};
		my $seq_id = $seq_ref->{header};
		my $id     = $seq_ref->{sequence_id};
		#print "<BR><BR> ############## SEQ ID: $seq_id  ($id)"; # DEBUG
		$refseq->compare_to_aligned_sequence($seq_ref);
		#$devtools->print_hash($seq_ref); die;
	}
	$self->{mutation_profiled} = 'true';
}

############################################################################
# Deriving statistics and frequencies
############################################################################

#**************************************************************************
# Subroutine:  calculate_coverage
# Description: calculate the proportion of nucleotides in MSA versus 
#               insertion table
#***************************************************************************
sub calculate_coverage {

	my ($self) = @_;

	my $sequences = $self->{sequences};
	my %coverage;
	my $total_nucs = 0;
	my $total_msa_nucs = 0;
	foreach my $seq_ref (@$sequences){
		
		#$devtools->print_hash($seq_ref); die;
		my $rawseqlength = $seq_ref->{raw_seq_length};
		unless ($rawseqlength) { die; }
		$total_nucs     += $rawseqlength;
		my $aln_seq		 = $seq_ref->{aln_seq};
		$aln_seq =~ s/-//g;
		my $alnlength	 = length($aln_seq);
		$total_msa_nucs += $alnlength;
		my $indels_array = $seq_ref->{indels};
		my $indels		 = $indels_array->[1];
		my $id			 = $seq_ref->{sequence_id};
		#$indels->[1][$x] $x begins in 1
		#my $ind_num		= scalar($indels->[1]);
		my $indel_count= $rawseqlength - $alnlength;
		my @chars =();
		#print "\n$id\t$rawseqlength\t$alnlength\t$indel_count\n";
		my $alncov	= ($alnlength/$rawseqlength);
		my $indelcov= ($indel_count/$rawseqlength);
		my %coverage_data;
		$coverage_data{align_coverage} = $alncov;
		$coverage_data{indel_coverage} = $indelcov;
	}
	#$devtools->print_hash(\%coverage); die;
	my $total_insertion_nucs = $total_nucs - $total_msa_nucs;
	$self->{total_nt} = $total_nucs;
	$self->{msa_nt}   = $total_msa_nucs;
	$self->{insertion_nt} = $total_insertion_nucs;
	$self->{coverage} = \%coverage;
}

#**************************************************************************
# Subroutine:  get_gene_counts
# Description: 
#***************************************************************************
sub get_gene_counts {

	my ($self) = @_;

	# Get data form self
	my $aligned_seqs_ref = $self->{hash_sequences};
	my $seq_ids_ref      = $self->{sequence_ids};
	my $refseq           = $self->{refseq};
	my $refseq_start     = $refseq->{start};
	my $refseq_stop      = $refseq->{stop};
	#print "\n\t start '$refseq_start', stop '$refseq_stop'";

	# Iterate through, profiling each sequence in the alignment
	my %gene_counts;
	my $lowest;
	my $highest;
	foreach my $seq_id (@$seq_ids_ref) {
		
		# Get this sequence
		my $align_seq_ref = $aligned_seqs_ref->{$seq_id};
		my $genes = $align_seq_ref->{gene_data};
		unless ($genes) { print "<br> ID: $seq_id no genes"; }
		#print "<br> ID: $seq_id";

		my $gene_name;
		my $gene_ref;
		while ( ( $gene_name, $gene_ref ) = each %$genes) {
			
			my $minimum_present = $gene_ref->{minimum_present};
			unless ($minimum_present) { next; }
			if ($gene_counts{$gene_name}) {
				$gene_counts{$gene_name}++;
			}
			else {
				$gene_counts{$gene_name} = 1;
			}
		}	
	}
	
	# Store data
	$self->{gene_counts} = \%gene_counts;
}

#**************************************************************************
# Subroutine:  get_codon_counts
# Description: 
#***************************************************************************
sub derive_codon_counts {

	my ($self, $field) = @_;

	# Create reference sequence fasta
	my $refseq           = $self->{refseq};
	my $seq_ids_ref      = $self->{sequence_ids};
	my $aligned_seqs_ref = $self->{hash_sequences};
	my $sequence_ids_ref = $self->{sequence_ids};
	my $header_id_link   = $self->{header_id_link},
	my $codon_alignment  = $self->{codon_alignment};
	unless ($refseq)         { die; }
	unless ($sequence_ids_ref) { die; }
	unless ($field) { $field = 'sequence_id'; }
	
	# Iterate through the sequence IDs
	my $codon_counts_ref = $self->{codon_counts};
	foreach my $seq_id (@$sequence_ids_ref) {
		
		#print "\n\t ### Sequence ID: $seq_id";
		# Get the aligned sequence
		my $header = $header_id_link->{$seq_id};
		my $align_seq_ref  = $aligned_seqs_ref->{$seq_id};
		#my $align_seq_ref  = $aligned_seqs_ref->{$header};
		unless ($align_seq_ref) { 
			$devtools->print_hash($aligned_seqs_ref);
			die; 
		}
		
		# Get codon counts for each gene
		my $refgenes = $refseq->{genes};
		foreach my $gene_ref (@$refgenes) {
		
			my $gene_name      = $gene_ref->{name};
			my $genes_data_ref = $align_seq_ref->{gene_data};	
			my $data_ref       = $align_seq_ref->{data};	
			my $gene_data_ref  = $genes_data_ref->{$gene_name};
			
			#print "\n\t ### Gene name: $gene_name";
			#$devtools->print_hash($align_seq_ref); die;
			#$devtools->print_hash($genes_data_ref); die;
			#$devtools->print_hash($gene_data_ref); die;
			#$devtools->print_hash($data_ref); die;
			
			my $positions_ref = $gene_data_ref->{codons};
			my @positions = keys %$positions_ref;
			foreach my $position_key (@positions) {
				
				# Sequence state jiggery pokery
				#print "\n\t<BR> ### KEY: $position_key";
				my $seq_state;
				if ($field eq 'sequence_id') {
					$seq_state = 'total';
				}
				else {
					$seq_state = $data_ref->{$field};
				}
				unless ($seq_state) { 
					$seq_state = 'NULL';
					#die; 
				}
				
				if ($codon_counts_ref->{$position_key}) {
					my $position_data = $codon_counts_ref->{$position_key};
					if ($position_data->{$field}) {
						my $pos_field_data = $position_data->{$field};
						if ($pos_field_data->{$seq_state}) {
							$pos_field_data->{$seq_state}++;
						}
						else {
							$pos_field_data->{$seq_state} = 1;
						}
					}
					else {
						my %position_field_data;
						$position_field_data{$seq_state} = 1;
						$position_data->{$field} = \%position_field_data;
					}
				}
				else {
					my %position_field_data;
					$position_field_data{$seq_state} = 1;
					my %position_data;
					$position_data{$field} = \%position_field_data;
					$codon_counts_ref->{$position_key} = \%position_data;
				}
			}
		}
	}

	# DEBUG
	#$devtools->print_array($fields_array_ref); die;
	#$devtools->print_hash($codon_counts_ref); die;
}


#***************************************************************************
# Subroutine:  derive_mutation_counts
# Description: count mutations relative to the reference squence 
#***************************************************************************
sub derive_mutation_counts {

	my ($self, $field) = @_;

	# Get data from self
	my $seq_ids_ref      = $self->{sequence_ids};
	my $seqs_ref         = $self->{hash_sequences};
	my $sequence_ids_ref = $self->{sequence_ids};
	my $codon_alignment  = $self->{codon_alignment};
	unless ($field) { $field = 'sequence_id'; }

	# Iterate through, profiling each sequence in the alignment
	my $mutation_counts_ref = $self->{mutation_counts};
	foreach my $seq_id (@$sequence_ids_ref) {
		
		# Get this sequence
		my $seq_ref  = $seqs_ref->{$seq_id};
		my $data_ref = $seq_ref->{data};
		#print "\n\t ### Sequence ID: $seq_id";
		
		my $mutations_ref = $seq_ref->{nonsyn_mutations};
		foreach my $mutation_ref (@$mutations_ref) {
			
			#$devtools->print_hash($mutation_ref);
			my $key = $mutation_ref->create_mutation_key();
			#print "<BR>\n\t SEQ ID $seq_id: key $key";
			#print "\n\t ### Mutation: $key";
			
			# Sequence state jiggery pokery
			my $state;
			if ($field eq 'sequence_id') {
				$state = 'total';
			}
			else {
				$state = $data_ref->{$field};
			}
			unless ($state) {
				#$devtools->print_hash($seq_ref);
				$state = 'NULL';
				#die "\n\t found no state for $field for $seq_id\n\n";
				#die;
			}
			
			if ($mutation_counts_ref->{$key}) {
				my $position_data = $mutation_counts_ref->{$key};
				if ($position_data->{$field}) {
					my $field_data = $position_data->{$field};
					$field_data->{$state}++;
				}
				else {
					my %field_data;
					$field_data{$state} = 1;
					$position_data->{$field} = \%field_data;
				}
			}
			else {	
				my %field_data;
				$field_data{$state} = 1;
				my %position_data;
				$position_data{$field} = \%field_data;
				$mutation_counts_ref->{$key} = \%position_data;
			}	
		}	
	}
	# DEBUG
	#$devtools->print_hash($mutation_counts_ref); die;
}

#***************************************************************************
# Subroutine:  derive_mutation_frequencies
# Description: determine frequencies of mutations (relative to refseq)  
#***************************************************************************
sub derive_mutation_frequencies {

	my ($self, $field) = @_;

	unless ($field) { $field = 'sequence_id'; }
	my $mutation_counts  = $self->{mutation_counts};
	my $codon_counts     = $self->{codon_counts};
	my $field_tracker    = $self->{field_tracker};
	unless ($field_tracker) { die; }
	
	unless ($codon_counts) {
		#$devtools->print_hash($codon_counts); die;
		die "\n\t need codon counts to calculate frequencies\n\n";
	}
	
	# Now calculate a freq for all positions
	my $mutation_freqs_ref = $self->{mutation_frequencies};
	my $mutation_key;
	my $mutation_data;
	while ( ( $mutation_key, $mutation_data ) = each %$mutation_counts) {
		
		# Create position key for gene from mutation key
		my @key = split (":", $mutation_key);
		my $gene     = $key[0];
		my $position = $key[2];
		my $position_key = $gene . ':' . $position;
		#print "\n\t MUTATION $mutation_key";
		
		# Get the number of valid codons (the denominator)
		my $position_data_ref = $codon_counts->{$position_key};
		unless ($position_data_ref) { die; }
		#my $total_data = $position_data_ref->{sequence_id};
		#my $total_aln_codons = $total_data->{total};;
		
		my $data2_ref = $position_data_ref->{$field};
		#$devtools->print_hash($position_data_ref); 
		
		# Get info on the character/field we are using to stratify the data
		my $field_data = $mutation_data->{$field};	
		unless ($field_data) {
			print "\n\t 1 no data found for '$field'\n\n";
			die;
		}
		#$devtools->print_hash($field_tracker); die;
		#$devtools->print_hash($field_data);
		#$devtools->print_hash($mutation_data); die;
		
		# Iterate through each 'state' for this character
		#print "\n\t field '$field'";
		my $states_hash = $field_tracker->{$field};
		unless ($states_hash) { die "\n\t No states for field '$field'"; }
		my @states = keys %$states_hash;
		foreach my $state (@states) {
		
			#print "\n\t field '$field' state '$state'";
			my $total_data = $position_data_ref->{$field};
			my $total_aln_codons = $total_data->{$state};
			my $field_data  = $mutation_data->{$field};
			my $state_count = $field_data->{$state};
			unless ($state_count) {
				#$devtools->print_hash($field_data); die;
				#die;
				$state_count = '0';
			}

			if ($field ne 'sequence_id') {
				my $field_counts = $mutation_data->{$field};
				unless ($field_counts) { die; }
				#$total_aln_codons = $field_counts->{$state};
			}
			
			# Work out frequency
			my $f_freq;
			if ($total_aln_codons) {
				my $freq = $state_count / $total_aln_codons;
				$f_freq = sprintf("%.3f", $freq);
				$f_freq .= " ($state_count/$total_aln_codons)";
				#print "\n\t frequency $f_freq";
			}
			else {
				$total_aln_codons = '0';
				$f_freq  = '0';
				$f_freq .= " ($state_count/$total_aln_codons)";
				# may indicate a problem need CHECK
				#$devtools->print_hash($position_data_ref);
				#$devtools->print_hash($states_hash);
				#$devtools->print_hash($total_data);
				#print "\n\t Position $position_key got no total for field";
				#print "\n '$field' state '$state'\n\n";
			} 

			if ($mutation_freqs_ref->{$mutation_key}) {
				my $position_data = $mutation_freqs_ref->{$mutation_key};
				if ($position_data->{$field}) {
					my $field_data = $position_data->{$field};
					$field_data->{$state} = $f_freq;
				}
				else {
					my %field_data;
					$field_data{$state} = $f_freq;
					$position_data->{$field} = \%field_data;
				}
			}
			
			else {
				my %field_data;
				$field_data{$state} = $f_freq;
				my %position_data;
				$position_data{$field} = \%field_data;
				$mutation_freqs_ref->{$mutation_key} = \%position_data;
			}
		}
	}
	#$devtools->print_hash($mutation_freqs_ref); die;
}

#***************************************************************************
# Subroutine:  create_mutation_frequency_table
# Description: create a table of mutation frequencies, ordered by gene and position 
#***************************************************************************
sub create_mutation_frequency_table {

	my ($self, $stratify_field) = @_;
	
	#print "\n\t Writing stratified mutation frequencies, please wait...\n";
	my $field_tracker       = $self->{field_tracker};
	my $mutation_freqs      = $self->{mutation_frequencies};
	my $refseq              = $self->{refseq};
	my $refseq_genes        = $refseq->{genes};

	# Create header row
	my @frequency_data;
	my @row;
	my @fields = qw [ sequence_id ];
	if ($stratify_field) {
		push (@fields, $stratify_field);
	}
	foreach my $field (@fields) {
		
		if ($field eq 'sequence_id') {
			#print "\n\t TOTAL";
			unshift (@row, "total");
			next;
		}
		else {
			#print "\n\t FIELD $field";
			my $states_ref = $field_tracker->{$field};
			my @states = keys %$states_ref;
			foreach my $state (@states) {
				my $column_name = "$field ($state)";
				#print "\n\t NAME ($column_name)";
				push (@row, $column_name);
			}
		}
	}
	unshift (@row, 'aa');
	unshift (@row, 'position');
	unshift (@row, 'ref_aa');
	unshift (@row, 'gene');
	my $row = join ("\t", @row);
	push (@frequency_data, "$row\n");

	# Get all mutations in gene order
	my @mutations;
	my @mutation_keys = keys %$mutation_freqs;
	my %gene_ordered_data;
	foreach my $key (@mutation_keys) {
		
		#print "\n<BR> KEY $key";
		my @key = split(':', $key);
		my $gene_name  = $key[0];
		my $position   = $key[2];
		my $aa         = $key[3];
	
		if ($gene_ordered_data{$gene_name}) {
		
			my $position_data = $gene_ordered_data{$gene_name};
			if ($position_data->{$position}) {
				my $changes_ref = $position_data->{$position};
				push (@$changes_ref, $aa);
			}
			else {
				my @changes;
				push (@changes, $aa);
				$position_data->{$position} = \@changes;
			}
		}
		else {
			my @changes;
			push (@changes, $aa);
			my %position_data;
			$position_data{$position} = \@changes;
		 	$gene_ordered_data{$gene_name} = \%position_data;
		}
	}

	my @sorted_keys;
	foreach my $gene (@$refseq_genes) {
		my $indexed_aas = $gene->{indexed_aas};
		my $gene_name   = $gene->{name};
		my $gene_positions_ref = $gene_ordered_data{$gene_name};
		my @positions = sort by_number keys %$gene_positions_ref;
		foreach my $position (@positions) {
			
			my $ref_aa = $indexed_aas->{$position};
			unless ($ref_aa) { 
				print "\n\t <BR>No ref AA at position $position in $gene_name";
				$ref_aa = 'X';
			}

			my $changes_ref = $gene_positions_ref->{$position};
			my @changes = sort @$changes_ref;
			foreach my $aa (@changes) {

				my @key;
				push (@key, $gene_name);
				push (@key, $ref_aa);
				push (@key, $position);
				push (@key, $aa);
				my $key = join(':', @key);
				#print "\n\t key = '$key'";
				push (@sorted_keys, $key);
			}
		}
	}
	#$devtools->print_array(\@sorted_keys); 

	# Create table data
	foreach my $key (@sorted_keys) {
	
		my @key = split (":", $key);
		my $gene     = $key[0];
		my $ref_aa   = $key[1];
		my $position = $key[2];
		my $aa       = $key[3];
		
		my @freqs;
		my $freq_data = $mutation_freqs->{$key};
		foreach my $field (@fields) {
			my $field_data = $freq_data->{$field};	
			my $states_ref = $field_tracker->{$field};
			my @states_array;
			if ($field eq 'sequence_id') {
				push (@states_array, 'total');	
			}
			else {
				@states_array = keys %$states_ref;
			}
			foreach my $state (@states_array) {
				#print "\n\t ### '$field', '$state'";	
				#$devtools->print_hash($field_data);
				my $freq_string = $field_data->{$state};		
				unless ($freq_string) { $freq_string = '-'; }
				if ($field eq 'sequence_id') {
					unshift (@freqs, $freq_string);	
				} 
				else {
					push (@freqs, $freq_string);	
				}
			}
		}
		#my $mutation_string = "$gene\t$ref_aa\t$position\t$aa";
		my $mutation_string = "$gene\t$ref_aa\t$position\t$aa";
		my $freq_string = join("\t", @freqs);
		my $data_row = $mutation_string . "\t" . $freq_string;
		push (@frequency_data, "$data_row\n");
	}
	
	# Store mutation frequencies and gene frequencies
	$self->{frequency_table} = \@frequency_data;
}

#***************************************************************************
# Subroutine:  derive_typical_list 
# Description: 
#***************************************************************************
sub derive_typical_list {
	
	my ($self, $typical_list_ref, $threshold) = @_;

	# Get stratifiers
	my $table_ref = $self->{frequency_table};
	unless ($table_ref and $threshold) { die; }	

	# Create table data rows
	my $head_row_string = shift @$table_ref;
	foreach my $row_string (@$table_ref) {
		
		chomp $row_string;
		my @row_data = split("\t", $row_string);
		my $gene     = shift @row_data;;
		my $ref_aa   = shift @row_data;
		my $position = shift @row_data;
		my $aa       = shift @row_data;
		my $mut_key  = $gene  . ':' . $ref_aa;
		   $mut_key .= ':' . $position . ':' . $aa;
		my $i = 1;
		my $typical = undef;
		my $freq_1 = $row_data[0];
		my @freq = split(/\s+/, $freq_1);
		my $freq = shift @freq;	
		if ($freq >= $threshold) {
			#print "\n\t adding to $mut_key to list ";
			$typical_list_ref->{$mut_key} = $freq; 	
		}
	}
	#$devtools->print_hash($typical_list_ref);
}

############################################################################
# Alignment manipulations
############################################################################

#***************************************************************************
# Subroutine:  do_cpg_correction 
# Description:  
#***************************************************************************
sub do_cpg_correction {

	my ($self, $threshold) = @_;

	unless ($threshold) {
		$threshold = 0.1;
	}
	$self->get_pattern_frequencies('dinucleotides', 1);
	$self->correct_cpg($threshold);
}

#***************************************************************************
# Subroutine:  correct_cpg 
# Description:  
#***************************************************************************
sub correct_cpg {

	my ($self, $threshold) = @_;
	
	# get data from self
	my $freqs_ref  = $self->{dinucleotides}->{frequencies};
	my $counts_ref = $self->{dinucleotides}->{counts};
	my $seqs_ref   = $self->{sequences};	
	
	# Set up the range 
	my $msa_length = $self->{msa_length};	
	my $count_stop = $msa_length - 1;	
	my @msa_positions = 1..$count_stop;
	
	# Make decisions on each position
	foreach my $position (@msa_positions) {
		
		my $position_freqs  = $freqs_ref->{$position};
		my $position_counts = $counts_ref->{$position};
		my $cpg_freq  = $position_freqs->{CG};
		my $cpg_count = $position_counts->{CG};
		print "\n\t Position $position: $cpg_freq ($cpg_count)";
		if ($cpg_count) {
			if ($cpg_freq > $threshold) {
				$devtools->print_hash($position_counts);
				$devtools->print_hash($position_freqs);
			}		
		}
	}
}

#***************************************************************************
# Subroutine:  get_pattern_frequencies 
# Description: TODO (not sure - find out)
#***************************************************************************
sub get_pattern_frequencies {

	my ($self, $pattern_name, $extend) = @_;
	
	# get data from self
	my $seqs_ref   = $self->{sequences};	
	my $seq_obj = @$seqs_ref[0];
	unless ($seq_obj) { die; }
	my $dinucs_ref = $seq_obj->{dinucs};	
	
	# Set up the range 
	my $msa_length = $self->{msa_length};	
	my $count_stop = $msa_length - 1;	
	my @msa_positions = 1..$count_stop;
	
	# Get freq and prop values for alignment dinucs
	my %freqs;
	my %counts;
	foreach my $position (@msa_positions) {
		
		my %freq;
		my %count;
		my $valid;
		print "\n\t # DOING position $position";
		foreach my $seq_ref (@$seqs_ref) {

			#$devtools->print_hash($seq_ref);
			#die;
			my $sequence = $seq_ref->{sequence};
			my $id       = $seq_ref->{sequence_id};
			
			#print "\n\t # ID: $id position $position";
			my @nucs = split ('', $sequence);
			my $i = $position - 1;
			my $pattern = $nucs[$i];
			unless ($pattern) { die; }
			if ($extend) {
				my $extended = 0;
				do {
					$i++;
					$extended++;
					my $char = $nucs[$i];
					unless ($char) { $char = '-'; }
					$pattern .= $char;
				} until ($extended eq $extend);
			}
			if ($pattern =~ '-' ) { next; }
			$valid++;
			# Increment count
			#print "\n\t$pattern";	
			if ($count{$pattern}) { $count{$pattern}++;   }
			else                  { $count{$pattern} = 1; }
		}
		#$devtools->print_hash(\%count);
		$counts{$position} = \%count;
		
		# convert counts to proportions/frequencies
		my @patterns = keys %counts;
		foreach my $pattern (@patterns) {
			my $count = $count{$pattern}; 
			my $freq;
			if ($valid and $count) {
				$freq = $count / $valid;
			}
			else {
				$freq = 0;
			}
			my $f_freq = sprintf("%.3f", $freq);
			$freq{$pattern} = $f_freq;
		}
		$freqs{$position} = \%freq;
		
		#$devtools->print_hash(\%freq);
	}
	
	# DEBUG
	#$devtools->print_hash(\%counts);
	#$devtools->print_hash(\%freqs);
	my %results;
	$results{frequencies} = \%freqs;
	$results{counts}      = \%counts;
	$results{$pattern_name} = \%results;
	$devtools->print_hash(\%results); die;
}

#***************************************************************************
# Subroutine:  change pattern
# Description:  
#***************************************************************************
sub change_pattern {

	my ($self, $search, $replace) = @_;

	# get data from self
	my $seqs_ref   = $self->{sequences};	
	
	# Set up the range 
	my $msa_length = $self->{msa_length};	
	my $count_stop = $msa_length - 1;	
	my @msa_positions = 1..$count_stop;
	
	# Get freq and prop values for alignment dinucs
	foreach my $position (@msa_positions) {
		print "\n\t # DOING position $position";
		foreach my $seq_ref (@$seqs_ref) {
			$seq_ref->do_pattern_conversion($position, $search, $replace);
		}
	}	
}

############################################################################
# Shannon entropy
############################################################################

#***************************************************************************
# Subroutine:  derive_shannon_entropy
# Description: evaluates shannon entropy of aligned sequences (DBM) 
#***************************************************************************
sub derive_shannon_entropy {
	
	my ($self) = @_;
	my $sequences = $self->{sequences};	#get the sequences array
	my $refseq = $self->{refseq};
	my %results;
	my @shannon_e;	
	my $number_aln=0; 
	my @line;
	my @genes;
	my $gene;
	my @characters;
	my @align;
	my $c=0;
	my $i;
	my $char;
	my $nchar;
	my $header;
	my $sequence;
	my %test_hash_obj;
	my $test_ref = \%test_hash_obj;;
	#### GET EACH CHARACTER IN ALIGNED SEQUENCES INTO AN ARRAY OF ARRAYS
	foreach my $aln (@$sequences){	#For each of the aln sequences
		$number_aln++;
		#$devtools->print_hash($aln); next;
		$header = $aln->{header}; 	#get the header
		$sequence = $aln->{aln_seq}; 	#get the aligned sequence
		my %translated_orfs;
		$refseq->{sequence}=$sequence;
		$refseq->get_translated_orfs(\%translated_orfs);
		#$devtools->print_hash(\%translated_orfs);
		#next;
		@genes = keys(%translated_orfs);
		my %aln_genes;
		foreach $gene (@genes){
			@characters = split(//, $translated_orfs{$gene});	#get each character separatelly
			$nchar = scalar(@characters);
			@{$aln_genes{$gene}} = @characters;
			#@{$align[$number_aln]} = @characters;	#put each array of characters in an array of arrays
			@line=();
			@characters=();
		}
		%{$align[$number_aln]} = %aln_genes;
	}
	#die;
	my $total;
	my %Count;
	my ($H, $p, $tc, $j);
	###FOR EACH GENE IN THE ALIGNMENT
	foreach $gene (@genes){
		$nchar=scalar(@{$align[1]{$gene}});
		###GO THROUGH THE ARRAY OF ARRAYS BY COLUMN
		for($i=0; $i < $nchar; $i++){
    		$total=0;
	    	#Initialize hash %Count
			foreach $char (keys %Count){
				delete $Count{$char};
			}   
			#####GO THROUGH THE ARRAY OF ARRAYS BY ROW
			for($j=1; $j<=$number_aln; $j++){
				$char = $align[$j]{$gene}[$i]; 	#get the character
				#print "\n$j\t$gene\t$nchar\t$i\t$char\n"; die;
				if($char =~ /(\-|\?|N)/){	#if the charcter is not defined or a gap get the next character
					next;
				}else{	
					$total++;	#counts the number of allowed characters
				}   
				if($Count{$char}){	#Count the appearance of each character
					$Count{$char}++;	
				}else{
					$Count{$char}=1;
				}   
			}   
			$H =0; 	
			####Calculates Shannon Entropy H = -SUM( PROBABILITY_OF_CHARACTER * LOG2(PROBABILITY_OF_CHARACTER))
			foreach $char (keys %Count){	
				$p =$Count{$char}/$total; 	#probability of each character
				$H +=$p *log($p); #sum of the entropy
			}   
			$H = -$H/log(2);	#transform to log base 2
			$tc= $i + 1;	#position in the alignment
			my $f_H = sprintf("%.2f", $H); 
			my $line = "$tc\t$f_H\t$total";  #POSITION SHANNON_ENTROPY TOTAL_NUMBER_OF_ACCEPTED_CHARACTERS_IN_COLUMN
			push (@shannon_e, $line);	#add each site H to the results array
		}
		@{$results{$gene}} = @shannon_e;
		@shannon_e = ();
	}
	$self->{shannon_results} = \%results;
	#$devtools->print_hash(\%results); 
}



############################################################################
# DESCRIBE & WRITE FXNs
############################################################################

#***************************************************************************
# Subroutine:  write_glue_msa
# Description: 
#***************************************************************************
sub write_glue_msa {

	my ($self, $path) = @_;
	
	# Get settings/paths/data from self	
	my $refseq       = $self->{refseq};
	my $refseq_name  = $self->{refseq}->{name};
	my $start        = $self->{lowest_start};
	my $stop         = $self->{highest_stop};
	my $aln_seqs     = $self->{sequences};

	# Create the GLUE MSA header
	unless ($refseq_name and $start and $stop) { die; }
	my $header = "#GLUE $refseq_name: $start-$stop\n\n";
	my @glue;
	push (@glue, $header);

	# Create the sequence data secion
	foreach my $seq_ref (@$aln_seqs) {
		
		# Get the sequence data
		#$devtools->print_hash($seq_ref);
		my $header  = $seq_ref->{header};
		#my $aln_seq = $seq_ref->{aln_seq};
		my $aln_seq = $seq_ref->{sequence};
	
		# Create FASTA
		my $fasta = "\n>$header\n$aln_seq";
		push(@glue, $fasta);
	}
	
	# Write the alignment (DEBUG)
	$fileio->write_file($path, \@glue);
	#die;	
}

#***************************************************************************
# Subroutine:  describe
# Description: show text output describing this GLUE MSA object
#***************************************************************************
sub describe {

	my ($self, $output_type) = @_;
	
	my $message = "# Alignment Summary";
	$io->show_output_message($message, $output_type);
	my $refseq       = $self->{refseq};
	my $refseq_name  = $refseq->{name};
	my $refseq_start = $refseq->{start};
	my $refseq_stop  = $refseq->{stop};
	my $seqs_ref     = $self->{hash_sequences};
	my $seq_ids_ref  = $self->{sequence_ids};
	
	$message = "# Sequences are profile aligned to $refseq_name ($refseq_start-$refseq_stop)";
	$io->show_output_message($message, $output_type);
	
	my $num_seqs = scalar @$seq_ids_ref;
	$message = "# Alignment contains $num_seqs sequences";
	$io->show_output_message($message, $output_type);

	# Store data as member variables		
	my $fields_array_ref = $self->{fields_array};

	# Iterate through, profiling each sequence in the alignment
	foreach my $seq_id (@$seq_ids_ref) {
		
		my $seq_ref = $seqs_ref->{$seq_id};
		#$devtools->print_hash($seq_ref);exit;

		my $sequence = $seq_ref->{sequence};
		#my $sequence = $seq_ref->{aln_seq};
		my $seq_len  = length $sequence;
		$message = "# Sequence $seq_id: length '$seq_len'";
		$io->show_output_message($message, $output_type);
		
		#$devtools->print_hash($seq_ref);
		my $inserts_ref = $seq_ref->{insertions};
		foreach my $insert (@$inserts_ref) {

			$message = "# Insertion $insert";
			$io->show_output_message($message, $output_type);
		}
		
		#$devtools->print_hash($seq_ref);
		foreach my $field (@$fields_array_ref) {
			my $value = $seq_ref->{$field};
			$message = "# $field: $value";
			$io->show_output_message($message, $output_type);
		}
		
		my $mutations_ref = $seq_ref->{nonsyn_mutations};
		foreach my $mutation_ref (@$mutations_ref) {
			my $mutation_key = $mutation_ref->create_mutation_key();
			$message = "# Mutation: $mutation_key";
			$io->show_output_message($message, $output_type);

		}
	}
	
	$message = "# Alignment Summary Ends";
	$io->show_output_message($message, $output_type);
}

#***************************************************************************
# Subroutine:  by number
# Description: by number - for use with perl 'sort'  (cryptic but works) 
#***************************************************************************
sub by_number { $a <=> $b }	

############################################################################
# EOF
############################################################################
