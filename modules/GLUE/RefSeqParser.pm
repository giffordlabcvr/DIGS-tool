#!/usr/bin/perl -w
############################################################################
# Script:       RefSeqParser 
# Description:  Class for parsing GLUE reference sequence files 
# History:      Rob Gifford, April 2010: Creation
#               RG updated May 2011 - support for spliced genes
############################################################################
package RefSeqParser;

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

############################################################################
# Globals
############################################################################

# Base objects
my $bioio    = BioIO->new();
my $fileio   = FileIO->new();
my $devtools = DevTools->new();
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create a new reference sequence parser object 
#***************************************************************************
sub new {

	my ($invocant, $file_path) = @_;
	my $class = ref($invocant) || $invocant;

	# Member variables
	my $self = {
	
	};

	bless ($self, $class);
	return $self;
}

############################################################################
# Parsing Functions
############################################################################

#***************************************************************************
# Subroutine:  parse_refseq_flatfile
# Description: upload references sequences from flat files
# Notes:  rvdb flatfiles for consensus sequences use a relaxed NEXUS format
#***************************************************************************
sub parse_refseq_flatfile {

	my ($self, $file_path, $refseq_data_ref, $start, $stop) = @_;

	# Read in the file
	my @refseq_file;
	my $status = $fileio->read_file($file_path, \@refseq_file);
	unless ($status and \@refseq_file) { 
		return 0; 
	}
	#$devtools->print_array(\@refseq_file); die;
	
	#=== Get the entire nucleic acid sequence	
	# Use the genbank fxn (the same start and end tokens are used: 'ORIGIN' & '//')
	my $sequence = $self->get_sequence(\@refseq_file);
	$sequence =~ s/\s+//g;    # IMPORTANT!!!!!
	$sequence = uc $sequence;
	$refseq_data_ref->{sequence} = $sequence;

	#=== Parse the metadata
	my %metadata;
	$self->parse_refseq_metadata(\@refseq_file, \%metadata);
	$refseq_data_ref->{metadata} = \%metadata;

	#=== Parse the alignment data block
	$self->parse_alignment_block(\@refseq_file, $refseq_data_ref);

	#=== Extract and parse the features block
	$self->parse_features(\@refseq_file, $refseq_data_ref);

	#=== Mutation lists
	$self->parse_mutlists_block(\@refseq_file, $refseq_data_ref);
	#$devtools->print_hash($refseq_data_ref);

	#=== Mutation frequencies
	#$self->parse_mutfreqs_block(\@refseq_file, $refseq_data_ref);
	# NOTE COPY TO FXN BELOW WHEN DONE

	#=== Sort out ORFs for alignment
	$self->set_defaults($refseq_data_ref, $start, $stop);
	if ($start and $stop) {
		
		# If we're dealing with a subsequence, extract it
		my $seq_obj = Sequence->new();
		my $sequence  = $refseq_data_ref->{sequence};
		my $subseq    = $seq_obj->extract_subsequence($sequence, $start, $stop);
		unless ($subseq) { die; }

		# Store the subsequence and its data
		$refseq_data_ref->{sequence} = $subseq;
		$refseq_data_ref->{start}    = $start;
		$refseq_data_ref->{stop}     = $stop;
		$refseq_data_ref->{seq_len}  = $stop - $start + 1;
	}

	#=== Adjust coordinates for the reference sequence based on data
	$self->set_coordinates($refseq_data_ref);

	# Create indexed sequences for main sequence and genes
	$self->create_indexed_sequences($refseq_data_ref, $start);
	#$devtools->print_hash($refseq_data_ref);

	return 1;
}

#***************************************************************************
# Subroutine:  parse_refseq_datastructure
# Description: same as above, but operates on a datastructure not a file
#***************************************************************************
sub parse_refseq_datastructure {

	my ($self, $refseq_filedata_ref, $refseq_data_ref, $start, $stop) = @_;

	#=== Get the entire nucleic acid sequence	
	# Use the genbank fxn (the same start and end tokens are used: 'ORIGIN' & '//')
	my $sequence = $self->get_sequence($refseq_filedata_ref);
	$sequence =~ s/\s+//g;    # IMPORTANT!!!!!!!!!
	$sequence = uc $sequence;
	$refseq_data_ref->{sequence} = $sequence;

	#=== Parse the metadata
	my %metadata;
	$self->parse_refseq_metadata($refseq_filedata_ref, \%metadata);
	$refseq_data_ref->{metadata} = \%metadata;
	#$devtools->print_hash($refseq_data_ref);
	
	#=== Parse the alignment data block
	$self->parse_alignment_block($refseq_filedata_ref, $refseq_data_ref);

	#=== Extract and parse the features block
	$self->parse_features($refseq_filedata_ref, $refseq_data_ref);
	#$devtools->print_array($refseq_filedata_ref);
	#$devtools->print_hash($refseq_data_ref);

	#=== Sort out ORFs for alignment
	$self->set_defaults($refseq_data_ref, $start, $stop);
	if ($start and $stop) {
		
		# If we're dealing with a subsequence, extract it
		my $seq_obj = Sequence->new();
		my $sequence  = $refseq_data_ref->{sequence};
		my $subseq    = $seq_obj->extract_subsequence($sequence, $start, $stop);
		unless ($subseq) { die; }
		
		# Store the subsequence and its data
		$refseq_data_ref->{sequence} = $subseq;
		$refseq_data_ref->{start}    = $start;
		$refseq_data_ref->{stop}     = $stop;
		$refseq_data_ref->{seq_len}  = $stop - $start + 1;
	}

	#=== Adjust coordinates for the reference sequence based on data
	$self->set_coordinates($refseq_data_ref);

	# Create indexed sequences for main sequence and genes
	$self->create_indexed_sequences($refseq_data_ref, $start);
	#$devtools->print_hash($refseq_data_ref);
}

############################################################################
# Load mutation data 
############################################################################

#***************************************************************************
# Subroutine:  parse_mutlists_block 
# Description: retrieve data from the metadata block of a refseq file
#***************************************************************************
sub parse_mutlists_block {

	my ($self, $refseq_file_ref, $refseq_data) = @_;
	
	# Extract the metadata block
	#$devtools->print_array($refseq_file_ref);
	my $start_mark = 'MutList;';
	my $end_mark   = 'Endblock';
	my @mut_lists;
	$fileio->extract_text_blocks($refseq_file_ref, \@mut_lists, 
	                                          $start_mark, $end_mark);
	#$devtools->print_array(\@mut_lists);
	#$devtools->print_hash($refseq_data); die;

	$self->load_mutation_lists(\@mut_lists, $refseq_data);
}

#***************************************************************************
# Subroutine:  load_mutation_list
# Description: load mutation list data as an array of hashes 
#***************************************************************************
sub load_mutation_lists {

	my ($self, $lists_ref, $data_ref) = @_;

	# Load list data
	my $got_typical = undef;
	my %hash_lists;
	my @listnames;
	my %list_colors;
	foreach my $list_ref (@$lists_ref) {
		
		#if ($list eq 'typical') { $got_typical = 'true'; }
		my %list;
		my @array_list;
		my %hash_list;
		# Iterate through list reading and inserting each row of data
		my $name;
		foreach my $line (@$list_ref) {
			
			if ($line =~ /^\s*$/) { next; } # discard blank line
			chomp $line;
			if ($line =~ /=/) { 
			
				my @bits = split('=', $line);
				my $field = $bits[0];
				$field =~ s/\s+//;
				my $value = $bits[1];
				$value =~ s/\s+//;
				$value =~ s/;//;
				$list{$field} = $value;
				if ($field eq 'name') {
					$name = $value;
				}
			}
			else {
				my @data = split("\t", $line);
				#$devtools->print_array(\@data);
				my %data;
				$data{listname}      = $data[0];
				$data{version}       = $data[1];
				$data{gene_name}     = $data[2];
				$data{reference_aa}  = $data[3];
				$data{position}      = $data[4];
				$data{aa}            = $data[5];
				$data{phenotype}     = $data[6];
				push (@array_list, \%data);

				my $gene_name     = $data[2];
				my $reference_aa  = $data[3];
				my $position      = $data[4];
				my $aa            = $data[5];
					
				unless ($gene_name and $reference_aa and $position and $aa) {
					die "$line";
				}
				
				my $key = $gene_name . ':' . $reference_aa . ':' . $position . ':' . $aa;
				$hash_list{$key} = \%data;
			}
		}
		
		# Store list data
		$list{hash_list}  = \%hash_list;
		$list{array_list} = \@array_list;
		$hash_lists{$name} = \%list;
		push (@listnames, $name);
	}

	# Store the list data
	$data_ref->{hash_lists}  = \%hash_lists;
	$data_ref->{listnames}   = \@listnames;
	if ($got_typical) {
		$data_ref->{typical} = 'defined';
	}

	# Check if typical polymorphisms have been defined
	#unless ($got_typical) {
		# If not derive from mutation prevalences file if present
		#if ($component_subdir{mutation_prevalences}) {
		#	my @prevalences;
		#	my $prevalences_path = $path . '/mutation_prevalences';		
		#	$fileio->read_file($prevalences_path, \@prevalences);
		#	#$self->derive_typical_list_from_prevalence_file($reference_name, \@prevalences)
		#}
	#}

}

#***************************************************************************
# Subroutine:  parse_mutfreqs_block 
# Description: retrieve data from the metadata block of a refseq file
#***************************************************************************
sub parse_mutfreqs_block {

	my ($self, $refseq_file_ref, $mutfreqs_ref) = @_;
	
	# Extract the metadata block
	my $start_mark = 'MutFreqs';
	my $end_mark   = 'Endblock;';
	$fileio->read_standard_field_value_block($refseq_file_ref, $start_mark, 
	                                          $end_mark, $mutfreqs_ref);
	#$devtools->print_hash($metadata_ref);
		#my $mutfreqs_ref = $refseq_data_ref->{mutfreqs};
		#$self->load_mutation_frequencies($mutfreqs_ref, $refseq_data_ref);
	die;
	
}

#***************************************************************************
# Subroutine:  parse_refseq_metadata
# Description: retrieve data from the metadata block of a refseq file
#***************************************************************************
sub parse_refseq_metadata {

	my ($self, $refseq_file_ref, $metadata_ref) = @_;
	
	# Extract the metadata block
	my $start_mark = 'Metadata';
	my $end_mark   = 'Endblock;';
	$fileio->read_standard_field_value_block($refseq_file_ref, $start_mark, 
	                                          $end_mark, $metadata_ref);
	#$devtools->print_hash($metadata_ref);
}

#***************************************************************************
# Subroutine:  parse_alignment_block 
# Description: parse the alignment details block of a RefSeq file
#***************************************************************************
sub parse_alignment_block {

	my ($self, $refseq_file_ref, $refseq_data_ref) = @_;

	# Extract the alignment block
	my @alignment_block; 
	$fileio->extract_text_block($refseq_file_ref, \@alignment_block, 
	                                          'Alignment', 'Endblock;');
	#$devtools->print_array(\@alignment_block);
	my %exclude_regions;
	foreach my $line (@alignment_block) {
		
		chomp $line;
		my @line = split("\t", $line);
		my $type = $line[0];
		if ($type eq 'Exclude') {
			my %exclude;
			$exclude{gene}  = $line[1];
			$exclude{start} = $line[2];
			$exclude{stop}  = $line[3];
			if ($line[4]) {
				$exclude{name}  = $line[4];
			}
			my $gene     = $exclude{gene};
			my $position = $exclude{start};
			my $stop     = $exclude{stop};
			do {
				my $position_key = $gene . ':' . $position;
				$exclude_regions{$position_key} = \%exclude;
				$position++;
			} until ($position eq $stop);
		}
	}
	$refseq_data_ref->{exclude} = \%exclude_regions;
	#$devtools->print_array($refseq_data_ref->{exclude});

}

#***************************************************************************
# Subroutine:  parse_features 
# Description: extract sequence feature details from the refseq file
#***************************************************************************
sub parse_features {

	my ($self, $refseq_file_ref, $refseq_data_ref) = @_;

	my $seq_len = $refseq_data_ref->{seq_len};
	my @features_block;
	$fileio->extract_text_block($refseq_file_ref, \@features_block, 'Features', 'Endblock;');
		
	# group by type
	my %by_type;
	my %by_name;
	my @names;
	foreach my $line (@features_block) {
	
		if  ($line =~ /^\s*$/)   { next; } # discard blank line
		if  ($line =~ /^#/)      { next; } # discard comment line
		#print $line;
		chomp $line;
		my @data = split("\t", $line);
		my $type      = shift @data;
		my $full_name = shift @data;
		my $name      = shift @data;
		push (@names, $name);
		#$devtools->print_array(\@data);

		# Iterate through the coordinates for this feature
		my @starts;
		my %start_indexed_exons;	
		my $i = 0;
		my $last_start = undef;
		my $exons = 0;
		foreach my $coordinate (@data) {
			$i++;
			my $remainder = $i % 2;	
			if ($remainder) {
				# this is an odd number so a start
				push (@starts, $coordinate);
				$last_start = $coordinate;
			}
			else {
				# this is an even number so a stop
				unless ($last_start) { 
					die "Error processing Refseq in features BLOCK, gene '$name' \n\n"; 
				}
				$start_indexed_exons{$last_start} = $coordinate;
				$exons++;
			}
		}
		
		my $first_start = undef;
		my $stop = undef;
		foreach my $start (@starts) {
			unless ($first_start) {
				$first_start = $start;
			}
			$stop = $start_indexed_exons{$start};
		}
		my $start = $first_start;
		my %data;
		$data{type}      = $type;
		$data{full_name} = $full_name;
		$data{name}      = $name;
		$data{start}     = $start;
		$data{stop}      = $stop;
		$data{starts}    = \@starts;
		$data{exons}     = \%start_indexed_exons;
		$data{num_exons} = $exons;
				
		#$devtools->print_hash(\%data);
		if ($by_type{$type}) {
			my $array_ref = $by_type{$type};
			push (@$array_ref, \%data); 
		}
		else {
			my @array;
			push (@array, \%data);
			$by_type{$type} = \@array;
		}
	}
	# DEBUG
	#$devtools->print_hash(\%by_type); die;	
	#my @sorted_orfs;
	#my @starts = sort by_number keys %by_start;
	#foreach my $start (@starts) {
	#	my $orf_ref = $by_start{$start};
	#	push (@sorted_orfs, $orf_ref);
	#}	

	#$refseq_data_ref->{genes}    = \@sorted_orfs;
	my $orf_ref = $by_type{ORF};
	#$devtools->print_array($orf_ref); # die;
	$refseq_data_ref->{genes}    = $orf_ref;
	#$refseq_data_ref->{genes}    = $by_type{ORF};
	$refseq_data_ref->{features} = $by_type{UTR};

	# make the genes hash
	my $genes_ref = $by_type{ORF};
	my %genes_hash;
	foreach my $gene_ref (@$genes_ref) {
		my $name = $gene_ref->{name};
		$genes_hash{$name} = $gene_ref;
	}
	$refseq_data_ref->{hash_genes} = \%genes_hash;

	# DEBUG
	#$devtools->print_hash($refseq_data_ref); die;	
	#$devtools->print_array(\@features_block); die;	
	#$devtools->print_hash(\%by_type); die;
	
}

#***************************************************************************
# Subroutine:  set_defaults 
# Description: Set default position coordinates 
#***************************************************************************
sub set_defaults {

	my ($self, $refseq_data_ref, $start, $stop) = @_;

	# Get the complete sequence and calculate its length
	my $sequence = $refseq_data_ref->{sequence};
	my $seq_len = length $sequence;
	
	# Set defaults (full length) if no start and stop params received
	if ($start and $stop) {
		$seq_len = $stop - $start + 1;
	}
	unless ($start)  {  $start = 1;       }
	unless ($stop)   {  $stop = $seq_len; }

	$refseq_data_ref->{start}   = $start;
	$refseq_data_ref->{stop}    = $stop;
	$refseq_data_ref->{seq_len} = $seq_len;
}

#***************************************************************************
# Subroutine:  set_coordinates 
# Description: 
#***************************************************************************
sub set_coordinates {
	
	my ($self, $refseq_data_ref) = @_;

	#$devtools->print_hash($refseq_data_ref);
	my $refseq_start = $refseq_data_ref->{start};
	my $refseq_stop  = $refseq_data_ref->{stop};
	my $frag_length  = $refseq_stop - $refseq_start + 1;
	my $genes_ref    = $refseq_data_ref->{genes};
	
	# Iterate through genes
	my @reset;
	foreach my $gene_ref (@$genes_ref) {
			
		# Get data from gene reference
		my $name   = $gene_ref->{name};
		my $start  = $gene_ref->{start};
		my $stop   = $gene_ref->{stop};
		unless ($start and $stop) {
			#$devtools->print_hash($refseq_data_ref);
			#$devtools->print_array($genes_ref);
			print "$name\n";
			die;
		}

		# Skip features not in range
		unless ($stop >= $refseq_start and $start <= $refseq_stop) {
		   next;
		}
		$self->reset_gene($refseq_data_ref, $gene_ref, $refseq_start, $frag_length);	

		push (@reset, $gene_ref);
	}
	
	#$devtools->print_array(\@reset);
	$refseq_data_ref->{genes} = \@reset;
}

#***************************************************************************
# Subroutine:  reset_gene 
# Description: 
#***************************************************************************
sub reset_gene {

	my ($self, $refseq_data_ref, $gene_ref, $refseq_start, $frag_length) = @_;

	#$devtools->print_hash($gene_ref);
	my $name   = $gene_ref->{name};
	my $start  = $gene_ref->{start};
	my $stop   = $gene_ref->{stop};
	my $starts = $gene_ref->{starts};
	my $exons  = $gene_ref->{exons};
	
	# Create variables for adjusted ranges
	my $adjust_start = $start - $refseq_start;
	my $adjust_stop  = $stop  - $refseq_start + 1;
	
	# DEBUG CHECK
	#print "\n\t HERE: $name: $adjust_start = $start - $refseq_start + 1";

	my @new_starts;
	my %new_exons;
	foreach my $start (@$starts) {
		
		my $stop = $exons->{$start};
		unless ($stop) { die; }
		my $exon_new_start = $start - $refseq_start + 1;
		my $exon_new_stop  = $stop  - $refseq_start  + 1;
		push (@new_starts, $exon_new_start);
		$new_exons{$exon_new_start} = $exon_new_stop;
	}
	$gene_ref->{starts} = \@new_starts; 
	$gene_ref->{exons}  = \%new_exons; 

	# Reset coordinates
	$gene_ref->{start} = $adjust_start; 
	$gene_ref->{stop}  = $adjust_stop;

	# If this ORF starts upstream of our refseq then 
	# set the coding start nt position (first nt first viable codon)
	my $coding_start = $adjust_start;
	if ($adjust_start < 1) {
		my $positive = $adjust_start * -1;
		my $remainder = $positive % 3;
		#print "\n\t # Remainder after dividing $positive by 3 = $remainder";
		$coding_start = 1 + ($remainder);
		#print "\n\t # $name coding start: 1 + (3 - $remainder) = $coding_start ";
	}
	$gene_ref->{coding_start} = $coding_start; 
	
	# If this ORF ends downstream of our refseq then 
	# set the coding end position (last nt of last viable codon)
	my $coding_stop = $adjust_stop;
	if ($adjust_stop >= $frag_length) {
		my $remainder = $frag_length % 3;
		$coding_stop =  $frag_length - $remainder ;
		#print "\n\t $name coding stop $coding_stop";
	}
	$gene_ref->{coding_stop} = $coding_stop; 
	
	#$devtools->print_hash($gene_ref);
	#my $refseq_seq = $refseq_data_ref->{sequence};
	#my $seq_obj = Sequence->new();
	#my $sub_sequence = $seq_obj->extract_subsequence($refseq_seq, $coding_start, $coding_stop);
	#$gene_ref->{sequence} = $sub_sequence; 

}

#***************************************************************************
# Subroutine:  create_indexed_sequences
# Description: 
#***************************************************************************
sub create_indexed_sequences {

	my ($self, $refseq_data_ref, $start) = @_;

	# Index the main sequence
	my %indexed_seq;
	my $seq_obj = Sequence->new();
	my $sequence = $refseq_data_ref->{sequence};
	$seq_obj->index_sequence($sequence, \%indexed_seq, $start);
	$refseq_data_ref->{indexed_seq} = \%indexed_seq;
	
	# Initialise hashes for indexed sequences
	my $genes_ref = $refseq_data_ref->{genes};
	foreach my $gene_ref (@$genes_ref) {

		my $gene_name  = $gene_ref->{name};
		my $gene_start = $gene_ref->{start};
		my $gene_stop  = $gene_ref->{stop};
		my $gene_coding_start = $gene_ref->{coding_start};
		my $gene_coding_stop  = $gene_ref->{coding_stop};
		#print "\n\t###\t $gene_name spans $gene_start-$gene_stop"; 
		#print "\n\t###\t coding spans $gene_coding_start-$gene_coding_stop"; 
		#$devtools->print_hash($gene_ref);

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
			unless ($exon_orf) {
				print "\n\t # $exon_start-$exon_stop";
			}
			else {
				$orf .= $exon_orf;
			}
		}
		
		my %aas;
		my %codons;
		my @aas = split ('', $orf);
		my $i=0;
		foreach my $aa (@aas) {
			$i++;
			$aas{$i} = $aa;
			$gene_coding_start++;
		}
		$gene_ref->{indexed_aas} = \%aas;
	}
}

#***************************************************************************
# Subroutine:  get_sequence
# Description: 
#***************************************************************************
sub get_sequence {

	my ($self, $file_ref, $f_width) = @_;

	# Set default width for sequence if width not set
	unless ($f_width) { $f_width = 80; }

	# Extract seq matrix
	my @sequence;
	$fileio->extract_text_block($file_ref, \@sequence, 'ORIGIN', '//');
	#$devtools->print_array(\@sequence);
	my $sequence = join("\n", @sequence);
	$sequence =~ s/\d//g; # remove numbers
    $sequence =~ s/[\s\/]//g;
	$sequence = uc $sequence; # We like sequences to be upper case

	# Format width
	my @f_sequence = split('', $sequence); 
	my $i = 0;
	my $f_sequence;
	foreach my $char (@f_sequence) {
		$i++;
		$f_sequence .= $char;
		my $modulo = $i % $f_width;
		if ($modulo eq 0) { # Insert a line break
			$f_sequence .= "\n";
		} 
			
	}
	return $f_sequence;
}

#***************************************************************************
# Subroutine:  load_frequencies
# Description:  
#***************************************************************************
sub load_frequencies {

	my ($self, $frequencies_path, $data_ref) = @_;

	# Read the file om
	my @frequencies_file;
	$fileio->read_file($frequencies_path, \@frequencies_file);
	
	# Get stratification data from the header line
	my $header_line  = shift(@frequencies_file);
	my @header_line  = split("\t", $header_line);
	my $gene_name    = shift @header_line;
	my $reference_aa = shift @header_line;
	my $position     = shift @header_line;
	my $aa           = shift @header_line;
	my $i = 0;
	my %stratifications;
	my @stratifications;
	foreach my $column (@header_line) {
		chomp $column;
		$stratifications{$i} = $column;
		push (@stratifications, $column);
		$i++;
	}

	# Iterate through the frequencies inserting data to DB
	my @all_mutations;
	my %frequency_data;
	foreach my $data_line (@frequencies_file) {
	
		chomp $data_line;
		my @data = split(/\s+/, $data_line);
		#$devtools->print_array(\@data);
		my $gene         = shift @data;
		my $reference_aa = shift @data;
		my $position     = shift @data;
		my $aa           = shift @data;
		my $key = $gene . ':' . $reference_aa . ':' . $position . ':' . $aa; 
		push (@all_mutations, $key);

		my $j = 0;
		foreach my $frequency (@data) {
		
			my $stratification = $stratifications{$j};
			if ($frequency_data{$stratification}) {
				my $hash_ref = $frequency_data{$stratification};
				$hash_ref->{$key}  = $frequency;
			}
			else {
				my %stratification_prevs;
				$stratification_prevs{$key} = $frequency;
				$frequency_data{$stratification} = \%stratification_prevs;
			}
			$j++;
		}
	}
	$data_ref->{frequencies}      = \%frequency_data;
	$data_ref->{frequency_muts}   = \@all_mutations;
	$data_ref->{frequency_strats} = \@stratifications;
}

############################################################################
# BASIC FXNS
############################################################################

sub by_number { $a <=> $b }	


############################################################################
# EOF 
############################################################################
