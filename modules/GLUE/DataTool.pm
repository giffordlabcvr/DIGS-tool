#!/usr/bin/perl -w
############################################################################
# Module:      DataTool.pm
# Description: GLUE-associated FASTA tool
# History:     January 2012: Created by Robert Gifford 
############################################################################
package DataTool;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::BioIO;
use Base::FileIO;
use Base::SeqIO;
use Base::DevTools;
use Base::Console;
use Base::Sequence; # For performing basic sequence manipulations
#use Base::HTML_Utilities;

############################################################################
# Globals
############################################################################

# Create base objects
my $fileio     = FileIO->new();
my $seqio      = SeqIO->new();
my $bioio      = BioIO->new();
my $devtools   = DevTools->new();
my $console    = Console->new();

my $header_limit  = 50;
my $default_delimiter =  "\t";  # default for delimited files
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Parameters
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {
		
		# Member variables
		#process_id           => $parameter_ref->{process_id},
		output_type          => $parameter_ref->{output_type},
		
		# Member classes
		blast_obj 			 =>	$parameter_ref->{blast_obj},
		
		# Paths
		#refseq_use_path      => $parameter_ref->{refseq_use_path}, 
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# SECTION: Data reformatting functions 
############################################################################

#***************************************************************************
# Subroutine:  fasta_to_delimited
# Description: FASTA to delimited 
#***************************************************************************
sub fasta_to_delimited {

	my ($self, $seqs_ref, $data_ref, $delimiter ) = @_;

	# Set the delimiter to default if unset
	unless ($delimiter)  { $delimiter = $default_delimiter; }

	# Reformat the sequences	
	foreach my $seq_ref (@$seqs_ref) {
		my $sequence_id  = $seq_ref->{sequence_id};
		my $header       = $seq_ref->{header};
		my $sequence     = $seq_ref->{sequence};
		my $seq_len      = length $sequence;
		$header =~ s/\s+//g; # Remove whitespace
		$sequence =~ s/\s+//g; # Remove whitespace
		chomp $header;
		chomp $sequence;
		push (@$data_ref, "$header\t$sequence\n");
		#push (@$data_ref, "$sequence_id\t$sequence\n");
	}	
	#my $header = "Sequence ID\tSequence\n";
	#unshift (@$data_ref, $header);	
}

#***************************************************************************
# Subroutine:  fasta_to_phylip
# Description: convert FASTA to PHYLIP 
#***************************************************************************
sub fasta_to_phylip {

	my ($self, $fasta_ref, $phylip_ref, $table_ref, $translation_ref) = @_;
	
	my $set_id_len = 20;

	my $phycount;
	my $last_seq_len;
	foreach my $sequence_ref (@$fasta_ref) {
		
		$phycount++;
		my $header      = $sequence_ref->{header};
		my $sequence_id = $sequence_ref->{sequence_id};
		my $sequence    = $sequence_ref->{sequence};
		chomp $sequence_id;
	
		my $seq_len = length $sequence;
		unless ($seq_len) {
			die;
		}
		unless ($last_seq_len) {
			$last_seq_len = $seq_len;
		}
		elsif ($seq_len eq $last_seq_len) {
			unless ($seq_len eq length $sequence) {
				print "\n\t sequence $sequence_id is a different length ($seq_len)\n
				       \t check FASTA sequences are aligned\n\n";
			}
		}
	
		# Create PHYLIP id
		#my $phy_id = $phycount . '_'; 
		my $phy_id = ''; 
		my @sequence_id = split ('', $sequence_id);
		my $id_len;
		foreach my $char (@sequence_id) {
			$phy_id .= $char;	
			$id_len = length($phy_id);
			if ($id_len eq 10) { last; }
		}
		my $spacer_len = ($set_id_len - $id_len);
		my $spacer = ' ' x $spacer_len;
		my $phy_seq = $phy_id . $spacer . $sequence . "\n";
		push(@$phylip_ref, $phy_seq);
		$translation_ref->{$phy_id} = $header;

		# store id relationship in translation table 
		my $id_pair = "$sequence_id\t$phy_id\n";
		push (@$table_ref, $id_pair);
	}
	
	# Create PHYLIP taxa and characters header
	my $num_taxa  = $phycount;
	my $num_chars = $last_seq_len;
	my $header_line = $num_taxa . '   ' . $num_chars . "\n";
	unshift(@$phylip_ref, $header_line);

}

#***************************************************************************
# Subroutine:  fasta_to_nexus 
# Description: Convert aligned FASTA to NEXUS format 
#***************************************************************************
sub fasta_to_nexus {

	my ($self, $seqs_ref, $data_ref) = @_;

	# Get the maximum sequence and taxa label length (number of chars), & number of taxa
	my $max_length       = 0; 
	my $max_label_length = 0;
	my %clean_sequences;
	my $num_taxa = 0;
	foreach my $seq_ref (@$seqs_ref) {
		
		$num_taxa++;

		# remove illegal chars from tax label
		my $tax_label = $seq_ref->{header};
		my $sequence  = $seq_ref->{sequence};
		unless ($tax_label and $sequence) { die; }
		
		$tax_label =~ s/-//g;
		#$tax_label =~ s/_/x/g;
		$tax_label =~ s/\t/_/g;
		
		my $label_length = length $tax_label;
		my $seq_length = length $sequence;
		if ($seq_length > $max_length) { 
			$max_length = $seq_length; 
		}
		if ($label_length > $max_label_length) { 
			$max_label_length = $label_length;
		}
		$clean_sequences{$tax_label} = $sequence;
	}

	# Put in the NEXUS hader and initial part of Data block
	my $num_char = $max_length; 
	push (@$data_ref, "#NEXUS\n");
	push (@$data_ref, "Begin DATA;\n");
	push (@$data_ref, "Dimensions ntax=$num_taxa nchar=$num_char;\n");
	push (@$data_ref, "Format Datatype=Nucleotide Gap=-;\n");
	push (@$data_ref, "Matrix\n");

	# Do the data matrix
	my $tax_label;
	my $sequence;
	while ( ( $tax_label, $sequence ) = each %clean_sequences) {
		chomp $sequence;
		my $seq_len = length $sequence;
		if ($seq_len < $max_length) {
		
			my $gap_len      = $max_length - $seq_len;
			my $trailing_gap = '-' x $gap_len;
			$sequence .= $trailing_gap;
		}
		
		my $label_length = length $tax_label;
		
		my $pad_length = ($max_label_length + 10) - $label_length;
		my $pad = ' ' x $pad_length;
		my $data_line = $tax_label . $pad . $sequence . "\n";
		push (@$data_ref, $data_line);
	}
	
	push (@$data_ref, "\n;\n");
	push (@$data_ref, "End;\n\n");
}

#***************************************************************************
# Subroutine:  delimited_to_fasta
# Description: convert a delimited text file to FASTA
# Assumptions: the last column in the text file contains the sequence 
#***************************************************************************
sub delimited_to_fasta {
	
	my ($self, $data_ref, $seqs_ref, $delimiter) = @_;
	
	# Set the delimiter to default if unset
	unless ($delimiter)  { $delimiter = $default_delimiter; }

	foreach my $line (@$data_ref) {
		#print "<BR> $line";
		chomp $line;
		my @fields = split(/$delimiter/, $line);
		my $sequence = pop @fields;
		my $header = join("\t", @fields);
		my $fasta = '>' . $header . "\n" . $sequence . "\n";
		push (@$seqs_ref, $fasta);	
	}
}

#***************************************************************************
# Subroutine:  genbank_to_fasta_and_data
# Description: convert GenBank file into a FASTA file and linked data file 
#***************************************************************************
sub genbank_to_fasta_and_data {
	
	my ($self, $infile, $seq_ref, $data_ref) = @_;

	# Declare variables
	my $i = 0;
	my $seq_id;
	my $seq_date;
	my $GBFILE = $infile;
	my $sequence       = '';
	my $accession      = 'NULL';
	my $iso_date       = 'NULL';
	my $iso_country    = 'NK';
	my $genotype       = undef;
	my $serotype       = undef;
	my $isolate        = undef;
	my $host           = undef;
	my $seq_flag       = undef;
	my $feature_flag   = undef;
	my $cds_flag       = undef;
	my $cd_coordinates = undef;
	my $matpep_flag    = undef;
	my $prev_line      = undef;
	my $matpep_coordinates = undef;
	unless (open(GBFILE, "$infile")) {
		print "\n\t Cannot open file \"$infile\"\n\n";
		exit;
	}

	# Iterate through file
	foreach my $line (<GBFILE>) {

	    chomp($line);
		$i++;
		#print "\n\t LINE $i\t $line";
	    # do line-by-line processing.
		if ($line =~ /^\/\//) {  # End of Genbank entry
			
			# Do some processing to get start and stop
			#print "\n\t ## Processed Genbank entry for '$seq_id' ($iso_country, $iso_date)";
			$seq_flag = undef;
			my %data;
			$data{sequence_id}         = $seq_id;
			$data{isolation_country}   = $iso_country;
			$data{isolation_date}      = $iso_date;
			$data{sequence_date}       = $seq_date;
			unless ($sequence and $seq_id) { 
				print "\n\t No sequence for '$seq_id"; sleep 2; exit;
			}
			$data{sequence}            = uc $sequence;
			my $genotype_method  = "genbank_annotation";
			unless ($genotype) { 
				$genotype = 'NULL';
				$genotype_method = 'NA';
			}
			$data{genotype} = $genotype;
			$data{genotype_method} = $genotype_method;
			unless ($serotype) {  $serotype = 'NULL'; }
			$data{serotype} = $serotype;
			unless ($isolate) {  $isolate = 'NULL'; }
			$data{isolate} = $isolate;
			unless ($host) {  $host = 'NULL'; }
			$data{host} = $host;
			#$devtools->print_hash(\%data); die;
			push(@$data_ref, \%data);
			
			$sequence = uc $sequence;
			my $header = $seq_id;
			#print "\n\t sequence ID = '$seq_id'";
			my $seq_obj = Sequence->new($sequence, $header, $seq_id);
			push(@$seq_ref, $seq_obj);

			# RESET
			$sequence    = '';
			$seq_id      = '';
			$seq_date    = '';
			$iso_country = 'NK';
			$genotype    = undef;
			$serotype    = undef;
			$host        = undef;
			$isolate     = undef;

		}
		elsif ($seq_flag) {
			$line =~ s/\s//g;
			$line =~ s/\d//g;
			$sequence .= $line; 
			#print "\n\t Adding to sequence '$seq_id'";
		}
		if ($line =~ /^FEATURES/) {
			#print "\n\tStart of features block";
			$feature_flag = 1;
		}
		if ($feature_flag) {
			if ($line =~ /country/) {
				my @line = split(/"/, $line);
				#$devtools->print_array(\@line);
				$iso_country = $line[1];
				$iso_country =~ s/'/_/g;
				my @iso_country = split (/:/, $iso_country);
				$iso_country = $iso_country[0];
				unless ($iso_country =~ /\w/) { $iso_country = 'NK'; }
				#print "\n\t # GOT country '$iso_country' for $seq_id";
			
			}
			if ($line =~ /collection_date/) {
				my @line = split(/"/, $line);
				$iso_date = $line[1];
				$iso_date =~ s/'/_/g;
				#print "\n\t # GOT date '$iso_date' for $seq_id";
			}
			if ($line =~ /note="genotype/) {
				my @line = split(/"/, $line);
				#$devtools->print_array(\@line);
				$genotype = $line[1];
				my @genotype = split (/;/, $genotype);
				$genotype = $genotype[0];
				$genotype =~ s/genotype://g;
				$genotype =~ s/ //g;
				$genotype =~ s/genotype//g;
				#print "\n\t # GOT genotype '$genotype' for $seq_id";
			}
			if ($line =~ /\/serotype/) {
				my @line = split(/=/, $line);
				#$devtools->print_array(\@line);
				$serotype = $line[1];
				my @serotype = split (/;/, $serotype);
				$serotype = $serotype[0];
				$serotype =~ s/"//g;
				$serotype =~ s/serotype//g;
				print "\n\t # GOT serotype '$serotype' for $seq_id";
			}
			if ($line =~ /\/isolate/) {
				my @line = split(/=/, $line);
				#$devtools->print_array(\@line);
				$isolate = $line[1];
				$isolate =~ s/"//g;
				print "\n\t # GOT isolate '$isolate' for $seq_id";
			}
			if ($line =~ /\/host/) {
				my @line = split(/=/, $line);
				$host = $line[1];
				$host =~ s/"//g;
				print "\n\t # GOT host '$host' for $seq_id";
			}
		}
		if ($line =~ /^ORIGIN/) {
			#print "\n\tStart of sequence for '$seq_id'";
			$seq_flag = 1;
			$feature_flag = undef;
		}
		elsif ($line =~ /^LOCUS/) {
			my @line = split(/\s+/, $line);
			#$devtools->print_array(\@line);
			$seq_id = $line[1];
			$seq_date= pop @line;
			print "\n\t ####### Processing entry ($seq_id) date '$seq_date'";
		}
		#elsif($line =~ /^ACCESSION/) {
		#	my @line = split (/\s+/, $line);
		#	$accession = pop @line;
		#}
		#elsif($compress_line =~ /^ORGANISM/) {
		#	$line =~ s/^\s*ORGANISM\s*//;
		#	$data_ref->{organism} = $line;
		#}
		#elsif(line =~ /^DEFINITION/) {
		#	$line =~ s/^\s*DEFINITION\s*//;
		#	$data_ref->{accession} = $line;
		#}
		#elsif($compress_line =~ /^VERSION/) {
		#	$line =~ s/^\s*VERSION\s*//;
		#	$data_ref->{version} = $line;
		#}
		my $prev_line = $line;
	}
}

############################################################################
# SECTION: Sequence sorting functions
############################################################################

#***************************************************************************
# Subroutine:  sort_seqs_by_length 
# Description: sort a fasta file of sequences in order of descending length
# Arguments: $file (a path)
#            $fasta_ref (array to store output)
#            $count_gaps (flag to count gap characters)
#***************************************************************************
sub sort_seqs_by_length {

	my ($self, $fasta_ref, $sorted_ref, $minimum, $count_gaps) = @_;
	
	# Index the sequences by size
	my %seq_lengths;
	my %sequences;
	my $exclude;
	foreach my $seq_ref (@$fasta_ref) {
		
		my $id           = $seq_ref->{sequence_id};
		my $sequence     = $seq_ref->{sequence};
		$sequences{$id}  = $seq_ref;
		my $adjusted_seq = $sequence;
		my $length;
		if ($count_gaps) {
			# Count any gaps if flag is set to do so
			$length = length $sequence;
		} 
		else {
			# Don't count gaps
			$adjusted_seq =~ s/-//g;
			$length = length $adjusted_seq;
		}
		my $store = undef;
		if ($minimum) {	
			#print "<BR> $minimum, $length";
			if ($length >= $minimum) { 
				$store = 1; 
			}
			else {
				$exclude++;
			}	
		}
		else { $store = 1; }
		if ($store) {
			if ($seq_lengths{$length}) {
				my $seqs_ref = $seq_lengths{$length};
				push (@$seqs_ref, $id);
			}
			else {
				my @seqs;
				push (@seqs, $id);
				$seq_lengths{$length} = \@seqs;
			}
		}
	}

	#my @fasta;
	#unless($self->{mode} eq 'html') {
	#print "\n\tSEQ LENGTHS IN DESCENDING NUMERIC ORDER:\n";
	my @lengths = sort by_number keys (%seq_lengths);
	foreach my $length (@lengths) {
		
		my $seqs_ref = $seq_lengths{$length};
		my $num_seq = scalar @$seqs_ref;
		#print "\n\tLength $length': $num_seq sequences";
		foreach my $id (@$seqs_ref) {
			my $seq_ref  = $sequences{$id};
			my $sequence = $seq_ref->{sequence};
			my $header   = $seq_ref->{header};
			my $fasta    = ">$header\n$sequence\n"; 
			unshift (@$sorted_ref, $fasta);
			#unshift (@$sorted_ref, "$sequence\n");
			#unshift (@$sorted_ref, ">$header\n");
		}
	}
	return $exclude;
}

#***************************************************************************
# Subroutine:  sort_sequences_on_data_column
# Description: Sort sequences using tab data file field
#***************************************************************************
sub sort_sequences_on_data_column {

	my ($self, $seqs_ref, $data_ref, $id_index, $sub_index, $delimiter) = @_;

	# sort the sequences
	my %output;
	foreach my $line (@$data_ref) {

		# get values
		my @line = split("$delimiter", $line);
		my $id         = $line[$id_index];
		my $sort_value = $line[$sub_index];
		$sort_value =~ s/\s//g;
		if ($sort_value eq '') {
			$sort_value = 'NULL';
		}
		chomp $id;         # always chomp for safety
		chomp $sort_value; # always chomp for safety
	
		my $list_ref = $output{$sort_value};
		unless ($list_ref) {
			my @seq_list;
			push (@seq_list, $id);
			$output{$sort_value} = \@seq_list;
		}
		else { 
			push (@$list_ref, $id); 
			$output{$sort_value} = $list_ref;
		}
	}

	my $sort_value;
	my $list_ref;
	while ( ( $sort_value, $list_ref ) = each %output) {
	
		my @fasta;
		foreach my $id (@$list_ref) {

			# get sequence
			my $sequence = $seqs_ref->{$id};
			if ($sequence) {
				my $fasta = ">$id\n$sequence\n";
				push(@fasta, $fasta);
			}
			else {
				print "\n\t Couldn't find sequence $id";
			}	
		}
		my $outfile = $sort_value . '.fas';	
		$outfile =~ s/\s//g;
		$fileio->write_file($outfile, \@fasta);
	}
}

############################################################################
# SECTION: Sequence filtering functions
############################################################################

#***************************************************************************
# Subroutine:  filter_seqs_by_header
# Description: filter a fasta file of sequences 
#***************************************************************************
sub filter_by_header {

	my ($self, $fasta_ref, $sorted_ref, $word, $rule, $ignore_case) = @_;

	my $removed = 0;	
	foreach my $seq_ref (@$fasta_ref) {

		my $header   = $seq_ref->{header};
		my $sequence = $seq_ref->{sequence};
		unless ($header) { die; }
		my $header_copy = $header;
		if ($ignore_case) {
			$word   = uc $word;
			$header_copy = uc $header;
		}

		if ($header_copy =~ /$word/)  {
			if ($rule eq 'include')  {
				push (@$sorted_ref, ">$header\n");
				push (@$sorted_ref, "$sequence\n");
			}
			else {
				$removed++;
			}
		}
		else {
			push (@$sorted_ref, ">$header\n");
			push (@$sorted_ref, "$sequence\n");
		}
	}
	return $removed;
}

#***************************************************************************
# Subroutine:  truncate_long_headers 
# Description: truncate long headers in a FASTA file 
#***************************************************************************
sub truncate_long_headers {

	my ($self, $file) = @_;
	
	my $count;
	my @fasta;
	$seqio->read_fasta($file, \@fasta);
	#$devtools->print_array(\@fasta); die;
	my @converted;
	my @table;
	foreach my $sequence_ref (@fasta) {
		
		$count++;
		print "\n\t #### DOING $count";
		my $sequence_id = $sequence_ref->{header};
		my $sequence    = $sequence_ref->{sequence};
		
		
		my @sequence_id = split ('', $sequence_id);
		my $truncate_id = $sequence_id; 
		my @bits = split (/\|/, $truncate_id);
		$truncate_id = pop @bits;
		
		my @bits2 = split (/\[/, $truncate_id);
		my $end  = shift @bits2;
		my $name = pop @bits2;
		$name =~ s/\]//g;
		$name =~ s/,//g;
		$truncate_id =~ s/\s+/_/g;
		$name =~ s/\s+/_/g;
		$end  =~ s/\s+/_/g;
		$end  =~ s/_\n//g;
		my $end_len = length $end;
		if ($end_len > 20) {
			my @end = split('', $end);
			my $e;
			my $i;
			foreach my $char (@end) {
				$i++;
				$e .= $char;
				if ($i >= 40) { last;}
			}
			$end = $e;
		}

		if ($name) {
			my $new_id = $name . $end . $count;
			$new_id =~ s/,//g;
			$new_id  =~ s/\//_/g;
			$new_id  =~ s/\'/-/g;
			#$truncate_id =~ s/\'//g;
			#$truncate_id =~ s/\"//g;
			#$truncate_id =~ s/\//g;
			print "\n\t OLD:  $sequence_id";
			print "\n\t BITS: $name + $end";
			print "\n\t NOW:  $new_id";
			#my $truncate_id = $phycount . '_'; 
			my $fasta = ">$new_id\n$sequence\n\n";
			push (@converted, $fasta);
			my $id_pair = "$new_id\t$sequence_id";
			push (@table, $id_pair);
		}
	}

	my $fasta_out = $file . '.converted.fas';
	my $table_out = $file . '.table.txt';
	$fileio->write_file($fasta_out, \@converted);
	$fileio->write_file($table_out, \@table);
}

#***************************************************************************
# Subroutine:  split_glue_refseq_file
# Description:  
#***************************************************************************
sub split_glue_refseq_file  {

	my ($self, $file, $refseqs_ref) = @_;

	my $initialised = undef;
	
	# Read in the file
	my @file;
	my $status = $fileio->read_file($file, \@file);
	unless ($status) { 
		die;
		return $status; 
	}

	#=== Get the entire nucleic acid sequence	
	# Use the genbank fxn (the same start and end tokens are used: 'ORIGIN' & '//')

	# iterate through the file
	my $i = 0;
	my $joined = join('', @file);
	#print "NT $joined\n\n\n\n";
	my @split = split("End;\n", $joined);
	foreach my $refseq (@split) {
		$i++;
		#print "\n\t ############### ERE $i \n\n$refseq\n\n"; #die;
		my @refseq = split("\n", $refseq);
		my @ffs;
		foreach my $line (@refseq) {
			push(@ffs, "$line\n");
		}
		push(@$refseqs_ref,  \@ffs); 
	}
	return 1;
}

############################################################################
# TAB-DELIMITED FILE MANIPULATION UTILITIES
############################################################################

#***************************************************************************
# Subroutine:  combine_data
# Description: Use shared IDs to combine lines of data in two files
#***************************************************************************
sub combine_data {
	
	my ($self, $file_a, $file_b, $a_index, $b_index) = @_;
	
	# create lookup hashes
	my %seen_in_a;
	foreach my $line (@$file_a) { 
		chomp $line;
		my @line_bits = split("\t", $line);
		my $id = $line_bits[$a_index];
		$seen_in_a{$id} = $line; 
	}
	
	my %seen_in_b;
	foreach my $line (@$file_b) { 
		chomp $line;
		my @line_bits = split("\t", $line);
		my $id = $line_bits[$b_index];
		$seen_in_b{$id} = $line; 
	}

	#$devtools->print_hash(\%seen_in_b);

	my @in_both;	
	my @in_first_only;
	foreach my $line (@$file_a) { 
		
		chomp $line;
		my @line_bits = split("\t", $line);
		my $id = $line_bits[$a_index];
		
		# Output for debugging etc
		print "\n\t Checking ID $id in file_a\t";
		
		if ($seen_in_b{$id}) { 
			my $combined_line = $line;
			$combined_line .= "\t$seen_in_b{$id}\n";
			#print $combined_line;
			push (@in_both, $combined_line); 
		}
		else { 
		
			# Output for debugging etc
			print "....couldn't find ID $id in file_b";
			# ----- HACK for if you want the line from 1st as well if no match
			push (@in_first_only, "$line\n");
		}
	}

	my $file1 = 'in_both_files.txt';
	$fileio->write_file($file1, \@in_both);

	my $file2 = 'in_first_only.txt';
	$fileio->write_file($file2, \@in_first_only);

}

#**************************************************************************
# Subroutine:  read_glue_msa 
# Description: read in a GLUE multiple sequence alignment file
#***************************************************************************
sub read_glue_msa {
	
	my ($self, $msa_path, $array_ref, $hash_ref) = @_;
	
	# Get the all-important first line
	my @file;
	$fileio->read_file($msa_path, \@file);
	my $header = shift @file;
	chomp $header;
	unless ($header =~ /^\s*#GLUE/) { die; } # is this a GLUE file?
	my @header = split(/\s+/, $header);
	#$devtools->print_array(\@header);
	my $refseq_name  = $header[1];
	$refseq_name =~ s/\s+//g;
	$refseq_name =~ s/://g;
	my $coordinates  = $header[2];
	my @coordinates  = split('-', $coordinates);
	my $start = $coordinates[0];
	my $stop  = $coordinates[1];
	$hash_ref->{refseq_name}  = $refseq_name;
	$hash_ref->{lowest_start} = $start;
	$hash_ref->{highest_stop} = $stop;

	#my @header = split("\s", $header);
	my @sequences;
	$seqio->read_fasta($msa_path, \@$array_ref);

   # Create arrays of IDs and headers and figure out minimum size of alignment
	my $msa_length;
	my $seqs_ref;
	my $last_length = undef;
	foreach my $seq_ref (@$array_ref) {
		
		# Get sequence details
		my $seq_id    = $seq_ref->{sequence_id};
		my $header    = $seq_ref->{header};
		my $sequence  = $seq_ref->{sequence};
		#$seq_ref->{aln_seq} = $seq_ref->{sequence};
		unless ($sequence and $seq_id and $header) { die; }
		
		# Get start andstop coordinates
		$seq_ref->{aln_start} = $start;
		$seq_ref->{aln_stop}  = $stop;
		$seq_ref->{padded}    = 'true';
		#$self->do_seq_trim($seq_ref);

		# Check sequence length matches rest of alignment
		my $seq_len = length $sequence;
		my $raw_seq = $sequence;
		$raw_seq =~ s/-//g;
		my $raw_seq_len = length $sequence;
		$seq_ref->{raw_seq_length} = $raw_seq_len;

		unless ($last_length) {
			$stop = $seq_len;
			$last_length = $seq_len;
		}
		elsif ($seq_len ne $last_length) { 
			print "\n\t  GLUE alignment parsing error - unequal sequence lengths";	
			print "\n\t  This sequence length $seq_len, last sequence length $last_length";	
			die; 
		}
		$last_length = $seq_len;
		#print "\n\t length $seq_len last: $last_length";
		if ($start and $stop) {
			$msa_length = ($stop - $start) + 1;
		}
		else {
			unless ($msa_length) { $msa_length = $seq_len; }
			elsif ($seq_len != $msa_length) { 
				print "\n\t UNEQUAL sequence lengths $seq_len ne $msa_length\n\n";
				die; 
			}
		}		
	}
	$hash_ref->{msa_length} = $msa_length;
}

#***************************************************************************
# Subroutine:  link_sequences_and_data
# Description: attempt to link a fasta file to a tab-delimited data file 
#            - first column of the data file should correspond to sequence
#              headers in the FASTA file. 
#***************************************************************************
sub link_sequences_and_data {

	my ($self, $seqs_ref, $data_ref, $results_ref) = @_;

	# Sanity checking
	unless ($seqs_ref and $data_ref) { die; }
	
	# Set up data structures
	my %hash_sequences;
	my @seqs_without_data;
	my @data_without_seqs;
	
	my $topline = shift @$data_ref;
	chomp $topline;
	my %fields;
	my @fields = split("\t", $topline);
	my $i =0;
	foreach my $field (@fields) {
		$fields{$i} = $field;
		$i++;
	}

	# Add column headings row to 'data without seqs' array
	push (@data_without_seqs, $topline);

	# Index all sequence IDs in data file
	my %indexed_data;
	foreach my $seq_line (@$data_ref) {
		chomp $seq_line;
		my @line = split("\t", $seq_line);
		my $seq_id = shift @line;
		$indexed_data{$seq_id} = 1;
	}
	#$devtools->print_hash(\%indexed_data); die;
	
	# Index all sequences 
	my %indexed_seqs;
	my %link_hash;
	foreach my $seq_ref (@$seqs_ref) {
		# Get seq data
		my $id       = $seq_ref->{sequence_id};
		my $header   = $seq_ref->{header};
		my $sequence = $seq_ref->{sequence};
		$indexed_seqs{$id}    = $seq_ref; # Store in ref hash
		$link_hash{$header} = $id;
	}
	
	# Check all sequences in the alignment are in data file
	my $message = "Checking linking sequences & data...";
	#$io->show_output_message($message, $results_ref->{output_type});
	my $align_linked = 0;
	foreach my $seq_ref (@$seqs_ref) {
		
		# Get data for this sequence
		my $id     = $seq_ref->{sequence_id};
		my $header = $seq_ref->{header};
		my $seq_id = $link_hash{$header};
		if  ($indexed_data{$header}) {
			$align_linked++;
		}
		else {
			my $message = " No data for '$header' ";
			print "\n$message";
			#$io->show_output_message($message, $results_ref->{output_type});
			push (@seqs_without_data, $seq_ref);
		}
	}
	
	# Check all sequences in data file are in the alignment
	my $data_linked = 0;
	$i = 0; 
	foreach my $seq_line (@$data_ref) {
		$i++;
		#print $seq_line;die;
		chomp $seq_line;
		my @line = split("\t", $seq_line);
		my $header = shift @line;
		my $seq_id  = $link_hash{$header};
		if ($seq_id) {	
			
			# Add data to sequences
			my %data;
			my $x = 0;
			my $seq_ref = $indexed_seqs{$seq_id};
			foreach my $value (@line) {
				$x++;
				my $field = $fields{$x};
				unless ($field) { die; }
				$data{$field} = $value;
			}
			$seq_ref->{data} = \%data;
			$data_linked++;
		
		}
		else {
			my $message = " Sequence '$header' not in the FASTA file";
			print "\n\t$message";
			#$io->show_output_message($message, $results_ref->{output_type});
			push (@data_without_seqs, "$seq_line\n");
		}
	}

	# Record the results of this attempt to link
	$results_ref->{linked}            = $data_linked; # Number of linked seqs
	$results_ref->{seqs_without_data} = \@seqs_without_data;
	$results_ref->{data_without_seqs} = \@data_without_seqs;
	$results_ref->{header_id_link}    = \%link_hash;

	# Sanity checking
	#my $linked_ref = $results_ref->{linked};
	#unless ($linked_ref) {
	#	my $error = "\n\t No sequences were linked;\n\n"; 
	#	die $error;
	#	return 0;
	#}
}

#***************************************************************************
# Subroutine:  by_number
# Description: by number - for use with perl 'sort'  (cryptic but works) 
#***************************************************************************
sub by_number { $a <=> $b }	

############################################################################
# EOF
############################################################################
