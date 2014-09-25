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
use Base::FileIO;
use Base::SeqIO;
use Base::DevTools;
use Base::Console;
use Base::Sequence; # For performing basic sequence manipulations

############################################################################
# Globals
############################################################################

# Create base objects
my $fileio     = FileIO->new();
my $seqio      = SeqIO->new();
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
		process_id           => $parameter_ref->{process_id},
		output_type          => $parameter_ref->{output_type},
		
		# Member classes
		blast_obj 			 =>	$parameter_ref->{blast_obj},
		
		# Paths
		output_path          => $parameter_ref->{output_path},
		header_path          => $parameter_ref->{header_path},
		refseq_use_path      => $parameter_ref->{refseq_use_path}, 
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# SECTION: Handler functions
############################################################################

#***************************************************************************
# Subroutine:  run_reformat_tools_cmd_line
# Description: hand off to FASTA utilities
#***************************************************************************
sub run_reformat_tools_cmd_line {

	my ($self, $infile, $mode) = @_;

	#====# FILTER
	#====# SORT	
	#===# CONVERT 	
	if ($mode eq 1) {  # FASTA to delimited txt table 
		print "\n\t # Converting file '$infile' from FASTA to data";
		my @fasta;
		$seqio->read_fasta($infile, \@fasta);
		my $delimiter;
    	my $question = "\n\t What is the delimiter (comma\/tab)";
		my @choices  = qw [ c t ];
		my $choice = $console->ask_simple_choice_question($question, \@choices);
		if ($choice eq 'c') { $delimiter = ',';  }
		if ($choice eq 't') { $delimiter = "\t"; }
		my @data;
		$self->fasta_to_delimited(\@fasta, \@data, $delimiter);
		my $outfile = $infile . '.tsv.txt';
		$fileio->write_file($outfile, \@data);
		print "\n\t # File '$outfile' created";
	}
	elsif ($mode eq 2) { # Delimited to FASTA
		my @data;
		$fileio->read_file($infile, \@data);
		my @fasta;
		my $delimiter;
    	my $question = "\n\t What is the delimiter (comma\/tab)";
		my @choices  = qw [ c t ];
		my $choice = $console->ask_simple_choice_question($question, \@choices);
		if ($choice eq 'c') { $delimiter = ',';  }
		if ($choice eq 't') { $delimiter = "\t"; }
		my $outfile = "$infile.fas";
		$self->delimited_to_fasta(\@data, \@fasta, $delimiter);
		$fileio->write_file($outfile, \@fasta);
		print "\n\t # File '$outfile' created";
		#$devtools->print_array(\@data); die;
	}
	elsif ($mode eq 3) {  # Genbank to FASTA + Data
		print "\n\t # Converting file '$infile' from GenBank format to FASTA+DATA";
		my @fasta;
		my @data;
		$self->genbank_to_fasta_and_data($infile, \@fasta, \@data);
		#$devtools->print_array(\@fasta);
		#$devtools->print_array(\@data);
		my $outseqs = $infile . '.fas';
		print "\n\t Writing sequences to file '$outseqs'";
		$seqio->write_fasta($outseqs, \@fasta);
		print "\n\t # File '$outseqs' created";
		my $outdata = $infile . '.txt';
		print "\n\t Writing data to file '$outdata'";
		$seqio->write_delimited($outdata, \@data);
		print "\n\t # File '$outdata' created";
	}
	elsif ($mode eq 4) {  # FASTA to NEXUS
		print "\n\t # Converting file '$infile' from FASTA to NEXUS format";
		my @fasta;
		$seqio->read_fasta($infile, \@fasta);
		my $num_taxa = scalar @fasta;
		unless ($num_taxa) { die "\n\t NO SEQUENCES in '$infile'\n\n\n"; }
		my @nexus;
		$self->fasta_to_nexus(\@fasta, \@nexus);
		my $outfile = $infile . '.nex';
		$fileio->write_file($outfile, \@nexus);
		print "\n\t # File '$outfile' created";
	}
	elsif ($mode eq 5) {  # FASTA TO PHYLIP
		print "\n\t # Converting file '$infile' from FASTA to PHYLIP format";
		my @fasta;
		$seqio->read_fasta($infile, \@fasta);
		my $num_taxa = scalar @fasta;
		unless ($num_taxa) { die "\n\t NO SEQUENCES in '$infile'\n\n\n"; }
		my @phylip;
		$self->fasta_to_phylip(\@fasta, \@phylip);
		my $outfile = $infile . '.phy';
		$fileio->write_file($outfile, \@phylip);
		print "\n\t # File '$outfile' created";
	}
}

#***************************************************************************
# Subroutine:  run_refseq_tools_cmd_line
# Description: hand off to RESEQ utilities
#***************************************************************************
sub run_refseq_tools_cmd_line {

	my ($self, $infile, $mode) = @_;

	if ($mode eq 1) {  # GenBank to REFSEQ
		my @genbank;
		my @refseq;
		$fileio->read_file($infile, \@genbank);
		$self->genbank_to_refseq(\@genbank, \@refseq);
	}
	elsif ($mode eq 2) {

	}
	elsif ($mode eq 3) {
		my $refseq_parser = RefSeqParser->new();
		my $refseq;
		$refseq;

	}
}

#***************************************************************************
# Subroutine:  run_extract_tools_cmd_line
# Description: hand off to data utilities
#***************************************************************************
sub run_extract_tools_cmd_line {

	my ($self, $infile, $mode) = @_;

	$self->extract_random_seqs();
	$self->extract_seqs();
	$self->extract_fasta_seqs();

}

#***************************************************************************
# Subroutine:  run_sort_tools_cmd_line
# Description: hand off to sort utilities
#***************************************************************************
sub run_sort_tools_cmd_line {

	my ($self, $infile, $mode) = @_;

	if ($mode eq 1) {  # GenBank to FASTA
		$self->sort_seqs_by_length();
	}
	elsif ($mode eq 2) {
		$self->sort_by_identity();
	}
	elsif ($mode eq 3) {
		$self->sort_sequences_using_tabdelim_table();
	}
}

#***************************************************************************
# Subroutine:  run_data_tools_cmd_line
# Description: hand off to data utilities
#***************************************************************************
sub run_data_tools_cmd_line {

	my ($self, $infile, $mode) = @_;

	$self->link_sequences_and_data();
	$self->combine_data();

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

	my ($self, $fasta_ref, $phylip_ref, $table_ref) = @_;
	
	my $set_id_len = 20;

	my $phycount;
	my $seq_len;
	foreach my $sequence_ref (@$fasta_ref) {
		
		$phycount++;
		my $sequence_id = $sequence_ref->{header};
		my $sequence    = $sequence_ref->{sequence};
		unless ($seq_len) {
			$seq_len = length $sequence;
		}
		else {
			unless ($seq_len eq length $sequence) {
				die "\n\t sequence $sequence_id is a different length\n
				       \t check FASTA sequences are aligned\n\n";
			}
		}
	
		# Create PHYLIP id
		my $phy_id = $phycount . '_'; 
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
		
		# store id relationship in translation table 
		my $id_pair = "$sequence_id\t$phy_id\n";
		push (@$table_ref, $id_pair);
	
	}
	
	# Create PHYLIP taxa and characters header
	my $num_taxa  = $phycount;
	my $num_chars = $seq_len;
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

			push(@$data_ref, \%data);
			$sequence = uc $sequence;
			my $header = $accession;
			my $seq_obj = Sequence->new($sequence, $header, $seq_id);
			push(@$seq_ref, $seq_obj);

			# RESET
			$sequence    = '';
			$seq_id      = '';
			$seq_date    = '';
			$iso_country = 'NK';
			$genotype    = undef;

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
# SECTION: REFSEQ functions
############################################################################

#***************************************************************************
# Subroutine:  genbank_to_refseq
# Description: GenBank to REFSEQ
#***************************************************************************
sub genbank_to_refseq {

	my ($self, $gb_ref, $refseq_ref) = @_;

	# TODO
	foreach my $line (@$gb_ref) {
		
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
}

#***************************************************************************
# Subroutine:  sort_by_identity 
# Description: 
#***************************************************************************
sub sort_by_identity {

	my ($self, $seqs_ref, $sorted_ref, $output_ref) = @_;

	# Get report directory
	my $report_dir  = $self->{report_dir};
	my $output_type = $self->{output_type};
	unless ($report_dir and $output_type) { die "Error - no report directory defined";}

	# Get the GLUE reference sequence information
	my $refseq_seq;   
	my $refseq_name;
	if ($self->{refseq_name}) {
		my $refseq_name = $self->{refseq_name};	
		my $parser_obj  = RefSeqParser->new();
		my $flat_dir = $self->{refseq_use_path};
		unless ($flat_dir)  {die; }
		# Read reference sequence 
		my $refseq_file = $flat_dir . $refseq_name;
		my %params;
		$parser_obj->parse_refseq_flatfile($refseq_file, \%params);
		my $refseq = RefSeq->new(\%params);
		$refseq_seq   = $self->{refseq}->{sequence};
	}
	elsif ($self->{refseq_filename}) {
		die;
	}
	elsif ($self->{query_file_refseq}) {
		my $seq_ref = shift(@$seqs_ref);
		$refseq_seq   = $seq_ref->{sequence};
		$refseq_name  = $seq_ref->{name};
	}
	else { die; }

	# Create FASTA file for BLAST input
	my $lib_path = $report_dir . 'refseq.txt';
	unless ($refseq_seq and $refseq_name) { die "Error - no reference sequence found";}
	my $refseq_fasta = ">$refseq_name\n$refseq_seq\n";	
	$fileio->write_text_to_file($lib_path, $refseq_fasta);
	
	# Configure the BLAST tool
	my $blast_tool = BLAST_Tool->new($self);
	$blast_tool->{result_path} = $report_dir;
	$blast_tool->{blast_alg} = 'blastn';
	$blast_tool->{blast_genome_lib_path} = $lib_path;

	# BLAST each sequence
	foreach my $seq_ref (@$seqs_ref) {
		
		# do the assign
		my $id = $seq_ref->{sequence_id};
		my @hits;
		my %settings;
		$settings{blast_type} = 'blast2seq';
		$settings{set_params} = 'TRUE';
		$settings{word_size} = 11; 
		$settings{evalue}    = 10;
		$settings{penalty}   = -3;
		$settings{reward}    = 2;
		$settings{gapopen}   = 5;
		$settings{gapextend} = 2;
		$settings{outfmt}    = 7;
		my $result_path = $blast_tool->assign($seq_ref, \@hits, \%settings);
		#$devtools->print_array(\@hits); #die;
		$seq_ref->{result_path} = $result_path;
		$seq_ref->{result} = \@hits;
	}
	
	# Index the sequences by bitscore
	my %bit_scores;
	my %sequences;
	foreach my $seq_ref (@$seqs_ref) {
	
		my $id       = $seq_ref->{sequence_id};
		my $header   = $seq_ref->{header};
		my $sequence = $seq_ref->{sequence};
		
		# Get the best match from this file
		my $hits_ref      = $seq_ref->{result};
		my $top_match     = shift @$hits_ref;
		my $bitscore;
		unless ($top_match) { 
			#print "\n\t<BR><BR> header $header NO MATCH"; 
			$bitscore = '0';
			$seq_ref->{assigned_to}   = '-';	
			$seq_ref->{subject_start} = '-';
			$seq_ref->{subject_end}   = '-';
		}
		else {
			$bitscore = $top_match->{bit_score};	
			$seq_ref->{assigned_to}   = $top_match->{scaffold};	
			$seq_ref->{subject_start} = $top_match->{aln_start};
			$seq_ref->{subject_end}   = $top_match->{aln_stop};
		}
		unless ($bitscore) { $bitscore = '0'; }
		$seq_ref->{bitscore} = $bitscore;
		if ($bit_scores{$bitscore}) {
			my $seqs_ref = $bit_scores{$bitscore};;
			push (@$seqs_ref, $seq_ref);
		}
		else {
			my @seqs;
			push (@seqs, $seq_ref);
			$bit_scores{$bitscore} = \@seqs;
		}
	}

	#unless($self->{output_type} eq 'html') {
	#print "\n\tSEQ LENGTHS IN DESCENDING NUMERIC ORDER:\n";
	my @bitscores = sort by_number keys (%bit_scores);
	#my @sorted;
	foreach my $score (@bitscores) {
		
		my $seqs_ref =  $bit_scores{$score};
		my $num_seq  = scalar @$seqs_ref;
		#print "\n\t $num_seq sequences with bitscore $score";
		foreach my $seq_ref (@$seqs_ref) {	
			#$devtools->print_hash($seq_ref);
			unshift (@$sorted_ref, $seq_ref);
		}
	}

	# Write the FASTA file (sorted seqs)
	my @fasta;
	my $outpath = $report_dir . 'sorted_seqs.fas';
	$fileio->write_file($outpath, \@fasta);

	# Create the table
	my @table;
	$self->create_sorted_by_identity_table($seqs_ref, $sorted_ref, \@table);
		
	# Add the links if in HTML
	my $filelink  = "<p>Click <a href='$outpath'>here</a> to retrieve the sorted ";
	$filelink .= "sequences <BR><BR><BR>";
	push (@$output_ref, $filelink);
	push (@$output_ref, @table);
	
	# Close off the page
	push (@$output_ref,'<div id="menubottom"></div>');
}

#***************************************************************************
# Subroutine:  sort_sequences_using_tabdelim_table
# Description: Sort sequences using tab data file field
#***************************************************************************
sub sort_sequences_using_tabdelim_table {

	my ($self, $seqs_ref, $data_ref) = @_;

	my $delimiter = "\t";
	
	# get the indices
	my $question = "\n\t Which column contains the sequence IDs \[1st column=0\]?";
	my $id_index = ask_int_question($question);
	$question = "\n\t Which column contains the data to sort on \[1st column=0\]?";
	my $sub_index = ask_int_question($question);

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

#***************************************************************************
# Subroutine:  create_sorted_seqcutoff_subset
# Description: sort a fasta file of sequences in order of descending length
#              and discard all below size cut-off
#***************************************************************************
sub create_sorted_seqcutoff_subset {

	my ($self, $file, $fasta_ref, $cutoff) = @_;

	# Ensure sanity
	unless ($cutoff) { $cutoff = 1; }

	# Read the sequence FASTA
	my %sequences;
	$seqio->read_fasta_file_to_hash($file, \%sequences);
	
	# Index the sequences by size
	my $id;
	my $sequence;
	my %seq_lengths;
	while ( ( $id, $sequence ) = each %sequences ) {
	
		my $adjusted_seq = $sequence;
		
		my $length;
		# Don't count gaps
		$adjusted_seq =~ s/-//g;
		$length = length $adjusted_seq;
		
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


	#my @fasta;
	my @lengths = sort by_number keys (%seq_lengths);
	foreach my $length (@lengths) {
		
		if ($length < $cutoff) { next; }

		my $seqs_ref = $seq_lengths{$length};
	   	my $num_seq = scalar @$seqs_ref;
		#print "\n\tLength $length': $num_seq sequences";
		foreach my $id (@$seqs_ref) {
			my $sequence = $sequences{$id};
			my $fasta    = ">$id\n$sequence\n"; 
			unshift (@$fasta_ref, $fasta);
		}
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

############################################################################
# SECTION: Sequence extraction functions
############################################################################

#***************************************************************************
# Subroutine:  extract_random_seqs
# Description: extract a random subset of sequences from a FASTA file 
#***************************************************************************
sub extract_random_seqs {

	# To do: show about box for this utility
	#show_extract_fxn_about_box();
	my ($self) = @_;
	
	# Open the FASTA file
	my @fasta;
	my $question = "\n\t What is the name of the fasta file?";
	my $fasta_file = $console->ask_input_file_question($question, \@fasta);
	unless ($fasta_file) { return; }

	
	# Create lookup hash
	my $i = 0;
	my $sequence = undef;
	my $header   = undef;
	
	my %sequences;
	foreach my $line (@fasta) { 

		chomp $line;	
		if    ($line =~ /^\s*$/)   { next; } # discard blank line
		elsif ($line =~ /^\s*#/)   { next; } # discard comment line 
		elsif ($line =~ /^>/) {
		
			if ($sequence) {
				$i++;
				$sequences{$i} = "$header\n$sequence\n\n";;
				# Reset variables
				$sequence = undef;
				$header   = $line;
			}
			else {
				$line =~ s/\s+//g; # Remove whitespace
				$header = $line;
			}
		}
		# Any line that isn't a header or blank should be sequence
		else { $sequence .= $line; }
	}
	
	# Get last sequence
	$i++;
	$sequences{$i} = "$header\n$sequence\n\n";;

	# get how many seqs to extract
	$question = "\n\t How many sequences do you want to extract?";
	my $number = $console->ask_int_with_bounds_question($question, 1, ($i - 1));
	
	# seed random number generator
	srand(time|$$);
	
	my @seqs;
	my %chosen_seqs;
	my $num_chosen;
	do {

		my $index;
		do {
			$index = int(rand(($number + 1)));
		} until ($index > 0);
		
		print "\n\t chose index $index";
		unless ($chosen_seqs{$index}) { 
			
			my $sequence = $sequences{$index};
			print "\n\t apparently capturing $sequence";
			push (@seqs, $sequence);
			
			$chosen_seqs{$index} = 'yes';
			$num_chosen++;
		}

	} until ($num_chosen eq $number);

	print "\n\t apparently $num_chosen = $number";

	# write out the sequences
	my $file = 'random_from_' . $fasta_file;
	$fileio->write_file($file, \@seqs);
}


#***************************************************************************
# Subroutine:  extract_fasta_seqs
# Description: 
#***************************************************************************
sub extract_fasta_seqs {

	my ($self) = @_;
	
	# show about box for this utility
	#show_extract_fxn_about_box();
	# create counting variables
	my $num_matched = 0;
	my $num_unmatched = 0;

	# open the first file
	my @ids;
	my $question = "\n\t What is the name of the file with the ID list?";
	my $id_file = $console->ask_input_file_question($question, \@ids);
	unless ($id_file) { return; }

	$question = "\n\t Which column contains the ID field?";
	my $index = $console->ask_int_question($question);
	
	# open the second file
	my @fasta;
	$question = "\n\t What is the name of the fasta file?";
	my $fasta_file =$console-> ask_input_file_question($question, \@fasta);
	unless ($fasta_file) { return; }

	# create lookup hash
	my %seen_in_id_list;
	foreach my $line (@ids) { 
	
		chomp $line;
		my @line_bits = split("\t", $line);
		my $id = $line_bits[$index];
		$id =~ s/\s//g;
		$seen_in_id_list{$id} = 1; 
	}

	#print_hash(\%seen_in_id_list);

	my $i = 0;
	my $header;
	my $sequence;
	my $initialised;
	my @captured_seqs;
	my @rejected_seqs;
	my @matched_ids;
	foreach my $line (@fasta) {

		chomp $line;
		if ($line =~ /^(\s)*$/) { next; };	# skip blank lines
		if ($line =~ /^>/)      { 
			
			$i++;
			if ($initialised) {
				my @header_elements = split("\t", $header);
				my $id = $header_elements[0];
				$id =~ s/^>//g;
				
				# if the id (and others?) match then capture
				print "\n$i:\t\t Checking to see if sequence $id should be extracted";
				if ($seen_in_id_list{$id}) { 
					print "\n\t Matched $id";

					my $captured_seq = $header . "\n" . $sequence . "\n\n";
					push (@captured_seqs, "$captured_seq\n");
					push (@matched_ids, "$id\twas matched to $header\n");
					$num_matched++;
				}
				else { 
					print "\n\t Couldn't find ID $id";
					my $rejected_seq = $header . "\n" . $sequence . "\n\n";
					push (@rejected_seqs, "$rejected_seq\n");
					$num_unmatched++; 
				}
			}
			
			# set state variables, and capture new header 
			unless ($initialised) {$initialised = 'true'; }
			$sequence = undef;
			$header   = $line;
		}

		# any line that isn't a header or blank should be sequence
		else { $sequence .= $line; }
	}

	# get last sequence
	my @header_elements = split("\t", $header);
	my $id = $header_elements[0];
	$id =~ s/^>//g;
	print "\n$i:\t\t Checking to see if sequence $id should be extracted";
	
	if ($seen_in_id_list{$id}) { 
		print "\n\t Matched $id";
		my $captured_seq = $header . "\n" . $sequence . "\n\n";
		push (@captured_seqs, "$captured_seq\n");
		push (@matched_ids, "$id\twas matched to $header\n");
		$num_matched++;
	}
	else { 
		my $rejected_seq = $header . "\n" . $sequence . "\n\n";
		push (@rejected_seqs, "$rejected_seq\n");
		$num_unmatched++; 
	}

	# show stats 
	print "\n\t ***************** RESULTS ******************";
	print "\n\t       Number matched   = $num_matched";
	print "\n\t       Number unmatched = $num_unmatched";
	print "\n\t ********************************************";

	# write out the extracted sequences
	my $output_file = 'extracted_from_' . $fasta_file;
	$fileio->write_file($output_file, \@captured_seqs);

	# write out the rejected sequences
	my $reject_file = 'unmatched_from_' . $fasta_file;
	$fileio->write_file($reject_file, \@rejected_seqs);

	##### SUPERFLUOUS
	# write out the extracted sequence ids
	my $verify = 'extracted.txt';
	$fileio->write_file($verify, \@matched_ids);

}

#***************************************************************************
# Subroutine:  extract_seqs 
# Description: 
#***************************************************************************
sub extract_seqs {

	my ($self, $id_file, $seq_file) = @_;

	open IDFILE, "<$id_file" or die "\n\tCan't open $id_file\n";
	my %id_keys;
	while ( <IDFILE> ) {
		my $line1 = $_;
		chomp $line1;
		print "\n\t $line1";
		$id_keys{$line1} = 1;
	}
	

	open SEQFILE, "<$seq_file" or die "\n\tCan't open $seq_file\n";
	my $id  = undef;
	my $seq = '';
	my @extracted;
	my $i = 0;
	my $j = 0;
	my $k = 0;
	while ( <SEQFILE> ) {
	
		my $line = $_;
		chomp $line;
		if ($line =~ /^>/) {
			
			$line =~ s/^>//g;
			if ($id and $seq) {
				my $fasta = ">$id\n$seq\n\n";
				$i++;
				print "\n\t got contig number $i: $id";
				#print "\n\t $fasta";
				push (@extracted, $fasta);
				$id  = undef;
				$seq = undef;
			}
			
			if ($id_keys{$line}) {
				$j++;
				$id = $line;	
				print "\n\t encountered contig number $j: $id";
			}
			
		}
		elsif ($id) {
			
			unless ($seq) {
				$seq = $line;
			}
			else {
				$seq .= $line;
			}
		}
	}
	
	if ($id and $seq) {
		my $fasta = ">$id\n$seq\n\n";
		push (@extracted, $fasta);
	}

	my $outfile = $id_file . '.extracted.fas';
	$fileio->write_file($outfile, \@extracted);
}

#***************************************************************************
# Subroutine:   split_vglue_ref_file
# Description:  
#***************************************************************************
sub split_vglue_ref_file  {

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
		print "\n\t Checking ID $id in $file_a\t";
		
		if ($seen_in_b{$id}) { 
			my $combined_line = $line;
			$combined_line .= "\t$seen_in_b{$id}\n";
			#print $combined_line;
			push (@in_both, $combined_line); 
		}
		else { 
		
			# Output for debugging etc
			print "....couldn't find ID $id in $file_b";
			# ----- HACK for if you want the line from 1st as well if no match
			#push (@in_both, "$line\n"); 
			push (@in_first_only, "$id\n");
		}
	}
}

#***************************************************************************
# Subroutine:  link_sequences_and_data
# Description: attempt to link a fasta file to a tab-delimited data file 
#            - first column of the data file should correspond to sequence
#              headers in the FASTA file. 
#***************************************************************************
sub link_sequences_and_data {

	my ($self, $seqs_ref, $data_ref, $results_ref) = @_;

	my $io;
	# Sanity checking
	unless ($seqs_ref and $data_ref) { die; }
	
	# Set up data structures
	my %hash_sequences;
	my @seqs_without_data;
	my @data_without_seqs;
	
	# Add column headings row to 'data without seqs' array
	my $topline = @$data_ref[0];
	push (@data_without_seqs, $topline);

	# Index all sequence IDs in data file
	my %indexed_data;
	my $i = 0;
	foreach my $seq_line (@$data_ref) {
		$i++;
		unless ($i > 1) { next; } # Skip top line (column headings)
		chomp $seq_line;
		my @line = split("\t", $seq_line);
		my $seq_id = shift @line;
		$indexed_data{$seq_id} = 1;
	}
	
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
	$io->show_output_message($message, $results_ref->{output_type});
	
	my $align_linked = 0;
	foreach my $seq_ref (@$seqs_ref) {
		
		# Get data for this sequence
		my $id       = $seq_ref->{sequence_id};
		my $header   = $seq_ref->{header};
		
		my $seq_id  = $link_hash{$header};
		if  ($indexed_data{$header}) {
			$align_linked++;
		}
		else {
			my $message = " No data for '$header' ";
			$io->show_output_message($message, $results_ref->{output_type});
			push (@seqs_without_data, $seq_ref);
		}
	}
	
	# Check all sequences in data file are in the alignment
	my $data_linked = 0;
	$i = 0; 
	foreach my $seq_line (@$data_ref) {
		$i++;
		unless ($i > 1) { next; } # Skip top line (column headings)
		chomp $seq_line;
		my @line = split("\t", $seq_line);
		my $header = $line[0];
		my $seq_id  = $link_hash{$header};
		if  ($indexed_seqs{$seq_id}) {
			$data_linked++;
		}
		else {
			my $message = " Sequence '$header' ('$seq_id') in the FASTA file";
			$io->show_output_message($message, $results_ref->{output_type});
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

	# DEBUG
	#$devtools->print_hash($results_ref); die;
	#unless ($data_linked eq $align_linked) {
	#	die "\n\t ERROR in matching algorithm\n\n";
	#}
}

############################################################################
# SECTION: Summarising data & sequences
############################################################################

#***************************************************************************
# Subroutine:  count_fasta_sequences
# Description: count the number of sequences in a fasta file
# Arguments:   $file: the name of the file to read
# Returns:     $num_seqs:  number of sequences
#***************************************************************************
sub count_fasta_sequences {

	my ($self, $file) = @_;
	
	my %sequences;
	$seqio->read_fasta_file_to_hash($file, \%sequences);
	my @sequences = keys %sequences;
	my $num_seqs  = scalar @sequences;

	return $num_seqs;
}

#***************************************************************************
# Subroutine:  show_fasta_sequence_lengths
# Description: show the lengths of a sequence in a FASTA file 
# Arguments:   $file: the name of the file to read
#***************************************************************************
sub show_fasta_sequence_lengths {

	my ($self, $file) = @_;

	my %sequences;
	$seqio->read_fasta_file_to_hash($file, \%sequences);
	my @sequences = keys %sequences;
	my $num_seqs  = scalar @sequences;

	my $id;
	my $sequence;
	print "\n\n";
	while ( ( $id, $sequence ) = each %sequences) {
		my $seq_length = length $sequence;
		print "\n\t #### Sequence $id\t\tLength: $seq_length";
	}
	print "\n\n";

	print "\n\t *** There are $num_seqs uniquely";
	print " identified sequences in the file $file";

}

############################################################################
# SECTION: Web interface
############################################################################

#***************************************************************************
# Subroutine:  web_run_data_tools
# Description: hand off to FASTA utilities
#***************************************************************************
sub web_run_data_tools {

	my ($self, $infile, $mode) = @_;

	# Create the folder  (todo: use call to fxn)
	#my $output_type = $self->{output_type};
	#my $process_id  = $self->{process_id};
	#my $output_path = $self->{output_path};
	#unless ($output_type and $process_id and $output_path) { die; }
	#my $report_dir  = $output_path . $process_id . '/';
	#my $mkdir_cmd   = "mkdir $report_dir";
	#system $mkdir_cmd;
	#$self->{report_dir} = $report_dir;
	#my $message = " created report directory $report_dir";
	#$io->show_output_message($message, $output_type);

	# Raw sequences
	#my $output_ref = $self->{output};
	#my $raw_path = $report_dir . 'raw.fas';
	#my $raw_sequences_ref = $self->{sequences};
	#$seqio->write_fasta($raw_path, $raw_sequences_ref);
	#my $rawlink  = "<p>Click <a href='$raw_path'>here</a> to retrieve the raw ";
	#$rawlink .= "sequences.<BR><BR>";
	#push (@$output_ref, $rawlink);

	# Check if we are converting only
	#my $skip_fasta_steps = undef;
	#if ($self->{convert_method}) {
	#	if ($self->{convert_method} eq 'delimited_to_fasta') {
	#		$skip_fasta_steps = 1;
	#	}
	#}

	# Convert sequences
	my $report_dir = $self->{report_dir};
	
	# Get sequence objects in array
	my $sequences_ref = $self->{sequences};
	unless ($sequences_ref) { die; }

	#my @converted;
	#if ($self->{convert_format} and $self->{convert_method} eq 'fasta_to_delimited') {  
	#	$self->fasta_to_delimited($sequences_ref, \@converted);		
	#}	
	#elsif ($self->{convert_format} and $self->{convert_method} eq 'delimited_to_fasta') {  
	#	
	#	if ($self->{output_type} eq 'command_line') {	
	#		my @delimited;
	#		my $tabdelim_path = $self->{tabdelim_path};
	#		my $delimiter     = $self->{delimiter};
	#		unless ($delimiter)      { die; }
	#		unless ($tabdelim_path)  { die "No input file defined\n\n"; }
	#		$fileio->read_input_file($tabdelim_path, \@delimited);		
	#		$self->delimited_to_fasta(\@delimited, \@converted, $delimiter);	
	#	}
	#	else {
	#		my $delimited_ref = $self->{delimited};  
	#		my $delimiter     = $self->{delimiter};
	#		unless ($delimiter and $delimited_ref) { die; }
	#		$self->delimited_to_fasta($delimited_ref, \@converted, $delimiter);		
	#	}
	#}
	#====# FILTER
	my $sequences = $self->{sequences};
	my $outpath;
	if ($self->{filter_by_sequence} or $self->{filter_by_header}) {  
		$self->filter_sequences();
		$self->{filtered} = 1;
	}
	#====# SORT	
	if ($self->{sort_sequences}) {
		$self->sort_sequences();
		$self->{sorted} = 1;
	}

	#===# CONVERT 	
	if ($self->{convert_format}) {
		$self->convert_format($infile, $mode);
	}

	# Write the output
	#if ($self->{output_type} eq 'command_line') {	
	#	my $fasta_path = $self->{fasta_path};	
	#	my @path = split(/\//, $fasta_path);
	#	my $filename = pop @path;
	#	my $outfile = $report_dir .  $filename . '.converted.fa';
	#	print "\n\t ## Converted sequences written to $outfile";
	#	$fileio->write_file($outfile, \@converted);
	#}
	#else {
	#	my $output_ref = $self->{output};
	#	unless ($output_ref) { die; }
	#	my $outfile = 'converted_seqs.txt';
	#	my $outpath = $report_dir . $outfile;
	#	$fileio->write_file($outpath, \@converted);
	#	my $filelink  = "<p>Click <a href='$outpath'>here</a> to retrieve the";
	#	if ($self->{convert_method} ne 'delimited_to_fasta') {
	#		if ($self->{filtered}) { $filelink .= " filtered,"; }
	#		if ($self->{sorted})   { $filelink .= " sorted,";   }
	#	}
	#	$filelink .= " converted sequences <BR><BR>";
	#	push (@$output_ref, $filelink);
	#}
}

#***************************************************************************
# Subroutine:  web_sort_sequences 
# Description: handler fxn handing off to specific sorting fxns 
#***************************************************************************
sub web_sort_sequences {

	my ($self) = @_;

	# Get sequence objects in array
	my $sequences_ref = $self->{sequences};
	unless ($sequences_ref) { die; }

	# Get report directory
	my $report_dir = $self->{report_dir};
	my @sorted;

	# Hand off to method
	#if ($self->{sort_sequences} and $self->{sort_seqs_method} eq 'sort_by_length') {  
	if ($self->{sort_sequences}) {
		# Sort by sequence length
		my $minimum = $self->{minimum_seqlen};
		unless ($minimum) { $minimum = 1; }
		$self->sort_seqs_by_length($sequences_ref, \@sorted, $minimum);
		#$devtools->print_web_array(\@sorted);
		#exit;
	}
	# Sort by identity
	#elsif ($self->{sort_sequences} and $self->{sort_seqs_method} eq 'sort_by_identity') {  
	#	my @sorted;
	#	my @output;
	#	my $ref;
	#	$self->sort_by_identity($sequences_ref, \@sorted);
	#}   
	# TODO: implement these 
	# Sort by header
	#elsif ($self->{sort_sequences} eq 'sort_by_header') {
	#	die "UNIMPLEMENTED";	
	#}
	# Sort by header element
	#elsif ($self->{sort_sequences} eq 'sort_by_header_element') {
	#	die "UNIMPLEMENTED";	
	#}
	#else {
	#	if ($self->{output_type} eq 'command_line') {
	#		die $self->{usage};
	#	}
	#}
	
	# Write the output
	if ($self->{output_type} eq 'command_line') {
		my $fasta_path = $self->{fasta_path};	
		my @path = split(/\//, $fasta_path);
		my $filename = pop @path;
		my $outfile = $report_dir .  'sorted_seqs.fas';
		#my $outfile = $report_dir .  $filename . '.sorted_seqs.fa';
		print "\n\t ## Sorted sequences written to $outfile";
		#$fileio->write_file($outfile, \@sorted);
	}
	else {
		my $output_ref = $self->{output};
		unless ($output_ref) { die; }
		my $outfile = 'sorted_seqs.fas';
		my $outpath = $report_dir . $outfile;
		$fileio->write_file($outpath, \@sorted);
		my $filelink  = "<p>Click <a href='$outpath'>here</a> to retrieve the";
		if ($self->{filtered}) {
			$filelink .= " filtered,";
		}
		$filelink .= " sorted sequences <BR><BR>";
		push (@$output_ref, $filelink);
	}
}

#***************************************************************************
# Subroutine:  web_filter_sequences 
# Description: handler fxn handing off to specific filtering fxns 
#***************************************************************************
sub web_filter_sequences {

	my ($self) = @_;

	# Get sequence objects in array
	my $output_ref = $self->{output};
	unless ($output_ref) { die; }

	# Set up 
	my $report_dir = $self->{report_dir};
	my $outfile = 'filtered_seqs.fa';

	# Create flags
	my $by_sequence = undef;
	if ($self->{filter_by_sequence} and $self->{minimum_seqlen}) {
		my $sequences_ref = $self->{sequences};
		unless ($sequences_ref) { die; }
		my $minimum = $self->{minimum_seqlen};
		my @filtered;
		$self->sort_seqs_by_length($sequences_ref, \@filtered, $minimum);

		# convert to sequence objects
		my @filtered_seqs; 	
		$seqio->convert_fasta(\@filtered, \@filtered_seqs);
		#$devtools->print_array(\@filtered_seqs); die;
		my $i = scalar @filtered_seqs;
		$self->{sequences} = \@filtered_seqs;	
	}	
	if ($self->{filter_by_header}) {
		my $sequences_ref = $self->{sequences};
		my $word  = $self->{exclude_word};
		my $rule  = $self->{exc_inc_rule};
		my $case  = $self->{ignore_case};
		my @filtered;
		$self->filter_by_header($sequences_ref, \@filtered, $word, $rule, $case);
		# convert to sequence objects
		my @filtered_seqs; 	
		$seqio->convert_fasta(\@filtered, \@filtered_seqs);
		$self->{sequences} = \@filtered_seqs;	
	}

	# Write the output
	if ($self->{output_type} eq 'command_line') {
		my $fasta_path = $self->{fasta_path};	
		my @path = split(/\//, $fasta_path);
		my $filename = pop @path;
		my $outfile = $report_dir .  $filename . '.filtered.fa';
		my $sequences_ref = $self->{sequences};
		$seqio->write_fasta($outfile, $sequences_ref);
		print "\n\t ## Filtered sequences written to $outfile";
	}
	else {
		my $outpath = $report_dir . $outfile;
		my $sequences_ref = $self->{sequences};
		$seqio->write_fasta($outpath, $sequences_ref);
		my $filelink  = "<p>Click <a href='$outpath'>here</a> to retrieve the filtered ";
		$filelink .= "sequences.<BR><BR>";
		push (@$output_ref, $filelink);
	}	
}

#***************************************************************************
# Subroutine:  run_refseq_function 
# Description: top level handler 
#***************************************************************************
sub run_refseq_function {

	my ($self, $mode, $file, $outpath) = @_;

	$self->show_title();
	
	# USAGE Statement
   	my $USAGE .= "\t  r = refseq utilities:\n";
   	$USAGE    .= "\n\t    1 = split a GLUE file with multiple sequences into individual files"; 
   	$USAGE    .= "\n\t    2 = convert GLUE refseq to FASTA ORFs as nucs";
   	$USAGE    .= "\n\t    3 = convert GLUE refseq to FASTA ORFs as aminos";
   	$USAGE    .= "\n\t    4 = create publication style formatted reference sequence";
 	$USAGE  .= "\n\n\t  usage: $0 [options] -i [infile]\n\n";
	
	# Create a sequence object
	my $seq_obj = Sequence->new();
	my $parser_obj = RefSeqParser->new();
	if ($mode eq 1) { # Split a VGLUE reference file
		unless ($file) { die $USAGE; }
		print "\n\t # Splitting refseq file '$file'";
		my @refseqs;
		$parser_obj->split_vglue_ref_file($file, \@refseqs);
		foreach my $refseq_ref (@refseqs) {
			# Extract and parse the features block
			my %refseq_data;
			$parser_obj->parse_refseq_metadata($refseq_ref, \%refseq_data);
			my $name = $refseq_data{name};
			unless ($name) { die "No name found for refseq"; }
			#$devtools->print_hash(\%refseq_data); die;
			# Write the file
			if ($outpath) {  $name = $outpath . '/' . $name; }
			$fileio->write_output_file($name, $refseq_ref);
		}
	}
	elsif ($mode eq 2 or $mode eq 3) { #  
		# Parse the refseq file
		unless ($file) { die $USAGE; }
		print "\n\t # Getting FASTA ORFs from refseq '$file'";
		my %data;
		$parser_obj->parse_refseq_flatfile($file, \%data);
		
		# Get the ORFs from the refseq
		my $refseq = RefSeq->new(\%data);
		my $name = $refseq->{name};
		unless ($name) { die "No name found for refseq"; }
		my %orfs;
		$refseq->get_orfs(\%orfs);
		my @orf_names = keys %orfs;
		foreach my $orf_name (@orf_names) {
			my $seq = $orfs{$orf_name};
			# Write the file
			my $file_name = $name . '_' . $orf_name;
			if ($outpath) { 
				$file_name = $outpath . '/' . $name . '_' . $orf_name . '.fas';
			}
			if ($mode eq 2) {
				my $fasta = ">$orf_name\n$seq\n";
				$fileio->write_text_to_file($file_name, $fasta);
			}
			elsif ($mode eq 3) {
				my $aa_seq = $seq_obj->translate($seq);
				my $fasta = ">$orf_name\n$aa_seq\n";
				$fileio->write_text_to_file($file_name, $fasta);
			}
		}
	}
	elsif ($mode eq 4) { 
		
		unless ($file) { die $USAGE; }
		print "\n\t # Writing formatted refseq for '$file'";
		
		# Create annotated sequence as a single string
		my %data;
		$parser_obj->parse_refseq_flatfile($file, \%data);
		my $refseq = RefSeq->new(\%data);
		my %linear;
		$refseq->create_linear_formatted(\%linear);
		
		# Write the file
		my $name = $refseq->{name};
		unless ($name) { die "No name found for refseq"; }
		my $file_name = $name;
		if ($outpath) { 
			$file_name = $outpath . '/' . $name . '.text';
		}
		$refseq->write_linear_formatted_seq(\%linear, $file_name);
	}
	elsif ($mode eq 5) { 
		my %typical;	
		my $threshold = 0.1; # Temp hard setting
		my $alignment = RefSeqAlignment->new();
		$alignment->derive_typical_list($file, \%typical, $threshold);
		$alignment->write_typical_list('typical_list.txt', \%typical);
	}
	else {
		die $USAGE;
	}	
}

############################################################################
# SECTION: Define reference sequence core (online version) 
############################################################################

#***************************************************************************
# Subroutine:  define_refseq
# Description: define a reference sequence
#***************************************************************************
sub define_refseq {

	my ($self, $query) = @_;

	# Get params
	my $stage = 1;
	my %refseq_data;
	$self->get_refseq_form_params($query, \%refseq_data);

	# Create the reference sequence and write it to the report directory
	my $process_id  = $self->{process_id};
	my $output_path = $self->{output_path};
	my $report_dir  = $output_path . $process_id . '/';
	my $mkdir_cmd   = "mkdir $report_dir";
	system $mkdir_cmd;
	
	my $sequence     = $refseq_data{sequence};
	my $metadata_ref = $refseq_data{metadata};
	my $refseq_name  = $metadata_ref->{name};
	unless ($refseq_name) { 
		$devtools->print_hash($metadata_ref); die; 
	}
	my $fasta  = ">$refseq_name\n$sequence\n\n";
	my $raw_path = $report_dir . 'raw.fa';
	$fileio->write_text_to_file($raw_path, $fasta);
	# Convert mac line breaks
	my $command = "perl -pi -e 's/\r/\n/g' $raw_path";
	system $command;
	my @fasta;
	$seqio->read_fasta($raw_path, \@fasta);
	my $seq_ref = shift @fasta;
	$sequence = $seq_ref->{sequence};
	
	# Get the params from the form
	my $seq_type     = $refseq_data{sequence_type};
	my $eve_type     = $refseq_data{eve_type};
	my $genome_type  = $refseq_data{genome_type};
	my $genes_ref    = $refseq_data{genes};
	$refseq_data{sequence} = $sequence;
	my $refseq = RefSeq->new(\%refseq_data);
	$refseq->write_self_to_text($report_dir);
	
	# Create the main block (left panel) for the refseq HTML
	my $views_obj = Views->new($self);
	# First create the form for adding further annotation
	my @refseq_main;
	my %data;
	$data{refseq_name} = $refseq->{name};
	$data{report_dir}  = $report_dir;
	$self->create_form_part(\@refseq_main, \%data);

	# Then add refseq HTML below
	$views_obj->create_refseq_view($refseq, \@refseq_main);
	#$devtools->print_array(\@refseq_main);
	# Write the main block (left panel) to the report directory
	my $page_path = $report_dir . $refseq_name . '.html'; 
	#print "\n\t page path $page_path";
	$fileio->write_output_file($page_path, \@refseq_main);
    # Define the page in the correct format and write it
	my %site;
	$site{report_dir} = $report_dir;
	my %page;
	$page{title} = "Paleovirology online: defining GLUE reference sequence for $refseq_name";
	$views_obj->define_paleo_page(\%site, \%page, $refseq_name,  3, 'tool');
	my $pages_ref = $site{pages_hash};
	my $page_ref = $pages_ref->{$refseq_name};
	my @html;
	#$writer->assemble_page($page_ref, \@html);
	#my $output = join("\n", @html);
	#print $output;
}

#***************************************************************************
# Subroutine:  append_to_refseq
# Description: 
#***************************************************************************
sub append_to_refseq {

	my ($self, $query) = @_;

	# Get params
	my $stage = 1;
	my %form_data;
	$self->get_refseq_append_form_params($query, \%form_data);
	my $metadata_ref = $form_data{metadata};

	# Get the params from the form
	my $feature_type = $form_data{feature_type};
	my $report_dir   = $form_data{report_dir};
	my $append_path  = $form_data{append_path};
	unless ($append_path) { die "NO PATHS FOR APPEND"; }
	#$devtools->print_web_hash($refseq); die;

	# Parse the reference sequence and add new data
	my @refseq;
	my $refseq_parse_obj = RefSeqParser->new();
	my %params; 
	my $parser = RefSeqParser->new();
	$parser->parse_refseq_flatfile($append_path, \%params);
	my $refseq = RefSeq->new(\%params);
	my $refseq_name  = $refseq->{name};
	unless ($refseq_name) { die; }
	$refseq->append_new_data(\%form_data);
	$refseq->write_self_to_text($report_dir);
	#$devtools->print_web_hash($refseq);
	
	# Create the main block (left panel) for the refseq HTML
	my $views_obj = Views->new();

	# First create the form for adding further annotation
	my @refseq_main;
	my %data;
	$data{refseq_name} = $refseq->{name};
	$data{report_dir}  = $report_dir;
	$self->create_form_part(\@refseq_main, \%data);

	# Then add refseq HTML below
	$views_obj->create_refseq_view($refseq, \@refseq_main);
	# Write the main block (left panel) to the report directory
	my $page_path = $report_dir . $refseq_name . '.html'; 
	$fileio->write_output_file($page_path, \@refseq_main);
    # Define the page in the correct format and write it
	my @html;
	my %site;
	$site{report_dir} = $report_dir;
	my %page;
	$page{title} = "Paleovirology online: defining GLUE reference sequence for $refseq_name";
	$views_obj->define_paleo_page(\%site, \%page, $refseq_name,  3, 'tool');
	my $pages_ref = $site{pages_hash};
	my $page_ref = $pages_ref->{$refseq_name};
	#$writer->assemble_page($page_ref, \@html);

	# Write the main block (left panel) to the refseq HTML directory
	$fileio->write_output_file($page_path, \@html);
	my $output = join("\n", @html);
	print $output;
}

#***************************************************************************
# Subroutine: create_form_part 
# Description: 
#***************************************************************************
sub create_form_part {
	
	my ($self, $output_ref, $hash_ref) = @_;

	my $refseq_name = $hash_ref->{refseq_name};
	my $report_dir  = $hash_ref->{report_dir};
	my $append_path = $report_dir . $refseq_name;

	my $html = '';
	$html .= '<h4>Annotate an existing EVE reference sequence</h4>';
	$html .= '<div class="separator"></div>';
	$html .= '<form method="post" action="http://saturn.adarc.org/vglue_sandbox/glue.cgi" 
				enctype="multipart/form-data" id="resform">
				<p>';
	
	# Reference sequence top line
	$html .= "\n\t\t<INPUT TYPE='hidden' NAME='mode' VALUE='REFSEQ_APPEND'>";
	$html .= "\n\t\t<INPUT TYPE='hidden' NAME='append_path' VALUE='$append_path'>";
	$html .= "\n\t\t<INPUT TYPE='hidden' NAME='report_dir'  VALUE='$report_dir'>";
	$html .= "\n\t\t<INPUT TYPE='hidden' NAME='family' VALUE='Retroviridae'>";
	
	#  Sequence feature
	$html .= '<br><b>Define a new sequence feature</b><br><br>';
	$html .= "Feature name &nbsp; <input type='text' name='feature_name' value='Pol' ";
	$html .= "size='4' maxlength='30' />";
	$html .= '&nbsp;&nbsp;';
	$html .= '<select name="feature_type">';
	$html .= "<option selected value='ORF'> ORF</option>";
	$html .= "<option value='UTR'> UTR</option>";
	$html .= '</select>&nbsp; &nbsp;';
	$html .= "Coordinates &nbsp; <input type='text' name='csv_coordinates' value='4744,7095'";
	$html .= ' size="18" maxleform/ngth="80"/><BR><BR>';	
	$html .= "\n\t\t<div><input type='hidden' name='action' value='ANALYZE'>";
	$html .= "\n\t\t<a href='javascript:document.forms[0].submit();'>";
	$html .= "\n\t\t<img src='http://saturn.adarc.org/vglue_sandbox/site/images/form/btn.png' alt='ANALYZE'/></a>";
	$html .= "\n\t<br></div></form><br>";
	$html .= '<div class="separator"></div>';	
	push (@$output_ref, $html);	

	my $link_text  = "<br> Click <u><a href='$append_path'>here</a></u>";
	$link_text    .= " to view the raw reference sequence <br><br>";
	push (@$output_ref, $link_text);	

}

############################################################################
# SECTION: Functions for getting values from web forms
############################################################################

#***************************************************************************
# Subroutine:  get_refseq_form_params 
# Description: get forms from the start page of the 'refseq define' process
#***************************************************************************
sub get_refseq_form_params {

	my ($self, $query, $data_ref) = @_;

	# Create writer utility obj
	my $writer = HTML_Utilities->new();	
	
	# Get pasted sequence data if its there
	my $referrer = 'vglue_configure';
	my %exons;
	my %feature;
	my %metadata;
	my $raw_data = $query->param('SequenceData');
	if ($raw_data) {
		
		# Create sequence
		my @dataset = split("\n", $raw_data); 
		my @f_dataset;
		foreach my $line (@dataset) { push (@f_dataset, "$line\n"); }
		my $sequence = join ('', @f_dataset);
		$sequence =~ s/\n//g;

		# Initialise CGI output
		my $append_path  = $query->param('append_path');
		my $report_dir   = $query->param('report_dir');
		$data_ref->{report_dir}    = $report_dir;
		$data_ref->{append_path}   = $append_path;
		$data_ref->{sequence}      = $sequence;
		$data_ref->{raw_data}      = \@f_dataset;
		
		# Metadata 
		my $supertribe = $query->param('tribe');
		my $tribe;
		if ($supertribe eq 'Alnilamretrovirinae') {
			$tribe = 'Alnilam';
		}
		elsif ($supertribe eq 'Alnitakretrovirinae') {
			$tribe = 'Alnitak';
		}
		elsif ($supertribe eq 'Mintakaretrovirinae') {
			$tribe = 'Mintaka';
		}
		else { die; }

		$metadata{name}             = $query->param('name');
		$metadata{full_name}        = $query->param('full_name');
		$metadata{family}           = $query->param('family');
		$metadata{virus_subfamily}  = $query->param('subfamily');
		$metadata{virus_supertribe} = $supertribe;
		$metadata{virus_tribe}      = $tribe;
		$metadata{virus_genus}      = $query->param('genus');
		$metadata{virus_subgroup}   = $query->param('subgroup');
		$metadata{host_sci_name}    = $query->param('host_sci_name');
		$metadata{host_common_name} = $query->param('host_common_name');
		$metadata{accession}        = $query->param('accession');
		$metadata{sequence_type}    = $query->param('sequence_type');
		$metadata{genome_type}      = $query->param('genome_type');
		$metadata{genome_coverage}  = $query->param('genome_coverage');
		
		# Features 
		my $feature_type = $query->param('feature_type');
		my $feature_name = $query->param('feature_name');
		my $coordinates  = $query->param('csv_coordinates');
		if ($feature_type and $feature_name and $coordinates) {
			
			my %starts;
			my @starts;
			my @coordinates = split(',', $coordinates);
			my $odd = undef;
			my $start = undef;
			my $coding_start;  # THIS MIGHT BE POINTLESS - need to consolidate
			my $first_start = undef;
			my $stop = undef;
			my $num_exons = 0;
			foreach my $position (@coordinates) {
				
				#print "<BR> $position";
				if ($odd) { 
					$stop = $position;
					$starts{$start} = $position;
					$start = undef;
					$num_exons++;
					$odd = undef;  
				}
				else { 
					$start = $position;
					push (@starts, $start);
					unless ($first_start) {
						$first_start = $start;
					}
					$odd = 'true'; 
					
				}	
			}
			if ($odd) { die; } # Not an even number of coordinates
			
			$feature{stop}         = $stop;
			$feature{coding_stop}  = $stop;
			$feature{starts}       = \@starts;
			$feature{start}        = $first_start;
			$feature{coding_start} = $first_start;
			$feature{type}         = $feature_type;
			$feature{name}         = $feature_name;
			$feature{full_name}    = $feature_name;
			$feature{num_exons}    = $num_exons;
			$feature{exons}        = \%starts;
			#$feature{csv_coordinates}  = $coordinates;
		}
	}
	# print the formatted input page if no sequence data has been received
	else {
		# Read  the VGLUE header
		my @html;
		$fileio->read_input_file('./site/html/refseq_define.html', \@html);
		my $html  = join('', @html);
		print $html;
		exit;
	}

	# Make the input structured the same as it would be from the refseq parser
	my @features;
	push (@features, \%feature);
	$data_ref->{genes}    = \@features;	
	$data_ref->{metadata} = \%metadata;	
}

#***************************************************************************
# Subroutine: get_refseq_append_form_params
# Description: 
#***************************************************************************
sub get_refseq_append_form_params {

	my ($self, $query, $data_ref) = @_;

	# Create writer utility obj
	my $writer = HTML_Utilities->new();	
	
	# Get pasted sequence data if its there
	my %exons;
	my %feature;
	my %metadata;
	
	# Initialise CGI output
	my $append_path  = $query->param('append_path');
	my $report_dir   = $query->param('report_dir');
	my $genome_type  = $query->param('genome_type');
	my $seq_type     = $query->param('sequence_type');
	$data_ref->{report_dir}    = $report_dir;
	$data_ref->{append_path}   = $append_path;
	
	# Features 
	my $feature_type = $query->param('feature_type');
	my $feature_name = $query->param('feature_name');
	my $coordinates  = $query->param('csv_coordinates');
	if ($feature_type and $feature_name and $coordinates) {
		
		my %starts;
		my @starts;
		my @coordinates = split(',', $coordinates);
		my $odd = undef;
		my $start = undef;
		my $coding_start;  # THIS MIGHT BE POINTLESS - need to consolidate
		my $first_start = undef;
		my $stop = undef;
		my $num_exons = 0;
		foreach my $position (@coordinates) {
			
			#print "<BR> $position";
			if ($odd) { 
				$stop = $position;
				$starts{$start} = $position;
				$start = undef;
				$num_exons++;
				$odd = undef;  
			}
			else { 
				$start = $position;
				push (@starts, $start);
				unless ($first_start) {
					$first_start = $start;
				}
				$odd = 'true'; 
				
			}	
		}
		if ($odd) { die; } # Not an even number of coordinates
		
		$feature{stop}         = $stop;
		$feature{coding_stop}  = $stop;
		$feature{starts}       = \@starts;
		$feature{start}        = $first_start;
		$feature{coding_start} = $first_start;
		$feature{type}         = $feature_type;
		$feature{name}         = $feature_name;
		$feature{full_name}    = $feature_name;
		$feature{num_exons}    = $num_exons;
		$feature{exons}        = \%starts;
		#$feature{csv_coordinates}  = $coordinates;
	}
	else { # print the formatted input page if no sequence data has been received
		# TO DO FIX THIS
		die;
		# Read  the VGLUE header
		my $views_obj = Views->new();
		my $site_ref = $views_obj->{paleo_site};
		my @page;
		#$writer->assemble_page($site_ref, $referrer, \@page);
		print @page;
		exit;
	}

	# Make the input structured the same as it would be from the refseq parser
	my @features;
	push (@features, \%feature);
	$data_ref->{genes}    = \@features;	
	$data_ref->{metadata} = \%metadata;	
}

############################################################################
# Command line title blurb 
############################################################################

#***************************************************************************
# Subroutine:  show_title
# Description: does what it says 
#***************************************************************************
sub show_title {

	$console->refresh();
	my $title       = 'GLUE Data tool';
	my $version     = '2.0';
	my $description = 'GLUE-associated data tool';
	my $author      = 'Robert J. Gifford';
	my $contact		= '<robert.gifford@glasgow.ac.uk>';
	$console->show_about_box($title, $version, $description, $author, $contact);
}

#***************************************************************************
# Subroutine:  by_number
# Description: by number - for use with perl 'sort'  (cryptic but works) 
#***************************************************************************
sub by_number { $a <=> $b }	

############################################################################
# EOF
############################################################################
