#!/usr/bin/perl -w
############################################################################
# Module:      DataToolCmd.pm
# Description: Command line interface to the GLUE framework data tool
# History:     November 2014: Created by Robert Gifford 
############################################################################
package DataToolCmd;

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
use Base::Sequence; 

############################################################################
# Globals
############################################################################

# Base class instantiations
my $fileio     = FileIO->new();
my $seqio      = SeqIO->new();
my $bioio      = BioIO->new();
my $devtools   = DevTools->new();
my $console    = Console->new();
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
		datatool_obj         =>	$parameter_ref->{datatool_obj},
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
# Subroutine:  show_help_page
# Description: show help page information
#***************************************************************************
sub show_help_page {

	my($HELP) = "\n\t usage: $0 -m[options] -s[options] -f[files] -d[files]
      \n\t # Convert between formats";
	  $HELP  .= "\n\t  -m=1        :   FASTA to Delimited";
	  $HELP  .= "\n\t  -m=2        :   Delimited to FASTA";
	  $HELP  .= "\n\t  -m=3        :   Genbank to FASTA+DATA";
	  $HELP  .= "\n\t  -m=4        :   Genbank to REFSEQ"; 
	  $HELP  .= "\n\t  -m=5        :   REFSEQ to linear formatted"; 
	  $HELP  .= "\n\t  -m=6        :   REFSEQ to FASTA";
	  $HELP  .= "\n\t  -m=7        :   Concatenate FASTA files";
	  $HELP  .= "\n\t  -m=x        :   FASTA to NEXUS";
	  $HELP  .= "\n\t  -m=x        :   FASTA to PHYLIP";
	  $HELP  .= "\n\t  -m=x        :   Expand GLUE MSA";
	  $HELP  .= "\n\t  -m=x        :   Compress GLUE MSA";
	  $HELP  .= "\n\n\t # Sorting, filtering, linking";
	  $HELP  .= "\n\t  -s=1        :   Shorten FASTA headers"; 
	  $HELP  .= "\n\t  -s=2        :   Sort sequences by length"; 
	  $HELP  .= "\n\t  -s=3        :   Filter sequences by keyword"; 
	  $HELP  .= "\n\t  -s=4        :   Link two data files using a shared ID";
	  $HELP  .= "\n\t  -s=2        :   Sort sequences by data column"; 
	print $HELP;

}

#***************************************************************************
# Subroutine:  run_reformat_tools_cmd_line
# Description: hand off to FASTA utilities
#***************************************************************************
sub run_reformat_tools_cmd_line {

	my ($self, $infiles, $mode) = @_;

	$self->show_title();

	shift @$infiles; # Not sure why there's an empty array element
	#$devtools->print_array($infiles); die;

	my $result;
	if ($mode eq 1) {  # FASTA to delimited txt table 
		$result = $self->convert_fasta_to_delimited($infiles);
	}
	elsif ($mode eq 2) { # Delimited to FASTA
		$result = $self->convert_delimited_to_fasta($infiles);
	}
	elsif ($mode eq 3) {  # Genbank to FASTA + Data
		$self->convert_genbank_to_fasta_and_tabdelim($infiles);
	}
	elsif ($mode eq 4) {  # GenBank to REFSEQ
		$self->convert_genbank_file_to_GLUE_refseq($infiles);
	}
	elsif ($mode eq 5) {
		$self->write_GLUE_refseq_in_linear_format($infiles);
	}
	elsif ($mode eq 6) {
		$self->write_GLUE_refseq_as_fasta($infiles);
	}
	elsif ($mode eq 7) {  # FASTA to NEXUS
		$self->concatenate_fasta($infiles);
	}
	elsif ($mode eq 8) {  # FASTA to NEXUS
		$self->convert_fasta_file_to_nexus($infiles);
	}
	elsif ($mode eq 9) {  # FASTA TO PHYLIP
		$self->convert_fasta_file_to_phylip($infiles);
	}
	elsif ($mode eq 10) {  # SPLIT multi-sequence GLUE file 
	
	}
	elsif ($mode eq 11) {  # 
		$self->create_mutation_list_from_glue_msa($infiles);	
	}
	elsif ($mode eq 13) {
		die;
		$seqio->expand_glue_msa($infiles);
	}
	elsif ($mode eq 14) {
		die;
		$seqio->compress_glue_msa($infiles);
	}
	else { die; }

	return $result;
}

#***************************************************************************
# Subroutine:  run_sort_tools_cmd_line
# Description: hand off to sorting utilities
#***************************************************************************
sub run_sort_tools_cmd_line {

	my ($self, $seqfiles, $datafiles, $sort) = @_;
	
	my $process_id  = $self->{process_id};
	my $output_path = $self->{output_path};
	unless ($output_path) { $output_path .= './'; }		
	
	$self->show_title();
	shift @$seqfiles; # Not sure why there's an empty array element
	
	if ($sort eq 1) {  # Shorten FASTA headers (capture data)

	}
	elsif ($sort eq 2) {  # Data + sequence tools
		$self->do_length_based_fasta_sort($seqfiles);
	}
	elsif ($sort eq 3) {
		$self->do_header_based_fasta_filter($seqfiles)
	}
	elsif ($sort eq 4) {
		my $builder = GLUE_Build->new();
		my $report_dir = $output_path . $process_id . '/';
		$fileio->create_unique_directory($report_dir, 'text');
		$self->{report_dir} = $report_dir;		
	}
	elsif ($sort eq 5) {  
		$self->link_sequences_and_data();
	}
	elsif ($sort eq 6) {

	}
	elsif ($sort eq 7) {
	
	}
	else { die; }
}

#***************************************************************************
# Subroutine:  convert_fasta_to_delimited
# Description: convert FASTA formatted filed to tab-deimited columns 
#              column 1 (ID), column 2 (sequence)
#***************************************************************************
sub convert_fasta_to_delimited {

	my ($self, $infiles) = @_;

	my $datatool= $self->{datatool_obj};

	my $result = undef;
	foreach my $infile (@$infiles) {

		print "\n\t # Converting file '$infile' from FASTA to data";
		my @fasta;
		$seqio->read_fasta($infile, \@fasta);
		my $num_seqs = scalar @fasta;
		unless ($num_seqs) { 
			print "\n\t # No sequences read from file '$infile'";
			next;
		}
		my $delimiter;
		my $question = "\n\t What is the delimiter (comma\/tab)";
		my @choices  = qw [ c t ];
		my $choice = $console->ask_simple_choice_question($question, \@choices);
		if ($choice eq 'c') { $delimiter = ',';  }
		if ($choice eq 't') { $delimiter = "\t"; }
		my @data;
		$datatool->fasta_to_delimited(\@fasta, \@data, $delimiter);
		my $outfile = $infile . '.tsv.txt';
		$result = $fileio->write_file($outfile, \@data);
		print "\n\t # File '$outfile' created";

	}
	return $result;
}

#***************************************************************************
# Subroutine:  convert_delimited_to_fasta
# Description: convert tab-delimited columns (seq in last column) to FASTA
#***************************************************************************
sub convert_delimited_to_fasta {

	my ($self, $infiles) = @_;

	my $datatool= $self->{datatool_obj};

	foreach my $infile (@$infiles) {
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
		$datatool->delimited_to_fasta(\@data, \@fasta, $delimiter);
		$fileio->write_file($outfile, \@fasta);
		print "\n\t # File '$outfile' created";
		#$devtools->print_array(\@data); die;
	}
}

#***************************************************************************
# Subroutine:  convert_genbank_to_fasta_and_tabdelim
# Description: convert genbank to (i) fasta and (ii) data in a tab-delimited file
# Notes:       can deal with multiple, concatenated Genbank files
#***************************************************************************
sub convert_genbank_to_fasta_and_tabdelim {

	my ($self, $infiles) = @_;

	my $datatool= $self->{datatool_obj};

	foreach my $infile (@$infiles) {
		print "\n\t # Converting file '$infile' from GenBank format to FASTA+DATA";
		my @fasta;
		my @data;
		$datatool->genbank_to_fasta_and_data($infile, \@fasta, \@data);
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
}

#***************************************************************************
# Subroutine:  convert_genbank_file_to_GLUE_refseq
# Description: convert genbank file to a GLUE refseq with features table
#***************************************************************************
sub convert_genbank_file_to_GLUE_refseq {

	my ($self, $infiles) = @_;

	my $datatool= $self->{datatool_obj};

	foreach my $infile (@$infiles) {
		print "\n\t # Converting file '$infile' from GenBank format to GLUE format";
		my @genbank;
		my @refseq;
		$bioio->parse_gb_to_refseq($infile, \@refseq);
		my $refseq_path = $infile . '.glu';
		$fileio->write_file($refseq_path, \@refseq);
		print "\n\t # File $refseq_path created\n\n";
	}
}

#***************************************************************************
# Subroutine:  write_GLUE_refseq_in_linear_format
# Description: write GLUE reference sequence in linear format with both
#              nucleotide sequence and ORF features displayed 
#***************************************************************************
sub write_GLUE_refseq_in_linear_format {

	my ($self, $infiles) = @_;

	my $datatool= $self->{datatool_obj};

	foreach my $infile (@$infiles) {
		print "\n\t # Writing GLUE reference sequence '$infile' in linear format ";
		my $parser = RefSeqParser->new();
		my %params;
		$parser->parse_refseq_flatfile($infile, \%params);
		my $refseq = RefSeq->new(\%params);
		my %linear;
		$refseq->create_linear_formatted(\%linear);
		my $outfile = $infile . '.formatted.txt';
		$refseq->write_linear_formatted_seq(\%linear, $outfile);
		print "\n\t file '$outfile' created\n\n";
	}
}

#***************************************************************************
# Subroutine:  write_GLUE_refseq_as_fasta
# Description: write GLUE reference sequence as FASTA, with separate files
#              for ORFs (NT), translated ORFs (AA) and UTRs
#***************************************************************************
sub write_GLUE_refseq_as_fasta {

	my ($self, $infiles) = @_;

	my $datatool= $self->{datatool_obj};

	foreach my $infile (@$infiles) {
	
		my $seq_obj = Sequence->new();    # Create a sequence object
		my $parser = RefSeqParser->new();
		my %params;
		$parser->parse_refseq_flatfile($infile, \%params);
		my $refseq = RefSeq->new(\%params);
		my $refseq_name = $refseq->{'name'};
		my %orfs;
		$refseq->get_orfs(\%orfs);
		my @orfs = sort keys %orfs;
		my $outfile = $infile . '.orfs.fa';
		my @outfile1;
		my @outfile2;
		foreach my $orf (@orfs) {
			my $sequence = $orfs{$orf};
			my $fasta = ">$refseq_name $orf\n$sequence\n";
			push (@outfile1, $fasta);
			my $aa_seq = $seq_obj->translate($sequence);
			my $aa_fasta = ">$refseq_name $orf\n$aa_seq\n";
			push (@outfile2, $aa_fasta);
		}
		my $outfile1 = $refseq_name . '.orfs.nt.fas';
		$fileio->write_file($outfile1, \@outfile1);
		my $outfile2 = $refseq_name . '.orfs.aa.fas';
		$fileio->write_file($outfile2, \@outfile2);
		#print "\n\t file '$outfile' created\n\n";
	}
}

#***************************************************************************
# Subroutine:  convert_fasta_file_to_nexus
# Description: command line interface for a FASTA-to-NEXUS conversion
#***************************************************************************
sub convert_fasta_file_to_nexus {

	my ($self, $infiles) = @_;

	my $datatool= $self->{datatool_obj};

	foreach my $infile (@$infiles) {
		
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
}

#***************************************************************************
# Subroutine:  convert_fasta_file_to_phylip 
# Description: command line interface for a FASTA-to-PHYLIP conversion
#***************************************************************************
sub convert_fasta_file_to_phylip {

	my ($self, $infiles) = @_;

	my $datatool= $self->{datatool_obj};

	foreach my $infile (@$infiles) {
			
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
# Subroutine:  concatenate_fasta
# Description: concatenate FASTA files
#***************************************************************************
sub concatenate_fasta {

	my ($self, $infiles) = @_;

	my $datatool= $self->{datatool_obj};

	my %lengths;
	my %concatenated;
	my %seq_ids;
	foreach my $infile (@$infiles) {

		print "\n\t infile '$infile'";
		my $seq_len;
		my %fileseqs;
		
		# Sequences
		my @sequences;
		$seqio->read_fasta($infile, \@sequences);
		foreach my $seq_ref (@sequences) {
			my $seq_id   = $seq_ref->{header};
			my $sequence = $seq_ref->{sequence};
			$seq_ids{$seq_id} = 1;
			$seq_len  = length $sequence;			
			if ($fileseqs{$seq_id}) {
				die;
			}
			else {
				$fileseqs{$seq_id} = $sequence;
			}
		}
		$concatenated{$infile} = \%fileseqs;
		$lengths{$infile}      = $seq_len;
	}
	#$devtools->print_hash(\%concatenated); die;

	
	my @segments;
	my @concatenated;
	my @seqids = keys %seq_ids;
	foreach my $seq_id (@seqids){

		my $segments = "";
		my $concatenated_seq = "";
		foreach my $infile (@$infiles) {
			
			my $concatenation_ref = $concatenated{$infile};
			my $sequence = $concatenation_ref->{$seq_id};
			unless ($sequence) {
				print "\n\t No sequence for '$seq_id' in file '$infile";
				my $seq_len = $lengths{$infile};
				my $pad = '-' x $seq_len;
				$concatenated_seq .= $pad;	
				$segments .= '-';	
			}
			else {
				$concatenated_seq .= $sequence;	
				$segments .= 'X';	
			}
		}
		my $fasta = ">$seq_id\n$concatenated_seq\n\n";
		push (@concatenated, $fasta);
		my $line = "$seq_id\t$segments\n";
		push (@segments, $line);
	}
	
	# Write files	
	$fileio->write_file('concatenated.fas', \@concatenated);
	$fileio->write_file('segment_matching.txt', \@segments);

}

#***************************************************************************
# Subroutine:  split_multisequence_glue_file
# Description: FASTA to delimited 
#***************************************************************************
sub split_multisequence_glue_file {

	my ($self, $file) = @_;

	print "\n\t # Splitting refseq file '$file'";
	my @refseqs;
	my $parser_obj = RefSeqParser->new();
	$parser_obj->split_glue_refseq_file($file, \@refseqs);

	foreach my $refseq_ref (@refseqs) {
		
		# Extract and parse the features block
		my %refseq_data;
		$parser_obj->parse_refseq_metadata($refseq_ref, \%refseq_data);
		my $name = $refseq_data{name};
		unless ($name) { die "No name found for refseq"; }
		#$devtools->print_hash(\%refseq_data); die;
		
		# Write the file
		$fileio->write_output_file($name, $refseq_ref);
	}
}

#***************************************************************************
# Subroutine:  create_mutation_list_from_glue_msa
# Description: cmd line interface to derive a mutation list from a GLUE MSA
#***************************************************************************
sub create_mutation_list_from_glue_msa {

	my ($self, $file) = @_;

	# Get the settings for the list
	my %list;	
	my $question2 = "\n\t Record mutations above or below threshold?";
	my @choices = qw [ above below ];
	my $rule  = $console->ask_simple_choice_question($question2, \@choices);
	my $question1 = "\n\t Specify threshold proportion as a percentage";
	my $threshold = $console->ask_int_with_bounds_question($question1, 0, 100);
	
	# Create the list
	my $alignment = RefSeqAlignment->new();
	$alignment->derive_typical_list($file, \%list, $threshold);
	$alignment->write_typical_list('typical_list.txt', \%list);
}

#***************************************************************************
# Subroutine:  convert_fasta_headers
# Description: write GLUE reference sequence in linear format with both
#***************************************************************************
sub convert_fasta_headers {

	my ($self, $seqfile) = @_;

	die;
	my @sequences;
	$seqio->read_fasta($seqfile, \@sequences);
	my @truncated;
	foreach my $seq_ref (@sequences) {
		$seq_ref->truncate_header();
		my $header = $seq_ref->{header};
		my $sequence = $seq_ref->{sequence};
		my  $fasta = ">$header\n$sequence\n";
		push (@truncated, $fasta);
	}
	my $outfile = $seqfile .= '.truncated.fas';
	$fileio->write_file($outfile, \@truncated);
	print "\n\t ### File '$outfile' created";
}

#***************************************************************************
# Subroutine:  do_length_based_fasta_sort
# Description: command line interface for a lenth based FASTA sort 
#***************************************************************************
sub do_length_based_fasta_sort {

	my ($self, $seqfiles) = @_;

	my $datatool = $self->{datatool_obj};

	foreach my $seqfile (@$seqfiles) {

		my @fasta;
		$seqio->read_fasta($seqfile, \@fasta);
		my $high;
		my $count = 0;
		foreach my $seq_ref (@fasta) {
			$count++;
			my $sequence   = $seq_ref->{sequence};
			my $seq_length = length $sequence;
			unless ($high) { $high = $seq_length; }
			if ($seq_length > $high) { $high = $seq_length; } 
		}
		unless ($high) { die; }
		my $question1 = "\n\t Specify sequence length minimum to use as a cutoff";
		my $minimum = $console->ask_int_with_bounds_question($question1, 1, $high);
		my $question2 = "\n\t Count gaps?";
		my $count_gaps = $console->ask_yes_no_question($question2);
		if ($count_gaps eq 'n') { $count_gaps = undef; }
		my @sorted;
		my $excluded = 0;
		$excluded = $datatool->sort_seqs_by_length(\@fasta, \@sorted, $minimum, $count_gaps);
		my $outfile = $seqfile .= '.sorted.fas';
		$fileio->write_file($outfile, \@sorted);
		print "\n\t ### File '$outfile' created";
		print "\n\t ### There were a total of $count sequences";
		my $included = $count;
		if ($excluded) {
			$included =  $count-$excluded;
			print "\n\t ### $excluded sequences shorter than $minimum nt were excluded";
		}
		print "\n\t ### New file containes $included sequences";
	}
}

#***************************************************************************
# Subroutine:  do_header_based_fasta_filter
# Description: command line interface for a FASTA header-based filter fxn
#***************************************************************************
sub do_header_based_fasta_filter {

	my ($self, $infile) = @_;

	my $sequences_ref = $self->{sequences};
	my $question1 = "\n\t Enter phrase to filter on:";
	my $word  = $console->ask_question($question1);
	my $question2 = "\n\t Exclude or include?";
	my @choices = qw [ exclude include ];
	my $rule  = $console->ask_simple_choice_question($question2, \@choices);
	my $question3 = "\n\t Ignore case?";
	my $case  = $console->ask_yes_no_question($question3);
	my @filtered;
	my $filtered = 0;
	$filtered = $self->filter_by_header($sequences_ref, \@filtered, $word, $rule, $case);
	print "\n\t filtered $filtered sequences with headers containing phrase '$word'";
	# convert to sequence objects
	my @filtered_seqs; 	
	$seqio->convert_fasta(\@filtered, \@filtered_seqs);
	$self->{sequences} = \@filtered_seqs;
}

#***************************************************************************
# Subroutine:  do_fasta_sort_using_delimited_file_column 
# Description: sort a FASTA file using a column of values in a delimited file
#***************************************************************************
sub do_fasta_sort_using_delimited_file_column {

	my ($self, $seqfile, $datafile) = @_;

	my @fasta;
	$seqio->read_fasta($seqfile, \@fasta);
	my @data;
	$fileio->read_file($datafile, \@data);
	my $delimiter;
	my $question = "\n\t What is the delimiter (comma\/tab)";
	my @choices  = qw [ c t ];
	my $choice = $console->ask_simple_choice_question($question, \@choices);
	if ($choice eq 'c') { $delimiter = ',';  }
	if ($choice eq 't') { $delimiter = "\t"; }
	# get the indices
	$question = "\n\t Which column contains the sequence IDs?";
	my $id_index = $console->ask_int_question($question);
	$id_index--;
	$question = "\n\t Which column contains the data to sort on?";
	my $sub_index = $console->ask_int_question($question);
	$sub_index--;
	die "\n\t UNFINISHdED SORRY\n\n"; #TODO
	$self->sort_sequences_on_data_column(\@fasta, \@data, $id_index, $sub_index, $delimiter);

}

#***************************************************************************
# Subroutine:  combine_data_from_two_files 
# Description: 
#***************************************************************************
sub combine_data_from_two_files {

	my ($self, $seqfile, $datafile) = @_;

	my @data1;
	$fileio->read_file($seqfile, \@data1);
	my @data2;
	$fileio->read_file($datafile, \@data2);

	my $delimiter;
	my $question = "\n\t What is the delimiter (comma\/tab)";
	my @choices  = qw [ c t ];
	my $choice = $console->ask_simple_choice_question($question, \@choices);
	if ($choice eq 'c') { $delimiter = ',';  }
	if ($choice eq 't') { $delimiter = "\t"; }
	
	# get the indices
	$question = "\n\t Which column contains the sequence IDs in file 1?";
	my $a_index = $console->ask_int_question($question);
	$a_index--;
	$question = "\n\t Which column contains the sequence IDs in file 2?";
	my $b_index = $console->ask_int_question($question);
	$b_index--;
	$self->combine_data(\@data1, \@data2, $a_index, $b_index);
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

############################################################################
# EOF
############################################################################
