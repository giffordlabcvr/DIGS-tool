#!/usr/bin/perl -w
############################################################################
# Module:      Data_tool.pm
# Description: GLUE-associated FASTA tool
# History:     January 2012: Created by Robert Gifford 
############################################################################
package Data_tool;

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
#use Base::Text_Utilities;

############################################################################
# Globals
############################################################################

# Create base objects
my $fileio     = FileIO->new();
my $seqio      = SeqIO->new();
my $devtools   = DevTools->new();
my $console    = Console->new();
#my $textutils  = Text_Utilities->new();
#my $io         = IO->new();

my $alignments_dir = "db/alignments/";  # Alignments directory
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
		tmp_path             => $parameter_ref->{tmp_path}, 
		refseq_use_path      => $parameter_ref->{refseq_use_path}, 
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# GET WEB FORM DATA
############################################################################

#***************************************************************************
# Subroutine:  do_fasta_function 
# Description: handler fxn for handing off to specific FASTA fxns 
#***************************************************************************
sub do_fasta_function {

	my ($self) = @_;

	# Create the folder  (todo: use call to fxn)
	my $output_type = $self->{output_type};
	my $process_id  = $self->{process_id};
	my $output_path = $self->{output_path};
	my $report_dir  = $output_path . $process_id . '/';
	my $mkdir_cmd   = "mkdir $report_dir";
	system $mkdir_cmd;
	$self->{report_dir} = $report_dir;
	#my $message = " created report directory $report_dir";
	#$io->show_output_message($message, $output_type);

	my $output_ref = $self->{output};
	my $raw_path = $report_dir . 'raw.fas';
	my $raw_sequences_ref = $self->{sequences};
	$seqio->write_fasta($raw_path, $raw_sequences_ref);
	my $rawlink  = "<p>Click <a href='$raw_path'>here</a> to retrieve the raw ";
	$rawlink .= "sequences.<BR><BR>";
	push (@$output_ref, $rawlink);


	# check if we are converting only
	my $skip_fasta_steps = undef;
	if ($self->{convert_method}) {
		if ($self->{convert_method} eq 'delimited_to_fasta') {
			$skip_fasta_steps = 1;
		}
	}

	# Apply FASTA only procedures
	unless ($skip_fasta_steps) {
		
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
	}	
	
	#===# CONVERT 	
	if ($self->{convert_format}) {
		$self->convert_format();
	}
}

#***************************************************************************
# Subroutine:  convert_format 
# Description: handler fxn handing off to specific format conversion fxns 
#***************************************************************************
sub convert_format {

	my ($self) = @_;

	# Convert sequences
	my $report_dir = $self->{report_dir};
	
	# Get sequence objects in array
	my $sequences_ref = $self->{sequences};
	unless ($sequences_ref) { die; }
	my @converted;
	if ($self->{convert_format} and $self->{convert_method} eq 'fasta_to_delimited') {  
		$self->fasta_to_delimited($sequences_ref, \@converted);		
	}	
	elsif ($self->{convert_format} and $self->{convert_method} eq 'delimited_to_fasta') {  
		
		if ($self->{output_type} eq 'command_line') {	

			my @delimited;
			my $tabdelim_path = $self->{tabdelim_path};
			my $delimiter     = $self->{delimiter};
			unless ($delimiter)      { die; }
			unless ($tabdelim_path)  { die "No input file defined\n\n"; }
			$fileio->read_input_file($tabdelim_path, \@delimited);		
			$self->delimited_to_fasta(\@delimited, \@converted, $delimiter);	
		}
		else {
			my $delimited_ref = $self->{delimited};  
			my $delimiter     = $self->{delimiter};
			unless ($delimiter and $delimited_ref) { die; }
			$self->delimited_to_fasta($delimited_ref, \@converted, $delimiter);		
		}
	}
	
	# Write the output
	if ($self->{output_type} eq 'command_line') {	
		my $fasta_path = $self->{fasta_path};	
		my @path = split(/\//, $fasta_path);
		my $filename = pop @path;
		my $outfile = $report_dir .  $filename . '.converted.fa';
		print "\n\t ## Converted sequences written to $outfile";
		$fileio->write_output_file($outfile, \@converted);
	}
	else {
		my $output_ref = $self->{output};
		unless ($output_ref) { die; }
		my $outfile = 'converted_seqs.txt';
		my $outpath = $report_dir . $outfile;
		$fileio->write_output_file($outpath, \@converted);
		my $filelink  = "<p>Click <a href='$outpath'>here</a> to retrieve the";
		if ($self->{convert_method} ne 'delimited_to_fasta') {
			if ($self->{filtered}) { $filelink .= " filtered,"; }
			if ($self->{sorted})   { $filelink .= " sorted,";   }
		}
		$filelink .= " converted sequences <BR><BR>";
		push (@$output_ref, $filelink);
	}
}

#***************************************************************************
# Subroutine:  fasta_to_delimited
# Description: 
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
# Subroutine:  delimited_to_fasta
# Description: 
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
# Subroutine:  filter_sequences 
# Description: handler fxn handing off to specific filtering fxns 
#***************************************************************************
sub filter_sequences {

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
# Subroutine:  sort_sequences 
# Description: handler fxn handing off to specific sorting fxns 
#***************************************************************************
sub sort_sequences {

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
		#$fileio->write_output_file($outfile, \@sorted);
	}
	else {
		my $output_ref = $self->{output};
		unless ($output_ref) { die; }
		my $outfile = 'sorted_seqs.fas';
		my $outpath = $report_dir . $outfile;
		$fileio->write_output_file($outpath, \@sorted);
		my $filelink  = "<p>Click <a href='$outpath'>here</a> to retrieve the";
		if ($self->{filtered}) {
			$filelink .= " filtered,";
		}
		$filelink .= " sorted sequences <BR><BR>";
		push (@$output_ref, $filelink);
	}
}

############################################################################
# FASTA TOOL SORTING FUNCTIONS
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
# Subroutine:  sort_seqs_by_length 
# Description: sort a fasta file of sequences in order of descending length
# Arguments: $file (a path)
#            $fasta_ref (array to store output)
#            $count_gaps (flag to count gap characters)
#***************************************************************************
sub sort_seqs_by_length2 {

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
		else {
			$store = 1; 
		}
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
		$self->load_reference_sequence($refseq_name);
		$refseq_seq   = $self->{refseq}->{sequence};
		$refseq_name  = $self->{refseq}->{name};
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
	$fileio->write_output_file($outpath, \@fasta);

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

############################################################################
# FILTERING SEQUENCES
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
# Subroutine:  sort sequences using tab data 
# Description: Sort Seqs by Any Category
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
		$fileio->write_output_file($outfile, \@fasta);
	}
}

############################################################################
# COUNTING SEQUENCES
############################################################################

#***************************************************************************
# Subroutine:  count fasta sequences
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
# REFORMAT HEADERS
############################################################################

#***************************************************************************
# Subroutine:  truncate_long_headers 
# Description: truncate long headers in a FASTA file 
# Arguments:  
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
	$fileio->write_output_file($fasta_out, \@converted);
	$fileio->write_output_file($table_out, \@table);
}

#***************************************************************************
# Subroutine:  truncate_long_headers 
# Description: truncate long headers in a FASTA file 
# Arguments:  
#***************************************************************************
sub truncate_long_headers2 {

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
		my $header   = $sequence_ref->{header};
		my $sequence = $sequence_ref->{sequence};
		my @bits = split (/=/, $header);
		my $truncate1 = shift @bits;
		my $truncate2 = shift @bits;
		my @bits2 = split (/\(/, $truncate2);
		my @bits3 = split (/\|/, $truncate1);
		my $truncate3 = shift @bits2;
		my $truncate4 = pop @bits3;
		$truncate3  =~ s/\s+/_/g;
		$truncate4  =~ s/\s+/_/g;
		print "\n\t TRUNCATE: $truncate3 AND $truncate4";
		my $new_id = $truncate3 . $truncate4;
		$new_id  =~ s/\//_/g;
		$new_id  =~ s/\'/-/g;
		print "\n\t NEW ID: $new_id";
		
		my $fasta = ">$new_id\n$sequence\n\n";
		push (@converted, $fasta);
		my $id_pair = "$new_id\t$header";
		push (@table, $id_pair);
	}
	my $fasta_out = $file . '.converted.fas';
	my $table_out = $file . '.table.txt';
	$fileio->write_output_file($fasta_out, \@converted);
	$fileio->write_output_file($table_out, \@table);
}

############################################################################
# EXTRACT & SPLIT
############################################################################

#***************************************************************************
# Subroutine:  split_fasta 
# Description: 
#***************************************************************************
sub split_fasta {

	my ($self) = @_;
	
	# get the file name and read
	my @fasta;
	my $question = "\n\t What is the name of the file? ";
	$console->ask_input_file_question($question, \@fasta);
	
	my $file;
	my @sequence;
	foreach my $line (@fasta) {

		if ($line =~ '>') {		
			if (@sequence) {
				$fileio->write_output_file($file, \@sequence);
				@sequence = ();
			}
			push(@sequence, $line);
			$file = $line;
			chomp $file;
			$file =~ s/^>//;
			$file =~ s/\s//;
			$file .= '.fa';
		
		}
		else { 
			push (@sequence, $line); 
		}
	}
	$fileio->write_output_file($file, \@sequence);
}

#***************************************************************************
# Subroutine:  extract_random_seqs
# Description: 
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
	$fileio->write_output_file($file, \@seqs);
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
	$fileio->write_output_file($output_file, \@captured_seqs);

	# write out the rejected sequences
	my $reject_file = 'unmatched_from_' . $fasta_file;
	$fileio->write_output_file($reject_file, \@rejected_seqs);

	##### SUPERFLUOUS
	# write out the extracted sequence ids
	my $verify = 'extracted.txt';
	$fileio->write_output_file($verify, \@matched_ids);

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
	$fileio->write_output_file($outfile, \@extracted);
}

############################################################################
# LOS ALAMOS HIV DB
############################################################################

#***************************************************************************
# Subroutine:  reformat los alamos headers
# Description: convert Los Alamos headers (see examples below) to 
#              HIV_Central formatted headers.
#    
#	Los Alamos headers look like: >B.FR.83.HXB2_K03455
#	Their general format is >subtype.country.year.patient_accession
#	Think this is consistent for all seqs but have included some validation 
#   in case not.
#***************************************************************************
sub reformat_LA_headers {
	
	my ($self) = @_;
	
	# open the file
	my @los_alamos;
	my $question = "\n\t Name of the FASTA file with Los Alamos sequences?";
	my $la_file = $console->ask_input_file_question($question, \@los_alamos);

	# we don't use subtypes in the header but keep that data for a separate
	# comparsion
	my @la_subtypes;
	
	foreach my $line (@los_alamos) {

		if ($line =~ /^>/)      { 
			
			chomp $line;
			
			print "\nOri line: $line";
			my $la_id = $line;
			
			# extract data elements from header
			my @values = split(/\./, $line);
			
			my $source     = 'Los_Alamos';
			my $subtype    = $values[0];
			my $country    = $values[1];
			my $year       = $values[2];
			my $patient_id = $values[3];
			
			# and year to sample date (take midpoint since we only have
			# the year value - i.e. XXXX-06-15
			print "\n\tthe year is $year";
			my $samp_date;
			if    ($year =~ /-/) { $samp_date = '0000-00-00';            }
			elsif ($year < 30)   { $samp_date = '20' . $year . '-06-15'; }
			else                 { $samp_date = '19' . $year . '-06-15'; }
			
			# rewrite the header in HIV Central format
			$line  = $la_id;
			$line .= " [source=$source]";
			$line .= " [country=$country]";
			$line .= " [patient_id=$patient_id]";
			$line .= " [sample_date=$samp_date]";
			$line .= "\n";
			print "\nNew line: $line";

			# sort and store the subtype
			$la_id   =~ s/>//g;
			$subtype =~ s/>//g;
			push (@la_subtypes, "$la_id\t$subtype\n");
		}
	}
	
	# write out the reformatted file
	my $output_file = $la_file . '.updated';
	$fileio->write_output_file($output_file, \@los_alamos);
	
	# write out the subtype data 
	my $subtype_file = $la_file . '.subtypes';
	$fileio->write_output_file($subtype_file, \@la_subtypes);
}

#***************************************************************************
# Subroutine:  convert_los_alamos_epitope_tables_to_fasta
# Description: 
#***************************************************************************
sub convert_los_alamos_epitope_tables_to_fasta {
	
	my ($self, $directory_path) = @_;
	
	print "\n\t  # # # # # WATCH OUT: file must be in unix NOT dos format # # # # #\n";
	my $x;
	my @fasta;
	my @directory;die;

	foreach my $file (@directory) {

		$x++;
		print "\n\t FILE $x";
		
		my @file_bits = split('_', $file);
		my $file_stem = $file_bits[0];	
	
		my @csv_data;
		my $file_path = $directory_path . '/' . $file;	
		$fileio->read_input_file($file_path, \@csv_data);
		my @array;
		die;
		#$core->convert_csv_to_array_of_hash(\@csv_data, \@array);
		#$devtools->print_array(\@array);
		
		# Write out as fasta
		foreach my $hash_ref (@array) {
			my $start   = $hash_ref->{Hxb2locstart};	
			my $end     = $hash_ref->{Hxb2locend};	
			my $epitope = $hash_ref->{Epitope};	
			my $line = ">LosAlamos_$file_stem" . "_$start" . "_$end\n";
			$line .= "$epitope\n\n";
			push(@fasta, $line);
		}
	
		$fileio->write_output_file('Los_alamos_fasta.fas.txt', \@fasta);
	}
}

#***************************************************************************
# Subroutine:  convert_los_alamos_fasta_to_tabdelim
# Description: 
#***************************************************************************
sub convert_los_alamos_fasta_to_tabdelim {
	
	my ($self, $file_path) = @_;
	
	print "\n\t  # # # # # WATCH OUT: file must be in unix NOT dos format # # # # #\n";
	my @sequences;
	$seqio->read_fasta($file_path, \@sequences);

	my @tab_delim;
	push (@tab_delim, "Seq_ID\tSubtype\tCountry\tYear\n"); # Create header line
	my @new_fasta;
	foreach my $seq_ref (@sequences) {

		#$devtools->print_hash($seq_ref);
		my $header = $seq_ref->{header};
		my $alias  = $seq_ref->{header};
		my $seq    = $seq_ref->{sequence};
		my @header = split (/\./, $header);
		#$devtools->print_array(\@header);
	
		my $subtype = $header[0];
		my $country = $header[1];
		my $year_raw = $header[2];
		my @year = split ('', $year_raw);
		my $i    = 0;
		my $year = '';
		foreach my $char (@year) {
			$i++;
			if ($i > 4) { next; }
			$year .= $char; 
		}
		my @row;
		push (@row, $alias);
		push (@row, $subtype);
		push (@row, $country);
		push (@row, $year);
		my $row = join("\t", @row);
		push (@tab_delim, "$row\n");
		my $fasta = ">$alias\n$seq\n";
		push(@new_fasta, $fasta); 
	}

	my $outfile = $file_path . ".txt";
	$fileio->write_output_file($outfile, \@tab_delim);

	#my $seqfile = $file_path . ".fas";
	#$fileio->write_output_file($seqfile, \@new_fasta);

}

############################################################################
# TAB-DELIMITD FILE MANIPULATION UTILITIES
############################################################################

#***************************************************************************
# Subroutine:  combine data
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
# Subroutine:  fasta_to_phylip
# Description: 
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
	
		# Create phylip id
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

############################################################################
# Format conversion: NEXUS
############################################################################

#***************************************************************************
# Subroutine:  fasta_to_nexus 
# Description: Convert aligned FASTA to NEXUS format 
# Arguments:   $data_ref: reference to an array with the fasta sequence data
#***************************************************************************
sub fasta_to_nexus {

	my ($self, $file, $data_ref, $table_ref) = @_;
	
	# extract fasta to hash
	my @sequences;
	$seqio->read_fasta($file, \@sequences);
	my $num_taxa = scalar @sequences;
	unless ($num_taxa) { die "\n\t NO SEQUENCES in '$file'\n\n\n"; }

	# Get the maximum sequence and taxa label length (number of chars), & number of taxa
	my $max_length       = 0; 
	my $max_label_length = 0;
	my %clean_sequences;
	foreach my $seq_ref (@sequences) {
		
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
	push (@$data_ref, "#NEXUS");
	push (@$data_ref, "Begin DATA;");
	push (@$data_ref, "Dimensions ntax=$num_taxa nchar=$num_char;");
	push (@$data_ref, "Format Datatype=Nucleotide Gap=-;");
	push (@$data_ref, "Matrix");

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

############################################################################
# Format conversion: PHYLIP
############################################################################

#***************************************************************************
# Subroutine:  phylip_to_fasta 
# Description: 
#***************************************************************************
sub phylip_to_fasta {

	my ($self, $infile) = @_;

	my @phylip;
	$fileio->read_input_file($infile, \@phylip);
	my $taxa_line;
	do {
		my $line  = shift @phylip;
		if    ($line =~ /^\s*$/)   { next; } # discard blank line
		elsif ($line =~ /^\s*#/)   { next; } # discard comment line 
		else  { 
			$taxa_line = $line;
		}
	} until ($taxa_line);

	my @fasta;
	my $inline;
	my $outfile = "$infile.fasta";
	foreach my $inline (@phylip) {
	
		chomp($inline);
		if    ($inline =~ /^\s*$/)   { next; } # discard blank line
		elsif ($inline =~ /^\s*#/)   { next; } # discard comment line 
		
		my @line = split(/\s+/, $inline);
		my $header;
		my $sequence;
		foreach my $bit (@line) {
			my $line_type = $seqio->determine_string_type($bit);	
			if ($line_type eq 'sequence') {
				$sequence .= $bit;

			}
			elsif ($line_type eq 'header') {
				$header .= $bit;
			}
			else {
				print "\n\t Line type '$bit' is type $line_type";
				die;
			}
		}
		my $fasta = ">$header\n$sequence\n\n";
		push (@fasta, $fasta);
	}
	$fileio->write_output_file($outfile, \@fasta);
}

#***************************************************************************
# Subroutine:  get header words
# Description: 
#***************************************************************************
sub get_header_words {

	my ($self, $seqs_ref, $data_ref, $delimiter ) = @_;

	# Reformat the sequences	
	foreach my $seq_ref (@$seqs_ref) {
		my $sequence_id  = $seq_ref->{sequence_id};
		my $header       = $seq_ref->{header};
		chomp $header;
	}	
}

#***************************************************************************
# Subroutine:  load_reference_data 
# Description: 
#***************************************************************************
sub load_reference_sequence {

	my ($self, $refseq_name) = @_;

	my $parser_obj  = RefSeqParser->new();
	my $flat_dir = $self->{refseq_use_path};
	unless ($flat_dir)  {die; }

	# Read reference sequence 
	my $refseq_file = $flat_dir . $refseq_name;
	my %params;
	$parser_obj->parse_refseq_flatfile($refseq_file, \%params);
	my $refseq = RefSeq->new(\%params);
	$self->{refseq} = $refseq;

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
# Subroutine:  by number
# Description: by number - for use with perl 'sort'  (cryptic but works) 
#***************************************************************************
sub by_number { $a <=> $b }	

############################################################################
# EOF
############################################################################

