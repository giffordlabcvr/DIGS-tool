#!/usr/bin/perl -w
############################################################################
# Module:      BLAST.pm 
# Description: Interface to BLAST
# History:     Rob Gifford January 2007: Creation
############################################################################
package BLAST;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::FileIO;

############################################################################
# Globals
############################################################################

my $fileio   = FileIO->new();
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Create a new BLAST.pm 'object'
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Member variables
	my $self = {
	    blast_bin_path  => $parameter_ref->{blast_bin_path},
	};

	bless ($self, $class);
	return $self;
}

############################################################################
# Top level fxns
############################################################################

#***************************************************************************
# Subroutine:  BLAST
# Description: Execute BLAST a search
#***************************************************************************
sub blast {

	my ($self, $method, $target_path, $probe_path, $result_path, $options_ref) = @_;
		
	# Get paths from self
	my $blast_path  = $self->{blast_bin_path};
	my $word_size   = $options_ref->{word_size}; 
	my $evalue		= $options_ref->{evalue};
	my $penalty     = $options_ref->{penalty};
	my $reward      = $options_ref->{reward};
	my $gapopen     = $options_ref->{gapopen};
	my $gapextend   = $options_ref->{gapextend};
	my $outfmt      = $options_ref->{outfmt};
	unless ($outfmt) { $outfmt = 7; } # default output format is tab-delimited ( -outfmt 7) 
	
	# Create BLAST command
	$blast_path .= $method;
	my $blast_type;  
	my $set_params;  
	if ($options_ref) {
		$blast_type = $options_ref->{blast_type};
		$set_params = $options_ref->{set_params};
	}
	my $command;
	if ($blast_type) {
		$command  = "$blast_path -query $probe_path -subject $target_path ";
		$command .= "-out $result_path ";
	}
	else {
		$command  = "$blast_path -query $probe_path -db $target_path ";
		$command .= "-out $result_path ";
	}

	if ($set_params) {  
		if ($word_size)  { $command .= " -word_size $word_size "; }
		if ($evalue)     { $command .= " -evalue $evalue "; } 
		if ($penalty and $reward and $gapopen and $gapextend)    { 
			$command .= " -penalty $penalty ";  
			$command .= " -reward $reward ";  
			$command .= " -gapopen $gapopen "; 
			$command .= " -gapextend $gapextend "; 
		} 
	}

	# Set the output format for BLAST
	$command .= "-outfmt $outfmt";

	# Execute the command 
	system $command;		
}

#***************************************************************************
# Subroutine:  extract_sequence
# Description: Interface to the BLAST+ sequence extraction functions 
#              see: http://www.ncbi.nlm.nih.gov/books/NBK1763/
#***************************************************************************
sub extract_sequence {

	my ($self, $target_path, $hit_ref) = @_;

	# Get path to BLAST binary
	my $blast_path  = $self->{blast_bin_path};
	
	# Get hit parameters
	my $start       = $hit_ref->{subject_start};
	my $end         = $hit_ref->{subject_end};
	my $orientation = $hit_ref->{orientation};
	my $scaffold    = $hit_ref->{scaffold};
	unless ($start and $end and $orientation and $scaffold and $target_path) { die; }

	# Parsing for blastdbcmd
	my @gi = split(/\|/,$scaffold);	
	if (scalar(@gi) > 1) {
		$scaffold = $gi[1];
	}

	# Create the command
	# Command example: 
	# /bin/blast/blastdbcmd -db hs_alt_HuRef_chrX.fa -entry 157734237 
	# -range 10-60 -strand minus
	my $command = $blast_path . "blastdbcmd -db $target_path";
	$command .= " -entry $scaffold ";
	$command .= " -range $start-$end ";
	if ($orientation eq '-') { $command .= ' -strand minus '; }
	
	# Execute the command
	my @sequence = `$command`;
	shift @sequence;  # Remove header
	my $sequence = join ('', @sequence);
	$sequence =~ s/\n//g;

	# Check we got sequences
	unless ($sequence) {
		print "\n\t ## FAILED TO EXTRACT sequence using command '$command'\n";
	}
	
	return $sequence;
}

############################################################################
# Lower level parsing fxns
############################################################################

#***************************************************************************
# Subroutine:  parse_tab_format_results
# Description: parse BLAST results
#***************************************************************************
sub parse_tab_format_results {

	my ($self, $file, $result_ref, $bitscore_cutoff) = @_;

	# Read the file into an array
	my $fileio = FileIO->new();
	my @file;
	$fileio->read_file($file, \@file);

	# Process the file
	my @matches;
	foreach my $line (@file) {

		if ($line =~ /^#/) { next; } # skip comments
		my @data = split("\t", $line);

		# GET THE DATA FOR THIS HIT
		my %match;
		$match{probe}         = $data[0];
		$match{scaffold}      = $data[1];
		$match{identity}      = $data[2];
		$match{align_len}     = $data[3];
		$match{mismatches}    = $data[4];
		my $gap_openings      = $data[5];
		unless ($gap_openings) {
			$gap_openings = '0';
		}
		$match{gap_openings}  = $gap_openings;
		$match{query_start}   = $data[6];
		$match{query_stop}    = $data[7];
		my $aln_start         = $data[8];
		my $aln_stop          = $data[9];
		my $e_value           = $data[10];
		my $bit_score         = $data[11];
		chomp $bit_score;
		
		# APPLY CONDITIONS
		if ($bitscore_cutoff) {
			if ($bit_score < $bitscore_cutoff) {
				# Skip matches with bitscores below threshold
				next;
			}
		}

		# Convert e value
		$self->convert_evalue($e_value, \%match);
		
		# Set orientaton
		my $orientation;
		if  ($aln_stop < $aln_start) {
			$orientation = '-';
			
			# switch the start and stop around if in -ve orientation
			my $switch_start = $aln_stop; 
			my $switch_stop   = $aln_start; 
			$match{aln_start} = $switch_start;
			$match{aln_stop}   = $switch_stop;
		}
		else {
			$orientation = '+';
			$match{aln_start} = $aln_start;
			$match{aln_stop}   = $aln_stop;
		}
		$match{orientation}   = $orientation; 
		$match{bit_score}     = $bit_score;
		push(@$result_ref, \%match);
	}		
}

#***************************************************************************
# Subroutine:  convert_evalue
# Description: split evalue string into number and exponential 
#***************************************************************************
sub convert_evalue { 

	my ($self, $e_value, $data_ref) = @_;
	
	my $e_value_num = 0;
	my $e_value_exp = 0;
	if ($e_value =~ /e/) {
		my @e_value_bits = split ("e-", $e_value);
		$e_value_num = $e_value_bits[0];
		$e_value_exp = $e_value_bits[1];
	}
	else {
		$e_value_num = $e_value;
		$e_value_exp = 1;
	}
	$data_ref->{e_value_num} = $e_value_num;
	$data_ref->{e_value_exp} = $e_value_exp;
}	

############################################################################
# EOF
############################################################################
