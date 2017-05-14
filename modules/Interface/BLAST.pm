#!/usr/bin/perl -w
############################################################################
# Module:      BLAST.pm 
# Description: A Perl interface to the BLAST executables
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
use Base::DevTools;

############################################################################
# Globals
############################################################################

my $fileio   = FileIO->new();
my $devtools = DevTools->new();
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
	    num_threads     => $parameter_ref->{num_threads},

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
	
	my $num_threads  = $self->{num_threads};
	if  ($num_threads) {
		$command .= "-num_threads $num_threads ";
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
	#print "\n\t COMMAND $command"; sleep 1; exit;
	system $command;
	#die;		
}

#***************************************************************************
# Subroutine:  extract_sequence
# Description: Interface to the BLAST+ sequence extraction functions 
#              see: http://www.ncbi.nlm.nih.gov/books/NBK1763/
#***************************************************************************
sub extract_sequence {

	my ($self, $target_path, $data_ref) = @_;

	# Get path to BLAST binary
	my $blast_path  = $self->{blast_bin_path};
	
	# Get extraction parameters
	my $start       = $data_ref->{start};
	my $end         = $data_ref->{end};
	my $orientation = $data_ref->{orientation};
	my $scaffold    = $data_ref->{scaffold};
	unless ($start and $end and $orientation and $scaffold and $target_path) {
		$devtools->print_hash($data_ref); die; 
	}

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

	my ($self, $file, $result_ref) = @_;

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
		my $bitscore         = $data[11];
		chomp $bitscore;
		
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
			$match{aln_stop}  = $switch_stop;
		}
		else {
			$orientation = '+';
			$match{aln_start} = $aln_start;
			$match{aln_stop}  = $aln_stop;
		}
		
		$match{orientation}   = $orientation; 
		$match{bitscore}     = $bitscore;
		push(@$result_ref, \%match);
	}		
}

#***************************************************************************
# Subroutine:  convert_evalue
# Description: split evalue string into number and exponential 
#***************************************************************************
sub convert_evalue { 

	my ($self, $e_value, $data_ref) = @_;
	
	my $evalue_num = 0;
	my $evalue_exp = 0;
	if ($e_value =~ /e/) {
		my @evalue_bits = split ("e-", $e_value);
		$evalue_num = $evalue_bits[0];
		$evalue_exp = $evalue_bits[1];
	}
	else {
		$evalue_num = $e_value;
		$evalue_exp = 1;
	}
	$data_ref->{evalue_num} = $evalue_num;
	$data_ref->{evalue_exp} = $evalue_exp;
}	

############################################################################
# EOF
############################################################################
