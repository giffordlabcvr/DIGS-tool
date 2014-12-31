#!/usr/bin/perl -w
############################################################################
# Script:       RefSeqLibrary.pm 
# Description:  Functions for working with GLUE alignments 
# History:      Rob Gifford, December 2014: Creation
############################################################################
package RefSeqLibrary;

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
# Description: Create a new reference sequence library object
#***************************************************************************
sub new {

	my ($invocant, $params_ref) = @_;
	
	my $class = ref($invocant) || $invocant;
		
	my $self = {
		

	};

	bless ($self, $class);
	return $self;
}

############################################################################
# Get alignment statistics
############################################################################

#***************************************************************************
# Subroutine:  summarise_reference_library
# Description: summarise a library of GLUE reference sequences
#***************************************************************************
sub summarise_reference_library {

	my ($self, $path) = @_;

	# Load the reference sequences
	my %reference_sequences;
	$self->load_glue_reference_library(\%reference_sequences);
	die;

}

#***************************************************************************
# Subroutine:  load_glue_reference_library
# Description: read GLUE reference sequences into a hash with name as key
# Libraries are flat files, they can be read in different ways:
# 	from a flat directory
# 	from a directory with subdirectories (reads all files)
# concatenated GLUE refseqs in individual files can be read 
# refseq names must be unique
#***************************************************************************
sub load_glue_reference_library {
	
	my ($self, $refseqs_ref) = @_;


	my $ref_vglue   = $self->{reference_vglue};
	my $parser_obj  = RefSeqParser->new();

	my @refseqs;





}

#***************************************************************************
# Subroutine:  load_glue_reference_library_old
# Description: 
#***************************************************************************
sub load_glue_reference_library_old {
	
	my ($self) = @_;

	my $ref_vglue   = $self->{reference_vglue};
	my $report_dir  = $self->{output_path};
	my $tmp_path    = $self->{tmp_path};
	
	# Set up library
	my $parser_obj  = RefSeqParser->new();
	my @refseqs;
	my $status = $parser_obj->split_vglue_ref_file($ref_vglue, \@refseqs);
	unless ($status) { die "\n\t Input error: couldn't open query VGLUE file\n\n"; }
	my $num_refseqs = scalar @refseqs;
	print "\n\t Total of '$num_refseqs' GLUE formatted references sequences to use as references";
	my $i = 0;
	foreach my $refseq_ref (@refseqs) {
		
		$i++;
		
		# Step 1: parse the refseq file and capture data in %params
		my %params;
		$parser_obj->parse_refseq_datastructure($refseq_ref, \%params);

		# Step 2: create the reference sequence object using %params
		my $refseq = RefSeq->new(\%params);
		my $virus_name = $refseq->{name};
		print "\n\t Loading query refseq $i:  '$virus_name'";
		
		# Get the translated ORFs
		my %orf_sequences;
		$refseq->get_translated_orfs(\%orf_sequences); # Get orfs nucleic acid seqs 
		
		# Get the translated ORFs
		my %utr_sequences;
		$refseq->get_utrs(\%utr_sequences); # Get orfs nucleic acid seqs 
		#$devtools->print_hash(\%utr_sequences); die;

		my @orf_names = keys %orf_sequences;
		my @ref_fasta_aa;
		foreach my $orf_name (@orf_names) {
			my $orf_seq = $orf_sequences{$orf_name};
			my $name = $virus_name . "_$orf_name";
			my $fasta = ">$name\n$orf_seq\n\n";
			push (@ref_fasta_aa, $fasta);
			my $lib_path = $report_dir . "/reference_lib_nt.fas";
			$fileio->write_output_file($lib_path, \@ref_fasta_aa); 
			my $formatdb_command = "./bin/blast/makeblastdb -in $lib_path";
			system $formatdb_command;
			#print "\n\t$formatdb_command\n";
			$self->{blast_orf_lib_path} = $lib_path; 
		}
		
		# Iterate through UTRs and add those
		my @utr_names = keys %utr_sequences;
		my @ref_fasta_nt;
		foreach my $utr_name (@utr_names) {
			my $utr_seq = $utr_sequences{$utr_name};
			my $name = $virus_name . "_$utr_name";
			my $fasta = ">$name\n$utr_seq\n\n";
			push (@ref_fasta_nt, $fasta);
			my $lib_path = $report_dir . "/reference_lib_nt.fas";
			$fileio->write_output_file($lib_path, \@ref_fasta_nt); 
			my $formatdb_command = "./bin/blast/makeblastdb -in $lib_path";
			system $formatdb_command;
			#print "\n\t$formatdb_command\n";
			$self->{blast_utr_lib_path} = $lib_path; 
		}
	}
}

#***************************************************************************
# Subroutine:  by number
# Description: by number - for use with perl 'sort'  (cryptic but works) 
#***************************************************************************
sub by_number { $a <=> $b }	

############################################################################
# EOF
############################################################################
