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
		
		# Flags
		process_id            => $params_ref->{process_id},
		reference_glue        => $params_ref->{reference_glue},
		query_glue            => $params_ref->{query_glue},

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

	my ($self) = @_;

	# Load the reference sequences
	my %refseq_library;
	my $path = $self->{reference_glue};
	$self->load_glue_reference_library(\%refseq_library, $path);

	# Write as simple list
	my $exists = $fileio->check_directory_exists($path);
	my @leaves;
	if ($exists) {
		$fileio->read_directory_tree_leaves_simple($path, \@leaves);
	}
	else { die "\n\t # Couldn't open directory '$path'\n\n\n"; }
	#$devtools->print_array(\@leaves); die;
	my @list;
	foreach my $hash (@leaves) {
		my $path = $hash->{path};
		$path =~ s/\//\t/g;
		push(@list, "$path\n");
	}
	$fileio->write_file("tab.txt", \@list);

	# Write features as FASTA
	$self->write_features_as_fasta(\%refseq_library);

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
	
	my ($self, $refseqs_ref, $path) = @_;

	unless ($path) { die; }

	my $exists = $fileio->check_directory_exists($path);
	my @leaves;
	if ($exists) {
		#print "\n\t Reading leaves for PATH $path";
		$fileio->read_directory_tree_leaves_simple($path, \@leaves);
	}
	else {
		die "\n\t # Couldn't open directory '$path'\n\n\n";
	}
	#$devtools->print_hash($self); die;

	my $parser_obj  = RefSeqParser->new();
	my $i = 0;
	foreach my $file_ref (@leaves) {
	
		#my $status = $parser_obj->split_glue_ref_file($query_glue, \@refseqs);
		#unless ($status) { die "\n\t Input error: couldn't open query GLUE file\n\n"; }
		#my $num_refseqs = scalar @refseqs;
		#$devtools->print_array(\@refseqs); die;
	
		# Parse the refseq file and capture data in %params
		$i++;
		my $refseq_path = $file_ref->{path};
		my $refseq_file = $file_ref->{file};
		my %params;
		$parser_obj->parse_refseq_flatfile($refseq_path, \%params);

		# Create the reference sequence object using %params
		my $refseq = RefSeq->new(\%params);
		my $virus_name = $refseq->{name};
		print "\n\t Loading GLUE reference sequence $i:  '$virus_name'";
		my $refseq_name = $refseq->{name};
		$refseqs_ref->{$refseq_name} = $refseq;
	}	
}

	
#***************************************************************************
# Subroutine:  write_features_as_fasta
# Description: summarise a library of GLUE reference sequences
#***************************************************************************
sub write_features_as_fasta {

	my ($self, $refseq_library) = @_;

	my @refseq_names = keys %$refseq_library;
	my %features;
	my %feature_names;
	foreach my $refseq_name (@refseq_names) {
	
		#$refseq_ref->describe();
		#$devtools->print_hash($refseq_ref); die;
		print "\n\t Getting features for '$refseq_name'";
		my $refseq = $refseq_library->{$refseq_name};
		my %utrs;
		my %na_orfs;
		my %aa_orfs;
		$refseq->get_utrs(\%utrs);
		$refseq->get_orfs(\%na_orfs);		
		$refseq->get_translated_orfs(\%aa_orfs);		
		my @utrs    = keys %utrs;
		my @na_orfs = keys %na_orfs;
		my @aa_orfs = keys %aa_orfs;
		my @features = (@utrs, @na_orfs, @aa_orfs);
		foreach my $feature (@features) {
			my $lc_feature = lc $feature;
			$feature_names{$feature} = $lc_feature;
		}
		#$devtools->print_hash(\%feature_names); die;	
		
		my %refseq_features;
		$refseq_features{utrs}    = \%utrs;
		$refseq_features{na_orfs} = \%na_orfs;
		$refseq_features{aa_orfs} = \%aa_orfs;
		my $refseq_name = $refseq->{name};
		$features{$refseq_name} = \%refseq_features;
	}

	my @features = keys %feature_names;
	my %utr_fasta;
	my %na_orf_fasta;
	my %aa_orf_fasta;
	foreach my $feature (@features) {
		
		#print "\n\t Storing feature $feature";
		foreach my $refseq_name (@refseq_names) {

			my $features = $features{$refseq_name};
			my $utrs_ref = $features->{utrs};
			my $na_orfs  = $features->{na_orfs};
			my $aa_orfs  = $features->{aa_orfs};
			#$devtools->print_hash($utrs_ref); die;	

			my $utr    = $utrs_ref->{$feature};
			my $na_orf = $na_orfs->{$feature};
			my $aa_orf = $aa_orfs->{$feature};
	
			my $lc_feature = lc $feature;	
			if ($utr) {
				my %seq_ref;
				$seq_ref{header} = $refseq_name . '_' . $lc_feature;
				$seq_ref{sequence} = $utr;

				if ($utr_fasta{$feature}) {
					my $fasta_array = $utr_fasta{$feature};
					push (@$fasta_array, \%seq_ref);
				}
				else {
					my @fasta_array;
					push (@fasta_array, \%seq_ref);
					$utr_fasta{$feature} = \@fasta_array; # No lower case for UTRS
				}
			}

			if ($na_orf) {
				my %seq_ref;
				$seq_ref{header} = $refseq_name . '_' . $lc_feature;
				$seq_ref{sequence} = $na_orf;

				if ($na_orf_fasta{$lc_feature}) {
					my $fasta_array = $na_orf_fasta{$lc_feature};
					push (@$fasta_array, \%seq_ref);
				}
				else {
					my @fasta_array;
					push (@fasta_array, \%seq_ref);
					$na_orf_fasta{$lc_feature} = \@fasta_array;
				}
			}
			if ($aa_orf) {
				my %seq_ref;
				$seq_ref{header} = $refseq_name . '_' . $lc_feature;
				$seq_ref{sequence} = $aa_orf;

				if ($aa_orf_fasta{$lc_feature}) {
					my $fasta_array = $aa_orf_fasta{$lc_feature};
					push (@$fasta_array, \%seq_ref);
				}
				else {
					my @fasta_array;
					push (@fasta_array, \%seq_ref);
					$aa_orf_fasta{$lc_feature} = \@fasta_array;
				}
			}
		}
	}

	# Create a structured directory for the output
	my $process_id = $self->{process_id};
	my $directory  = 'fasta_dump_' . $process_id;
	my $output_cmd = "mkdir $directory";
	system $output_cmd;

	my $na_subdir = "$directory/na_orfs";
	my $subdir_cmd = "mkdir $na_subdir";
	system $subdir_cmd;

	my $aa_subdir = "$directory/aa_orfs";
	$subdir_cmd = "mkdir $aa_subdir";
	system $subdir_cmd;

	my $utr_subdir = "$directory/utrs";
	$subdir_cmd = "mkdir $utr_subdir";
	system $subdir_cmd;

	# Write the files
	foreach my $feature (@features) {

		#print "\n\t Write loop: $feature";	
		if ($utr_fasta{$feature}) {
			#print "\n\t Writing FASTA for UTR feature '$feature'";
			my $utrs_ref = $utr_fasta{$feature};
			my $outpath = $utr_subdir . '/' . $feature . '_utr.fa';
			$seqio->write_fasta($outpath, $utrs_ref);
		}

		if ($aa_orf_fasta{$feature}) {
			#print "\n\t Writing amino acid FASTA for ORF '$feature'";
			my $aa_orfs_ref = $aa_orf_fasta{$feature};
			my $outpath = $aa_subdir . '/' . $feature . '_aa.fa';
			$seqio->write_fasta($outpath, $aa_orfs_ref);
		}

		if ($na_orf_fasta{$feature}) {
			#print "\n\t Writing nucleic acid FASTA for ORF '$feature'";
			my $na_orfs_ref = $na_orf_fasta{$feature};
			my $outpath = $na_subdir . '/' . $feature . '_na.fa';
			$seqio->write_fasta($outpath, $na_orfs_ref);
		}
	}
	print "\n\t # Done!\n\n";

}

############################################################################
# Deprecated
############################################################################

############################################################################
# EOF
############################################################################
