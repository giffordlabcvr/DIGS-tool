#!/usr/bin/perl -w
############################################################################
# Script:       Water.pm 
# Description:  Interface to Water
# History:      August 2012: Creation
############################################################################
package Water;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::Core;
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
# Description: 
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Member variables
	my $self = {
	    water_bin_path  => $parameter_ref->{water_bin_path},
	};

	bless ($self, $class);
	return $self;
}

############################################################################
# Public member functions
############################################################################

#***************************************************************************
# Subroutine:  run_water
# Description: 
#***************************************************************************
sub run_water {
	
	my ($self, $seq1_path, $seq2_path, $result_path) = @_;

	# Get data from self and show title
	my $water_path = $self->{water_bin_path};
	unless ($water_path) { die; }
	my $command = "$water_path  $seq1_path $seq2_path";
	   $command .= " -gapopen 5 -gapextend 2 $result_path &> /dev/null";
	#print "\n\t $command \n\n"; die;
	system $command;
}

#***************************************************************************
# Subroutine:  parse_water
# Description: 
#***************************************************************************
sub parse_water {
	
	my ($self, $result_path, $data_ref) = @_;
	
	my @waterout_lines;
	$fileio->read_input_file($result_path, \@waterout_lines);
	
	# Should only be one entry in array
	#my $line_with_identity = $waterout_lines_with_identity[0]; 
	my $identity;
	my $length_of_match;
	my $query_start;
	foreach my $line (@waterout_lines) {
		chomp $line;
		if ($line =~ /Identity/) {
			#print "\nIdentity '$line'";
			my @split_line = split  (/\s+/, $line); 
			$identity = $split_line[3];
			$identity =~ s/%//g;
			$identity =~ s/\)//g;
			$identity =~ s/\(//g;
			$split_line[2] =~ m/\d+\/(\d+)/;
			$length_of_match = $1;
		}
		if ($identity) {
			unless (($line =~ /^#--/) or ($line =~ /^\s*$/)) {
				my @split_line = split  (/\s+/, $line); 
				$query_start = $split_line[1];
			}
		}
	}

	# divide by whitespaces and place into array
	$data_ref->{identity}        = $identity;
	$data_ref->{length_of_match} = $length_of_match;
	$data_ref->{query_start}     = $query_start;

}

#***************************************************************************
# Subroutine:  sort_water_pbs_results
# Description: 
#***************************************************************************
sub sort_water_pbs_results {
	
	my ($self, $data_ref, $lowest_ref) = @_;
	
	# Following copied from another script - just sorts according to div score
	my $lowest_id;
	my $lowest_div;
	my %lowest;
	my $i = 0;
	foreach my $result_ref (@$data_ref) {
		$i++;
		my $div      = $result_ref->{"div"};
		$lowest{$i} = $result_ref;
		unless ($lowest_div)       {
			$lowest_div = $div;
			$lowest_id  = $i; 
		}
		elsif ($div < $lowest_div) { 
			$lowest_div = $div;
			$lowest_id  = $i; 
		}
	}
	#my $lowest = $lowest{$lowest_id};
	my $this_lowest = $lowest{$lowest_id};
	$lowest_ref->{pbs_assign}  = $this_lowest->{"PBS"};
	$lowest_ref->{divergence}  = $this_lowest->{"div"};
	$lowest_ref->{identity}    = $this_lowest->{"identity"};
	$lowest_ref->{aln_len}     = $this_lowest->{"length_of_match"};
	$lowest_ref->{locus_id}    = $this_lowest->{locus_id};
	$lowest_ref->{assigned_to} = $this_lowest->{assigned_to};
	$lowest_ref->{query_start} = $this_lowest->{query_start};
	#$devtools->print_hash($this_lowest);
	$lowest_ref = $this_lowest;
	#$devtools->print_hash($lowest_ref);
}

############################################################################
# EOF
############################################################################
