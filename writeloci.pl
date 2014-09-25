#!/usr/bin/perl -w
############################################################################
# Script:      writeloci.pl 
# Description: apply Missillac rules to Loci table from DIGS screening 
# History:     Version 1.0 Creation: Rob J Gifford 2014
############################################################################

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;
use Getopt::Long;

############################################################################
# Set up USAGE statement
############################################################################

# Initialise usage statement to print if usage is incorrect
my($USAGE) = "#### writeloci.pl -i=[infile] -m=[mode] (default=1)\n";

############################################################################
# Main program
############################################################################

# Run script
main();

# Exit program
exit;

############################################################################
# Subroutines
############################################################################

#***************************************************************************
# Subroutine:  main
# Description: top level handler fxn
#***************************************************************************
sub main {

	# Read in options using GetOpt::Long
	my $infile       = undef;
	my $mode		 = undef;
	GetOptions ('infile|i=s'  => \$infile, 
				'mode|m=i'    => \$mode,
	) or die $USAGE;
	unless ($infile and $mode) { die $USAGE; }
	
	# Read the file 
	my @raw_file;
	read_file($infile, \@raw_file);

	# Process the file line by line
	my @output;
	my $erv_number = 0;
	foreach my $line (@raw_file) {
	
		chomp $line; # remove newline
		#print "\n\t LINE: $line";

		my @line = split("\t", $line); # Split on tabs
		
		# Get the data from the line 
		# (note this assumes a certain column order)
		my $record_id     = $line[0];
		my $organism      = $line[1];
		my $erv_name      = $line[2];
		my $extract_start = $line[3];
		my $extract_end   = $line[4];
		my $orientation   = $line[5];
		my $chunk_name    = $line[6];
		my $genome_struc  = $line[7];

		$erv_number++; # Increment counter by one

		# Assign the group part of the name
		my $group = get_group_from_ERV_name($erv_name);

		# Assign the genome part of the name
		#my $genome = get_genome_name_from_ERV_organism($organism);
		#HACK - only using human loci
		my $genome = 'HoSa';

		# Create the locus text
		my @locus;
		push (@locus, 'ERV');
		push (@locus, $group);
		push (@locus, $erv_number);
		my $locus_name = join("-", @locus);
		$line .= "\t$locus_name\n";
		push (@output, $line);
	}
	
	# Write the output
	my $outfile = $infile . '.loci.txt';
	write_file($outfile, \@output);

	print "\n\n\t # Finished!\n\n\n";
}

#***************************************************************************
# Subroutine:  get_group_from_ERV_name
# Description: self explanatory
#***************************************************************************
sub get_group_from_ERV_name {

	my ($erv_name) = @_;

	# In line with the Missilac ERV edict
	# ERVs will be assigned to one of the following three groups
	# Gamma-like
	# Beta-like
	# Spuma-like

	# Define variables for the group names
	my $group1 = 'Gamma-like';
	my $group2 = 'Beta-like';
	my $group3 = 'Spuma-like';

	my $group = undef;
	if    ($erv_name eq 'CERV-1')      { $group = $group1; }
	elsif ($erv_name eq 'CERV-2')      { $group = $group1; }
	elsif ($erv_name eq 'HERV-T')      { $group = $group1; }
	elsif ($erv_name eq 'HERV-E')      { $group = $group1; }
	elsif ($erv_name eq 'RR-HERV-I')   { $group = $group1; }
	elsif ($erv_name eq 'HERV-Fb')     { $group = $group1; }
	elsif ($erv_name eq 'HERV-XA')     { $group = $group1; }
	elsif ($erv_name eq 'HERV-H')      { $group = $group1; }
	elsif ($erv_name eq 'HERV-I')      { $group = $group1; }
	elsif ($erv_name eq 'HERV-R')      { $group = $group1; }
	elsif ($erv_name eq 'HERV-W')      { $group = $group1; }
	elsif ($erv_name eq 'ERV-9')       { $group = $group1; }
	elsif ($erv_name eq 'HERV-30')      { $group = $group1; }
	elsif ($erv_name eq 'HERV-K-14C')   { $group = $group2; }
	elsif ($erv_name eq 'HERV-K-HML2')  { $group = $group2; }
	elsif ($erv_name eq 'HERV-K-HML4')  { $group = $group2; }
	elsif ($erv_name eq 'HERV-K-HML5')  { $group = $group2; }
	elsif ($erv_name eq 'HERV-K-HML6')  { $group = $group2; }
	elsif ($erv_name eq 'HERV-K-HML7')  { $group = $group2; }
	elsif ($erv_name eq 'HERV-K-HML8')  { $group = $group2; }
	elsif ($erv_name eq 'HERV-K-HML9')  { $group = $group2; }
	elsif ($erv_name eq 'RhERV-K-HML2') { $group = $group2; }
	elsif ($erv_name eq 'HERV-L')      { $group = $group3; }
	elsif ($erv_name eq 'HERV-S')      { $group = $group3; }
	
	return $group;
}

#***************************************************************************
# Subroutine:  read_file
# Description: read an input file to an array
# Arguments:   $file: the name of the file to read
#              $array_ref: array to copy to
#***************************************************************************
sub read_file {

	my ($file, $array_ref) = @_;

	unless (-f $file) {
		if (-d $file) {
			print "\n\t Cannot open file \"$file\" - it is a directory\n\n";
			return 0;
		}
		else {
			print "\n\t Cannot open file \"$file\"\n\n";
			return 0;
		}

 	}
	unless (open(INFILE, "$file")) {
		print "\n\t Cannot open file \"$file\"\n\n";
		return 0;
	}
	@$array_ref = <INFILE>;
	close INFILE;

	return 1;
}


#***************************************************************************
# Subroutine:  write_file
# Description: write an array to an ouput file
# Arguments:   $file: the name of the file to write to 
#              $array_ref: array to copy
#***************************************************************************
sub write_file {

	my ($self, $file, $array_ref) = @_;
	unless (open(OUTFILE, ">$file")) {
		print "\n\t Couldn't open file \"$file\" for writing\n\n";
		return 0;
	}
	print OUTFILE @$array_ref;
	close OUTFILE;
	print "\n\t File \"$file\" created!\n\n";
}

############################################################################
# EOF
############################################################################
