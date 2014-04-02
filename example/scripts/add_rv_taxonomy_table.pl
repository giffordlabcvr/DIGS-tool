#!/usr/bin/perl -w
############################################################################
# Script:      add_rv_taxonomy_table.pl 
# Description: worked example script for DIGS - adding an auxillary table
#              and populating with data derived from a tab-delimited file
# History:     Version 1.0 Creation: Rob J Gifford 2014
############################################################################

unless ($ENV{DIGS}) {
	print  "\n\n\t Environment variable '\$DIGS' is undefined\n";
	print  "(path to genome data directory)\n\n\n";
	exit;
}
# Include a local library of PERL modules 
use lib ($ENV{DIGS}) . '/modules/'; 

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;
use Getopt::Long;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base modules
use Base::Console;
use Base::DevTools;
use Base::FileIO;

############################################################################
# Paths & Globals
############################################################################


############################################################################
# Instantiations for program 'classes' (PERL's Object-Oriented Emulation)
############################################################################

# Base utilites
my $seqio      = SeqIO->new();
my $fileio     = FileIO->new();
my $devtools   = DevTools->new();
my $console    = Console->new();

############################################################################
# Set up USAGE statement
############################################################################

# Initialise usage statement to print if usage is incorrect
my ($USAGE)  = "\n\t #### add_rv_taxonomy_table.pl :\n";
    $USAGE  .= "\n\t usage: $0 -i=[data file path]\n";
 	$USAGE  .= "\n\n";

############################################################################
# Main program
############################################################################

# Run script
main();

# Exit program
print "\n\n\t DONE\n\n\n";
exit;

############################################################################
# Subroutines
############################################################################

#***************************************************************************
# Subroutine:  main
# Description: main fxn of this script
#***************************************************************************
sub main {
	
	# Show title
	show_title();

	# Read the input
	my @rv_data;
	read_input(\@rv_data);

	# Create the table
	create_rv_taxonomy_table();

	# Enter the data
	insert_rv_taxonomy_data(\@rv_data);

}

#***************************************************************************
# Subroutine:  read_input
# Description: 
#***************************************************************************
sub read_input {
	
	# Read the input


}

#***************************************************************************
# Subroutine:  create_rv_taxonomy_table
# Description: 
#***************************************************************************
sub create_rv_taxonomy_table {
	
	# Create RV taxonomy table


}

#***************************************************************************
# Subroutine:  insert_rv_taxonomy_data
# Description: 
#***************************************************************************
sub insert_rv_taxonomy_data {
	
	# Insert RV taxonomy data


}

############################################################################
# End of file 
############################################################################
