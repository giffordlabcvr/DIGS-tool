#!/usr/bin/perl -w
############################################################################
# Script:      data_tool.pl 
# Description: a collection of tools for working with sequences + data
# History:     Version 1.0 Creation: Rob J Gifford 2014
############################################################################

# use a local modules for this program (for portability)
use lib './modules/'; 

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;
use CGI;
use Getopt::Long;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base modules
use Base::Console;
use Base::DevTools;
use Base::FileIO;
use Base::SeqIO;

# GLUE program modules
use GLUE::Data_tool;

############################################################################
# Paths
############################################################################

############################################################################
# Globals
############################################################################

# Process ID and time - used to create a unique ID for each program run
my $pid  = $$;
my $time = time;

# Create a unique ID for this process
my $process_id   = $pid . '_' . $time;
my $user = $ENV{"USER"};

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
my($USAGE) = "\n #### Data tool - a collection of utilities for working with sequences + data\n";
  $USAGE  .= "\n\t\t  -m=1        : FASTA to delimited";
  $USAGE  .= "\n\t\t  -m=2        : FASTA to NEXUS";
  $USAGE  .= "\n\t\t  -m=3        : FASTA to PHYLIP";
  $USAGE  .= "\n\t\t  -m=4        : Delimited to FASTA";
  $USAGE  .= "\n\t\t  -m=5        : Extract by identifier";
  $USAGE  .= "\n\t\t  -m=6        : Filter by identifier";
  $USAGE  .= "\n\n  usage: $0 [options] \n\n";

############################################################################
# Main program
############################################################################

# Run script
$console->refresh();
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
	my $mode = undef;
	GetOptions ('mode|m=i' => \$mode,
	) or die $USAGE;
	unless ($mode) { die $USAGE; }

	# Hand off to the appropriate object
	if ($mode eq 1) { 
	
	}
	else { 
		die $USAGE; 
	}
	print "\n\n\t # Finished!\n\n\n";
}

############################################################################
# EOF 
############################################################################
