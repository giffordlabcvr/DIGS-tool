#!/usr/bin/perl -w
############################################################################
# Script:      digs_tool.pl database-integrated genome screening (DIGS) tool
# Description: A tool for exploring genomes 'in silico' using BLAST and 
#              a relational database
# History:     Version 1.12 December 2016
############################################################################

# Check the required environment variable are defined
unless ($ENV{DIGS_GENOMES}) {
	print  "\n\n\t Required environment variable '\$DIGS_GENOMES' is undefined\n";
	exit;
}
unless ($ENV{DIGS_HOME}) {
	print  "\n\n\t Required environment variable '\$DIGS_HOME' is undefined\n";
	exit;
}

# Include the PERL module library for DIGS 
use lib ($ENV{DIGS_HOME}) . '/modules/'; 

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
use Base::FileIO;

# Third party program interface modules
use Interface::BLAST;   # Interface to BLAST 
use Interface::MySQLtable;   # Interface to BLAST 

# DIGS framework modules
use DIGS::DIGS;
use DIGS::ScreenBuilder;
use DIGS::TargetDB;
use DIGS::ScreeningDB;

############################################################################
# Paths & Globals
############################################################################

# Paths
my $genome_use_path = $ENV{DIGS_GENOMES} . '/';    
my $blast_bin_path  = '';  # leave blank if BLAST+ programs are in your path 

# Version number	
my $program_version = '1.1';

# Create a unique process ID for this DIGS screening process
my $pid  = $$;
my $time = time;
my $process_id  = $pid . '_' . $time;

############################################################################
# Instantiations
############################################################################

# Base utilites
my $fileio     = FileIO->new();
my $console    = Console->new();

# Interface to BLAST
my %blast_params;
$blast_params{blast_bin_path} = $blast_bin_path;
my $blast_obj = BLAST->new(\%blast_params);

# Instantiate main program classes using global settings
my %params;
$params{program_version}    = $program_version;
$params{process_id}         = $process_id;
$params{blast_bin_path}     = $blast_bin_path; 
$params{genome_use_path}    = $genome_use_path;
$params{blast_obj}          = $blast_obj;
my $digs_tool_obj = DIGS->new(\%params);

############################################################################
# Set up USAGE statement
############################################################################

# Initialise usage statement to print if usage is incorrect
my ($USAGE) = "\n\t  usage: $0 m=[option] -i=[control file] -h=[help]\n\n";

############################################################################
# Main program
############################################################################

# Run script
main();
exit;

############################################################################
# Subroutines
############################################################################

#***************************************************************************
# Subroutine:  main
# Description: top level handler fxn
#***************************************************************************
sub main {
	
	# Options that require a file path
	my $infile   = undef;
	
	# Options that require a numerical value
	my $mode     = undef;
	my $utility  = undef;
	
	# Options that don't require a value
	my $help       = undef;
	my $extra_help = undef;
	my $version    = undef;
	my $analysis   = undef;
	my $test       = undef;
	
	# Read in options using GetOpt::Long
	GetOptions ('mode|m=i'      => \$mode,
			    'utility=i'     => \$utility,
			    'infile|i=s'    => \$infile,    
			    'analysis|a=i'  => \$analysis,
			    'help'          => \$help,
			    'extra_help'    => \$extra_help,
			    'version'       => \$version,
			    'test'          => \$test,
	) or die $USAGE;
	
	# Hand off to functions
	if ($version) { 
		print "\n\t DIGS tool version $program_version\n\n"
	}
	elsif ($help) { # Show help page
		$digs_tool_obj->show_help_page();
		exit;
	}
	elsif ($test) { # Run inbuilt tests
		$digs_tool_obj->validate();
		exit;
	}
	elsif ($extra_help) {
		$digs_tool_obj->show_utility_help_page();
		exit;
	}
	elsif ($mode) { # Main DIGS tool functions 
		$digs_tool_obj->run_digs_process($mode, $infile); 
	}
	elsif ($utility) { # Utility functions
		$digs_tool_obj->run_utility_process($utility, $infile); 
	}

	else { die $USAGE; }

	# Exit script
	print "\n\n\t # Exit\n\n";
}

############################################################################
# End of file 
############################################################################
