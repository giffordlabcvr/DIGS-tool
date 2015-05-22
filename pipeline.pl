#!/usr/bin/perl -w
############################################################################
# Script:      pipeline.pl 
# Description: control script for database-integrated genome screening (DIGS)
# History:     Version 1.0 Creation: Rob J Gifford 2014
############################################################################

unless ($ENV{DIGS_HOME2}) {
	print  "\n\n\t Environment variable '\$DIGS_HOME2' is undefined\n";
	exit;
}
unless ($ENV{DIGS_GENOMES2}) {
	print  "\n\n\t Environment variable '\$DIGS_GENOMES2' is undefined\n";
	exit;
}

# Include the PERL module library for DIGS 
use lib ($ENV{DIGS_HOME2}) . '/modules/'; 

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

# Third party program interface modules
use Interface::BLAST;   # Interface to BLAST 
use Interface::MySQLtable;   # Interface to BLAST 

# DIGS framework modules
use DIGS::Pipeline;
use DIGS::ScreenBuild;
use DIGS::TargetDB;
use DIGS::ScreeningDB;
use DIGS::Consolidation;
# GLUE framework modules
use GLUE::RefSeqLibrary; # Functions to set up screen
use GLUE::RefSeqParser;  # GLUE RefSeq parsing
use GLUE::RefSeq;        # GLUE RefSeq

############################################################################
# Paths & Globals
############################################################################

# Paths
my $genome_use_path = $ENV{DIGS_GENOMES2} . '/';    
my $blast_bin_path  = '';  # leave blank if BLAST+ programs are in your path 

# Version number	
my $program_version = '1.0';

# Create a unique process ID for this DIGS screening process
my $pid  = $$;
my $time = time;
my $process_id  = $pid . '_' . $time;

############################################################################
# Instantiations
############################################################################

# Base utilites
my $seqio      = SeqIO->new();
my $fileio     = FileIO->new();
my $devtools   = DevTools->new();
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
my $pipeline_obj = Pipeline->new(\%params);

############################################################################
# Set up USAGE statement
############################################################################

# Initialise usage statement to print if usage is incorrect
my ($USAGE) = "\n\t  usage: $0 m=[option] -i=[control file]\n\n";

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
	
	# Read in options using GetOpt::Long
	my $mode    = undef;
	my $infile  = undef;
	my $utility = undef;
	my $help    = undef;
	my $version = undef;
	GetOptions ('mode|m=i'    => \$mode, 
			    'infile|i=s'  => \$infile,
			    'utility|u=i' => \$utility,
			    'help'        => \$help,
			    'version'     => \$version,
	) or die $USAGE;

	# Sanity checking for input 
	if ($help)    { # Show help page
		$pipeline_obj->show_help_page();
	}
	elsif ($version) {
		print "\n\t DIGS tool version '$program_version'\n\n"
	}
	elsif ($mode) { # Hand off to Pipeline.pm
		$pipeline_obj->run_digs_process($mode, $infile); 
	}
	elsif ($utility) { # Hand off to Pipeline.pm utility function
		$pipeline_obj->run_utility_function($utility, $infile); 
	}
	else {
		die $USAGE;
	}
	print "\n\n\t # Exit\n\n";
}

############################################################################
# End of file 
############################################################################
