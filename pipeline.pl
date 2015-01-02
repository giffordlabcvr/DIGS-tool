#!/usr/bin/perl -w
############################################################################
# Script:      pipeline.pl 
# Description: control script for database-integrated genome screening (DIGS)
# History:     Version 1.0 Creation: Rob J Gifford 2014
############################################################################

unless ($ENV{DIGS}) {
	print  "\n\n\t Environment variable '\$DIGS' is undefined\n";
	print  "(path to this directory)\n\n\n";
	exit;
}
unless ($ENV{GENOMES}) {
	print  "\n\n\t Environment variable '\$GENOMES' is undefined\n";
	print  "(path to directory containing target sequence files)\n\n\n";
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

# Third party program interface modules
use Interface::BLAST;   # Interface to BLAST 

# DIGS framework modules
use DIGS::Pipeline;
use DIGS::ScreenBuild;
use DIGS::TargetDB;
use DIGS::ScreeningDB;

# GLUE framework modules
use GLUE::RefSeqLibrary; # Functions to set up screen
use GLUE::RefSeqParser;  # GLUE RefSeq parsing
use GLUE::RefSeq;        # GLUE RefSeq

############################################################################
# Paths & Globals
############################################################################

# Paths
#my $blast_bin_path       = $ENV{DIGS} . '/bin/blast/';  # Example
my $blast_bin_path        = '';              # Path to directory with BLAST+ programs
                                             # leave blank if BLAST+ programs in path 
#my $genome_use_path       = $ENV{DIGS} . '/targets/';   # genome data directory
my $genome_use_path       = $ENV{GENOMES} . '/';        # genome data directory
my $output_path           = $ENV{DIGS} . '/proc/';      # default process directory
	
# Process ID and time - used to create a unique ID for each program run
my $pid  = $$;
my $time = time;
my $process_id  = $pid . '_' . $time;
my $program_version = '1.0';

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
$params{output_path}        = $output_path; 
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
		$pipeline_obj->run_digs_function($mode, $infile); 
		print "\n\t # Exit\n\n";
	}
	else {
		die $USAGE;
	}
}

############################################################################
# End of file 
############################################################################
