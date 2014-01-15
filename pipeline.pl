#!/usr/bin/perl -w
############################################################################
# Script:      pipeline.pl 
# Description: control script for performing database-guided genome screening
# History:     Version 1.0 Creation: Rob J Gifford 2014
############################################################################

# Include a local library of PERL modules 
use lib './modules/'; 
unless ($ENV{GENOMES}) {
	print  "\n\n\t PLEASE DEFINE '\$GENOMES' (Path to genome data directory)\n\n\n";
	exit;
}

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

# Third party program interface modules
use Interface::BLAST;   # Interface to BLAST 

# Paleovirology program modules
use DGGS::Pipeline;
use DGGS::ScreenBuild;
use DGGS::GenomeControl;
use DGGS::DB;

############################################################################
# Paths & Globals
############################################################################

# Paths
my $blast_bin_path        = '';                  # path to directory with BLAST+ programs
                                                 # leave blank if BLAST+ programs in path 
#my $blast_bin_path       = './bin/blast/';      
my $genome_use_path       = $ENV{GENOMES};       # Genome data directory
my $output_path           = './process/';        # Process directory

# Database connection 
my $mysql_server   = 'localhost'; # SET ON INSTALL IF YOU WANT TO USE A REMOTE DB
my $mysql_username = 'root';      # SET ON INSTALL IF YOU WANT A NON-ROOT USER
my $mysql_password = '';          # SET THIS ON INSTALL
	
# Process ID and time - used to create a unique ID for each program run
my $pid  = $$;
my $time = time;
my $process_id   = $pid . '_' . $time;

############################################################################
# Instantiations for program 'classes' (PERL's Object-Oriented Emulation)
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
$params{process_id}         = $process_id;
$params{blast_bin_path}     = $blast_bin_path; 
$params{genome_use_path}    = $genome_use_path;
$params{output_path}        = $output_path; 
$params{server}             = $mysql_server;
$params{username}           = $mysql_username;
$params{password}           = $mysql_password;
$params{blast_obj}          = $blast_obj;
my $pipeline_obj = Pipeline->new(\%params);

############################################################################
# Set up USAGE statement
############################################################################

# Initialise usage statement to print if usage is incorrect
my($USAGE) = "\n #### Screening:";
  $USAGE  .= "\n\t\t  -s=[ctl file]             : Run screen using ctl file";
  $USAGE  .= "\n\t\t  -u=[option] -i=[ctl file] : Utility functions";
  $USAGE  .= "\n\n  usage: $0 [options] \n\n";

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
	my $screen  = undef;
	my $utils   = undef;
	my $infile  = undef;
	GetOptions ('screen|s=s' => \$screen, 
	            'utils|u=s'  => \$utils, 
			    'infile|i=s' => \$infile,
	) or die $USAGE;

	# Hand off to the appropriate object
	if ($screen) {
		$pipeline_obj->run_screen($screen); 
	}
	elsif ($utils) {
		$pipeline_obj->run_screen_utility_function($utils, $infile); 
	}
	else { # command line script called without arguments
		$console->refresh();
		die $USAGE; 
	}
}

############################################################################
# EOF 
############################################################################
