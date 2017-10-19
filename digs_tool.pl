#!/usr/bin/perl -w
############################################################################
# Script:      digs_tool.pl database-integrated genome screening (DIGS) tool
# Description: A tool for exploring genomes 'in silico' using BLAST and 
#              a relational database
# History:     Updated: February 2017
############################################################################

# Capture variables set in the environment 
unless ($ENV{'DIGS_GENOMES'}) {
	print  "\n\n\t Required environment variable '\$DIGS_GENOMES' is undefined\n";
	exit;
}
unless ($ENV{'DIGS_HOME'}) {
	print  "\n\n\t Required environment variable '\$DIGS_HOME' is undefined\n";
	exit;
}
unless ($ENV{'DIGS_MYSQL_USER'}) {
	print  "\n\n\t Required environment variable '\$DIGS_MYSQL_USER' is undefined\n";
	exit;
}
unless ($ENV{'DIGS_MYSQL_PASSWORD'}) {
	print  "\n\n\t Required environment variable '\$DIGS_MYSQL_PASSWORD' is undefined\n";
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
use DIGS::Test;
use DIGS::Utility;
use DIGS::Nomenclature;

############################################################################
# Paths & Globals
############################################################################

# Paths
my $genome_use_path = $ENV{DIGS_GENOMES} . '/'; 
$genome_use_path =~ s/\/\//\//g; # Remove any double backslashes
my $blast_bin_path  = '';  # left empty if BLAST+ programs are in your path 

# Version number	
my $program_version = '1.13.2';

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

my $mysql_username = ($ENV{DIGS_MYSQL_USER}); 
my $mysql_password = ($ENV{DIGS_MYSQL_PASSWORD}); 

# Instantiate main program classes using global settings
my %params;
$params{program_version} = $program_version;
$params{process_id}      = $process_id;
$params{blast_bin_path}  = $blast_bin_path; 
$params{genome_use_path} = $genome_use_path;
$params{blast_obj}       = $blast_obj;
$params{mysql_username}  = $mysql_username ; 
$params{mysql_password}  = $mysql_password; 

my $digs_tool_obj = DIGS->new(\%params);

############################################################################
# Set up USAGE statement
############################################################################

# Initialise usage statement to print if usage is incorrect
my $USAGE = "\n\t ### DIGS version $program_version";
  $USAGE .= "\n\t ### usage: $0 m=[option] -i=[control file] -h=[help]\n\n";

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

	# Set version
	
	# Options that require a file path
	my $infile      = undef;
	
	# Options that require a numerical value
	my $mode         = undef;
	my $database     = undef;
	my $utility      = undef;
	my $genomes      = undef;
	
	# Options that don't require a value
	my $help         = undef;
	my $extra_help   = undef;
	my $verbose      = undef;
	my $force        = undef;
	my $test         = undef;
	my $create_ids   = undef;
	
	# Read in options using GetOpt::Long
	GetOptions ('infile|i=s'     => \$infile,
	
	            'mode|m=i'       => \$mode,
			    'database|d=i'   => \$database,
			    'utility=i'      => \$utility,
			    'genomes=i'      => \$genomes,
			    'create_ids'     => \$create_ids,			      
			    'extra_help'     => \$extra_help,
			    'verbose'        => \$verbose,
			    'force'          => \$force,
			    'help'           => \$help,
			    'test'           => \$test,
			    			    
	) or die $USAGE;

	# Set flags based on options received
	if ($verbose) {  $digs_tool_obj->{verbose} = 'true'; }
	if ($force)   {  $digs_tool_obj->{force} = 'true';   }

	# Hand off to functions based on options received
	if ($help) { # Show help page
		$digs_tool_obj->show_help_page();
		exit;
	}
	elsif ($test) { # Run inbuilt tests
		my $test_obj = Test->new($digs_tool_obj);
		$test_obj->show_test_validation_options();
	}
	elsif ($extra_help) { # Show help for utility options
		my $utility_obj = Utility->new($digs_tool_obj);
		$utility_obj->show_utility_help_page();
		exit;
	}
	elsif ($mode) { # Main DIGS tool functions 
		$digs_tool_obj->run_digs_process($infile, $mode); 
	}
	elsif ($create_ids) { # Apply nomenclature rules to loci
		my $nomenclature_obj = Nomenclature->new($digs_tool_obj);
		$nomenclature_obj->create_standard_locus_ids($infile);
		
	}
	elsif ($database or $utility or $genomes or $utility) { # Utility functions
		my $utility_obj = Utility->new($digs_tool_obj);
		$utility_obj->run_utility_process($infile, $database, $genomes, $utility); 
	}	
	else { die $USAGE; }

	# Exit script
	print "\n\n\t # Exit\n\n";
}

############################################################################
# End of file 
############################################################################
