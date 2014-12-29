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
#use Base::HTML_Utilities;

# BLAST interface
use Interface::BLAST;

# GLUE program modules
use GLUE::DataTool;
use GLUE::DataToolCmd;
#use GLUE::DataToolWeb;
use GLUE::RefSeq;
use GLUE::RefSeqParser;

############################################################################
# Paths
############################################################################

############################################################################
# Globals
############################################################################

# Version number
my $program_version = '1.0 beta';

# Process ID and time - used to create a unique ID for each program run
my $pid  = $$;
my $time = time;

# Create a unique ID for this process
my $process_id   = $pid . '_' . $time;
my $user = $ENV{"USER"};

# Paths
my $output_path           = './site/reports/';     # Reports
my $refseq_use_path       = './db/refseq_flat/';   # RefSeq flat directory

############################################################################
# External program paths
my $blast_bin_path = ''; # BLAST

############################################################################
# Instantiations for program 'classes' (PERL's Object-Oriented Emulation)
############################################################################

# Base utilites
my $seqio      = SeqIO->new();
my $fileio     = FileIO->new();
my $devtools   = DevTools->new();
my $console    = Console->new();

# Interface to BLAST
#my %blast_params;
#$blast_params{blast_bin_path} = $blast_bin_path;
#my $blast_obj = BLAST->new(\%blast_params);

# Data tool
my %tool_params;
$tool_params{process_id}  = $process_id;
$tool_params{output_type} = 'text';  # Default is text
$tool_params{output_path} = $output_path;
$tool_params{refseq_use_path} = './db/refseq/';
my $datatool = DataTool->new(\%tool_params);
$tool_params{datatool_obj} = $datatool;
my $cmd_line_interface = DataToolCmd->new(\%tool_params);

############################################################################
# Set up USAGE statement
############################################################################

# Initialise usage statement to print if usage is incorrect
my ($USAGE) = "\n\t  usage: $0 -[options]\n\n";

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

	# Define options
	my $help         = undef;
	my $version      = undef;
	my $mode		 = undef;
	my $sort		 = undef;
	my @seqfiles     = undef;
	my @datafiles    = undef;

	# Read in options using GetOpt::Long
	GetOptions ('help!'         => \$help,
                'version!'      => \$version,
				'mode|m=i'      => \$mode,
				'sort|s=i'      => \$sort,
				'seqfiles|f=s'  => \@seqfiles, 
				'datafiles|d=s' => \@datafiles, 
	);

	if ($help)    { # Show help page
		$cmd_line_interface->show_help_page();  
	}
	elsif ($version)  { 
		print "\n\t # GLUE datatool.pl version $program_version\n\n";  
	}
	elsif ($mode) { # Data reformatting tools
		my $result = $cmd_line_interface->run_reformat_tools_cmd_line(\@seqfiles, $mode);
	}
	elsif ($sort) { # Data sorting tools
		$cmd_line_interface->run_sort_tools_cmd_line(\@seqfiles, \@datafiles, $sort);
	}
	else {
		die $USAGE;
	}
	print "\n\n\t # Exit\n\n\n";
}

############################################################################
# EOF
############################################################################
