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

# BLAST interface
use Interface::BLAST;

# GLUE program modules
use GLUE::DataTool;

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

# External program paths
my $blast_bin_path = './bin/blast/'; # BLAST

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

# Data tool
my %tool_params;
$tool_params{process_id}  = $process_id;
$tool_params{output_type} = 'text';
$tool_params{output_path} = './proc/';
$tool_params{refseq_use_path} = './db/refseq/';
my $datatool_obj = DataTool->new(\%tool_params);

############################################################################
# Set up USAGE statement
############################################################################

# Initialise usage statement to print if usage is incorrect
my($USAGE) = "#### Data_tool.pl: A collection of utilities for working with sequences + data\n";
  $USAGE  .= "\n # Convert between formats";
  $USAGE  .= "\n\t\t  -m=1        :     FASTA to Delimited";
  $USAGE  .= "\n\t\t  -m=2        : Delimited to FASTA";
  $USAGE  .= "\n\t\t  -m=3        :   Genbank to FASTA+DATA";
  $USAGE  .= "\n\t\t  -m=4        :     FASTA to NEXUS";
  $USAGE  .= "\n\t\t  -m=5        :     FASTA to PHYLIP";
#  $USAGE  .= "\n\t\t  -m=6        :    PHYLIP to FASTA*";
#  $USAGE  .= "\n\t\t  -m=7        :     NEXUS to FASTA*";
#  $USAGE  .= "\n # GLUE reference sequence utilities";
#  $USAGE  .= "\n\t\t  -g=1        :   Genbank to REFSEQ"; 
#  $USAGE  .= "\n\t\t  -g=2        :    REFSEQ to FASTA+DATA"; 
#  $USAGE  .= "\n\t\t  -g=3        :    REFSEQ to 'Pretty'"; 
#  $USAGE  .= "\n\t\t  -g=3        :    Get REFSEQ ORFs"; 
#  $USAGE  .= "\n # Manage data + sequences";
#  $USAGE  .= "\n\t\t  -d=1        :   Extract/filter/split sequences";
#  $USAGE  .= "\n\t\t  -d=2        :   Sort sequences"; 
#  $USAGE  .= "\n\t\t  -d=3        :   Data utilities";
  $USAGE  .= "\n\n  usage: $0 -m=[options] -i=[infile]\n\n";

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
	my $infile       = undef;
	my $mode		 = undef;
	my $glue		 = undef;
	my $data		 = undef;
	GetOptions ('infile|i=s'  => \$infile, 
				'mode|m=i'    => \$mode,
				'glue|g=i'    => \$glue,
				'data|d=i'    => \$data,
	) or die $USAGE;
	unless ($infile and $mode) { die $USAGE; }
	
	$datatool_obj->show_title();

	if ($mode) {    # Reformatting tools  

		if ($mode > 7) { die $USAGE; }
		$datatool_obj->run_reformat_tools_cmd_line($infile, $mode);
	}
	elsif ($glue) { # GLUE refseq tools

		if ($glue > 3) { die $USAGE; }
		$datatool_obj->run_refseq_tools_cmd_line($infile, $glue);
	}
	elsif ($data) { # Data tools

		if ($data eq 1)    {  # Sequence extraction fxns
			$datatool_obj->run_extract_tools_cmd_line($infile, $data);
		}
		elsif ($data eq 2) {  # Sequence sorting tools
			$datatool_obj->run_sort_tools_cmd_line($infile, $data);
		}
		elsif ($data eq 3) {  # Data + sequence tools
			$datatool_obj->run_data_tools_cmd_line($infile, $data);
		}
		else {
			die $USAGE;
		}
	}
	print "\n\n\t # Finished!\n\n\n";
}

############################################################################
# EOF
############################################################################
