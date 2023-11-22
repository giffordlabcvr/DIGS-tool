#!/usr/bin/perl -w
############################################################################
# Script:      get_digs_headers.pl
# Creator:     R.J. Gifford
# Description: DIGS for EVEs perl tools
# History:     Version 1.0
############################################################################

# Include the PERL module library for DIGS 
use lib ($ENV{DIGS_HOME} . '/modules/'); 

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;
use Getopt::Long;
use DBI;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base modules
use Base::Console;
use Base::FileIO;
use Base::DevTools;

############################################################################
# Paths & Globals
############################################################################

# Create a unique process ID for this DIGS screening process
my $pid  = $$;
my $time = time;
my $process_id  = $pid . '_' . $time;
my $version = '1.0';

############################################################################
# Instantiations
############################################################################

# Base utilites
my $fileio     = FileIO->new();
my $console    = Console->new();
my $devtools   = DevTools->new();

############################################################################
# Set up USAGE statement
############################################################################

# Initialise usage statement to print if usage is incorrect
my ($USAGE) = "\n\t  usage: $0 m=[option] -h=[help]\n\n";

############################################################################
# Main program
############################################################################

# Run script
main();

# Exit script
print "\n\n\t # Exit\n\n";
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
	my $help     = undef;

	# Read in options using GetOpt::Long
	GetOptions ('mode|m=i'       => \$mode, 
			    'help'           => \$help,
			    'infile|i=s'     => \$infile,
	) or die $USAGE;

	show_title();
	
	# Set flags based on options received
	if ($mode) { 
		if ($mode eq 1) {
			#format_ncbi_fasta_for_digs();
			format_ncbi_fasta_for_digs2($infile);
		}
	}
	elsif ($help) { # Show help page
		show_help_page();
		exit;
	}
	else {
		die $USAGE;
	}
	
}

#***************************************************************************
# Subroutine:  format_ncbi_fasta_for_digs2
# Description: simple header conversion routine
#***************************************************************************
sub format_ncbi_fasta_for_digs2 {

	my ($infile) = @_;

	my $species_name;
	my $species_name_question = "\n\t What is the name of the species?";
	$species_name = $console->ask_question($species_name_question);

	my @fasta;
	my @digs_fasta;
	$fileio->read_fasta($infile, \@fasta);
	my $i;
	foreach my $seq_hash (@fasta) {
		
		# Extract information about this sequence
		$i++;
		my $header   = $seq_hash->{header};
		my $sequence = $seq_hash->{sequence};
		$header=~ s/\//-/g;
		my @split1 = split(/\[/, $header);
		#$devtools->print_array(\@split1); #die;
		shift @split1;
		my $gene_element = shift @split1;
		$gene_element =~ s/]//g;
		$gene_element =~ s/\s+$//;
		my @gene_element = split('=', $gene_element);
		my $gene = pop @gene_element;
	   	      
		# Normalise gene name
		#my $normalised_gene = normalise_gene_names($name, $gene);

		# Create fasta
		#print "\n\t # $i DIGS $name ($normalised_gene)";
		print "\n\t # $i DIGS NAME: $species_name ($gene)";
		my $digs_header = $species_name . '_' . $gene;
		my $digs_fasta = ">$digs_header\n$sequence\n";	
		push (@digs_fasta, $digs_fasta);

	}

	my $outfile = $infile . '.out.txt';
	$fileio->write_file($outfile, \@digs_fasta);
	
}

############################################################################
# Title and help display
############################################################################

#***************************************************************************
# Subroutine:  show_title
# Description: show command line title blurb 
#***************************************************************************
sub show_title {

	$console->refresh();
	my $title       = 'get_digs_headers.pl';
	my $description = "Reformat NCBI 'coding sequence' FASTA for use with DIGS";
	my $author      = 'Robert J. Gifford';
	my $contact	    = '<robert.gifford@glasgow.ac.uk>';
	$console->show_about_box($title, $version, $description, $author, $contact);
}

#***************************************************************************
# Subroutine:  show_help_page
# Description: show help page information
#***************************************************************************
sub show_help_page {

	# Initialise usage statement to print if usage is incorrect
	my ($HELP)  = "\n\t Usage: $0 -m=[option] -i=[input file]\n";
        $HELP  .= "\n\t ### Main functions\n"; 
        $HELP  .= "\n\t -m=1  parse ncbi's 'coding sequence' FASTA to DIGS input FASTA\n\n"; 

	print "\n\t ### The input file for this script should be FASTA-formatted coding sequences downloaded from NCBI GenBank\n";

	print $HELP;
}

############################################################################
# End of file 
############################################################################
