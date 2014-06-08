#!/usr/bin/perl -w
############################################################################
# Script:      scaffold_profiler.pl 
# Description: worked example script for DIGS
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

use DIGS::ScreenBuild;

############################################################################
# Paths & Globals
############################################################################


############################################################################
# Instantiations for program 'classes' (PERL's Object-Oriented Emulation)
############################################################################

# Base utilites
my $fileio     = FileIO->new();
my $devtools   = DevTools->new();
my $console    = Console->new();

############################################################################
# Set up USAGE statement
############################################################################

# Initialise usage statement to print if usage is incorrect
my ($USAGE)  = "\n\t #### scaffold_profiler.pl :\n";
    $USAGE  .= "\n\t usage: $0 -i=[control file path]\n";
 	$USAGE  .= "\n\n";

############################################################################
# Main program
############################################################################

# Read in options using GetOpt::Long
my $ctl_file   = undef;
GetOptions ('infile|i=s' => \$ctl_file,
) or die $USAGE;
unless ($ctl_file) { die $USAGE; }

# Run script
main($ctl_file);

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

	my ($ctl_file) = @_;
	
	# Show title
	show_title();

	# Try opening control file first
	my @ctl_file;
	my $valid = $fileio->read_file($ctl_file, \@ctl_file);
	unless ($valid) {
		print "\n\t ### Couldn't open control file '$ctl_file'\n\n\n ";
		exit;
	}	

	print "\n\t ### Reading control file\n";
	my $loader_obj = ScreenBuild->new();

	# Read input file
	$fileio->read_file($ctl_file, \@ctl_file);

	# Parse the 'SCREENDB' block
	$loader_obj->parse_screendb_block(\@ctl_file);
	
	# Load screening database (includes some MacroLineage Tables
	my $db_name = $loader_obj->{db_name};
	unless ($db_name) { die "\n\t Error: no DB name set \n\n\n"; }
	my $db = DB->new($loader_obj);
	$db->load_screening_db($db_name);

	# Do BLAST results chain
	#do_blast_results($db);

	# Do extracted chain
	do_extracted($db);

}

#***************************************************************************
# Subroutine:  do_extracted
# Description:
#***************************************************************************
sub do_extracted {

	my ($db) = @_;
	
	# Show title
	show_title();

	# Get hits by scaffold
	my $extracted_table = $db->{extracted_table};

	# Set the fields to get values for
	my @fields = qw [ record_id target_name scaffold orientation
	                  extract_start extract_end
                      query_start query_end sequence_length
                      assigned_name assigned_gene ];

	# Order by ascending start coordinates within each scaffold in the target file
	my $where = "ORDER BY organism, data_type, version, target_name, scaffold, extract_start ";
	
	# Get the relevant loci
	my @hits;
	$extracted_table->select_rows(\@fields, \@hits, $where);

	# Iterate through consolidating as we go
	my $i = undef;
	my $first_target = 1;
	my %last_hit;
	my @output;
	my $extract_count = 0;
	foreach my $hit_ref (@hits)  {

		my $consolidated;
		$extract_count++;
	
		# Get hit values
		my $record_id     = $hit_ref->{record_id};
		my $target_name   = $hit_ref->{target_name};
		my $scaffold      = $hit_ref->{scaffold};
		my $orientation   = $hit_ref->{orientation};
		my $extract_start = $hit_ref->{extract_start};
		my $extract_end   = $hit_ref->{extract_end};
		my $assigned_name = $hit_ref->{assigned_name};
		my $assigned_gene = $hit_ref->{assigned_gene};
		
		# Get last hit values
		my $last_record_id     = $last_hit{record_id};
		my $last_scaffold      = $last_hit{scaffold};
		my $last_target_name   = $last_hit{target_name};
		my $last_orientation   = $last_hit{orientation};
		my $last_extract_start = $last_hit{extract_start};
		my $last_extract_end   = $last_hit{extract_end};

		if ($first_target) {
			my $line  = "\n\n ### Target $target_name; Scaffold: '$scaffold'\n";
			push (@output, $line);
			$first_target = undef;
		}
		
		my $gap = 0;
		if ($last_extract_start) {
			$gap = $extract_start - $last_extract_end;
		}
		unless ($last_target_name) {
			$last_target_name = $target_name;
		}	
		
		# Is this a new target file?
		if ($target_name ne $last_target_name) {
			$i=1;
			my $line  = "\n\n ### Target $target_name; Scaffold: '$scaffold'";
			push (@output, $line);	
		}
		# Is this a new scaffold
		elsif ($last_scaffold) {	
			if ($scaffold ne $last_scaffold) {
				$i=1;
				my $line  = "\n\n ### Target $target_name; Scaffold: '$scaffold'\n";
				push (@output, $line);	
			}
			else { $i++; 
				my $line  = "\t\t Gap of $gap nucleotides";
				push (@output, $line);	
			}
		}

		# Compose the line	
		my $line  = "\nExtracted seq $extract_count at\t $target_name, $scaffold";
		   $line .= ",$extract_start,$extract_end";
		   $line .= "\t\t $assigned_name,$assigned_gene ($orientation)";
		
		if ($i) {
			if ($i > 1 and $gap < 1000) { 
				$line .= "\t  ### ARRAY HIT\n"; 
			}
		}
		#print $line;
		push (@output, $line);	

		# Update last hit data
		$last_hit{record_id}     = $record_id;
		$last_hit{target_name}   = $target_name;
		$last_hit{scaffold}      = $scaffold;
		$last_hit{orientation}   = $orientation;
		$last_hit{extract_start} = $extract_start;
		$last_hit{extract_end}   = $extract_end;
		$last_hit{assigned_name} = $assigned_name;
		$last_hit{assigned_gene}  = $assigned_gene;
	}
	
	my $db_name = $db->{db_name};
	my $outfile = $db_name . '_scaffold_profile.txt';
	$fileio->write_file($outfile, \@output);

}


#***************************************************************************
# Subroutine: do_blast_results 
# Description: main fxn of this script
#***************************************************************************
sub do_blast_results {

	my ($db) = @_;
	
	# Get hits by scaffold
	my $blast_results_table = $db->{blast_results_table};
	
	# Set the fields to get values for
	my @fields = qw [ record_id target_name scaffold orientation
	                  subject_start subject_end
                      query_start query_end hit_length ];

	# Order by ascending start coordinates within each scaffold in the target file
	my $where = "ORDER BY organism, target_name, scaffold, subject_start ";
	
	# Get the relevant loci
	my @hits;
	$blast_results_table->select_rows(\@fields, \@hits, $where);

	# Iterate through consolidating as we go
	my $i;
	my %last_hit;
	my @output;
	my $blast_count= 0;
	foreach my $hit_ref (@hits)  {

		# Get hit values
		my $record_id     = $hit_ref->{record_id};
		my $target_name   = $hit_ref->{target_name};
		my $scaffold      = $hit_ref->{scaffold};
		my $orientation   = $hit_ref->{orientation};
		my $subject_start = $hit_ref->{subject_start};
		my $subject_end   = $hit_ref->{subject_end};
		my $query_start   = $hit_ref->{query_start};
		my $query_end     = $hit_ref->{query_end};
		
		# Get last hit values
		my $last_record_id     = $last_hit{record_id};
		my $last_scaffold      = $last_hit{scaffold};
		my $last_target_name   = $last_hit{target_name};
		my $last_orientation   = $last_hit{orientation};
		my $last_subject_start = $last_hit{subject_start};
		my $last_subject_end   = $last_hit{subject_end};
		my $gap = $subject_start - $last_subject_end;
		$blast_count++;

		if ($target_name ne $last_target_name) {
			$i=1;
			my $line  = "\n\n ### Target $target_name; Scaffold: '$scaffold'";
			#print $line;
			push (@output, $line);	
		}

		elsif ($last_scaffold) {	
			if ($scaffold ne $last_scaffold) {
				$i=1;
				my $line  = "\n\n ### Target $target_name; Scaffold: '$scaffold'";
				#print $line;
				push (@output, $line);	
			}
			else {
				$i++; 
				my $line  = "\t\t Gap of $gap nucleotides";
				#print $line;
				push (@output, $line);	
			}
		}

		my $line  = "\n Hit $blast_count at:\t $target_name, $scaffold";
		   $line .= ",$subject_start,$subject_end ($orientation)";
		   $line .= ": query: $query_start, $query_end";
		if ($i) {
			if ($i > 1 and $gap < 1000) { 
				$line .= "\t  ### ARRAY HIT\n"; 
			}
		}
		push (@output, $line);	
		#print $line;

		# Update last hit data
		$last_hit{record_id}     = $record_id;
		$last_hit{target_name}   = $target_name;
		$last_hit{scaffold}      = $scaffold;
		$last_hit{orientation}   = $orientation;
		$last_hit{subject_start} = $subject_start;
		$last_hit{subject_end}   = $subject_end;
	}


	my $db_name = $db->{db_name};
	my $outfile = $db_name . '_scaffold_profile_BLAST-table.txt';
	$fileio->write_file($outfile, \@output);

}

#***************************************************************************
# Subroutine:  show_title
# Description: show command line title blurb 
#***************************************************************************
sub show_title {

	$console->refresh();
	my $title       = 'DIGS Scaffold Profiler';
	my $version     = '1.0';
	my $description = 'Profile hits in a screening database generated by DIGS';
	my $author      = 'Robert J. Gifford';
	my $contact	    = '<robert.gifford@glasgow.ac.uk>';
	$console->show_about_box($title, $version, $description, $author, $contact);
}

############################################################################
# EOF
############################################################################
