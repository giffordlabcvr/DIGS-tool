#!/usr/bin/perl -w
############################################################################
# Script:      add_rv_taxonomy_table.pl 
# Description: worked example script for DIGS - adding an auxillary table
#              and populating with data derived from a tab-delimited file
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
my ($USAGE)  = "\n\t #### add_rv_taxonomy_table.pl :\n";
    $USAGE  .= "\n\t usage: $0 -m=[mode] -d=[data file path] -i=[ctl file path]\n";
  	$USAGE  .= "\n\t -m=1  create table"; 
  	$USAGE  .= "\n\t -m=2  flush table"; 
  	$USAGE  .= "\n\t -m=3  drop table"; 
 	$USAGE  .= "\n\n";

############################################################################
# Main program
############################################################################

# Read in options using GetOpt::Long
my $ctl_file   = undef;
my $data_file  = undef;
my $mode       = undef;
GetOptions ('data|d=s'   => \$data_file, 
		    'infile|i=s' => \$ctl_file,
		    'mode|m=i'   => \$mode,
) or die $USAGE;

unless ($data_file and $ctl_file) { die $USAGE; }

# Run script
main($ctl_file, $data_file);

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
	
	my ($infile) = @_;

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

	# Parse the 'SCREENDB' block
	$loader_obj->parse_screendb_block(\@ctl_file);
	
	# Load screening database (includes some MacroLineage Tables
	$loader_obj->set_screening_db();
	my $db = $loader_obj->{db};
	
	if ($mode eq 1) {
		# Create the table
		create_rv_taxonomy_table($db);
		load_rv_taxonomy_table($db);
	
		# Insert data
		insert_rv_taxonomy_data($db, $data_file);
	}

	elsif ($mode eq 2) {
		load_rv_taxonomy_table($db);
 		my $taxonomy_table = $db->{rv_taxonomy};
		$taxonomy_table->flush();
		$taxonomy_table->reset_primary_keys();
	}

	elsif ($mode eq 3) {

		# Get connection variables from self
		my $db_name  = $db->{db_name};
		my $server   = $db->{server};
		my $username = $db->{username};
		my $password = $db->{password};
		unless ($server and $username and $password and $db_name) { die; }
		my $dbh = DBI->connect("dbi:mysql:$db_name:$server", $username, $password);
		unless ($dbh) {	die "\n\t # Couldn't connect to $db_name database\n\n"; }
		
		# Execute the query
		my $drop = "DROP TABLE 'RV_taxonomy' ";
		my $sth = $dbh->prepare($drop);
    	unless ($sth->execute()) { print $drop; exit;}	
	}
	else { die $USAGE; }

}

#***************************************************************************
# Subroutine:  insert_rv_taxonomy_data
# Description: 
#***************************************************************************
sub insert_rv_taxonomy_data {
	
	my ($db_obj, $data_file) = @_;
	
	# Insert RV taxonomy data
	my @rv_data;
	my $data_valid = $fileio->read_file($data_file, \@rv_data);
	unless ($data_valid) {
		print "\n\t ### Couldn't open data file '$data_file'\n\n\n ";
		exit;
	}	

 	my $taxonomy_table = $db_obj->{rv_taxonomy};
	$taxonomy_table->flush();
	$taxonomy_table->reset_primary_keys();
	
	my $header_line = shift @rv_data; # remove column headings
	foreach my $line (@rv_data) {

		my @data = split("\t", $line);
		my %data;
		$data{name}       = $data[0];
		$data{family}     = $data[1];
		$data{subfamily}  = $data[2];
		$data{supertribe} = $data[3];
		$data{tribe}      = $data[4];
		$data{genus}      = $data[5];

 		my $taxonony_table = $db_obj->{rv_taxonomy};
		$taxonony_table->insert_row(\%data);
	}
}

#***************************************************************************
# Subroutine:  create_rv_taxonomy_table
# Description: 
#***************************************************************************
sub create_rv_taxonomy_table {
	
	my ($db_obj) = @_;

	#$devtools->print_hash($db_obj); 
	my $dbh = $db_obj->{dbh}; 

	# Create RV taxonomy table
	my $rv_taxonomy = "CREATE TABLE `RV_taxonomy` (
	  `Record_ID`     int(11) NOT NULL auto_increment,

	  `Name`          varchar(100) NOT NULL default '0',
	  `Family`        varchar(100) NOT NULL default '0',
	  `Subfamily`     varchar(100) NOT NULL default '0',
	  `Supertribe`    varchar(100) NOT NULL default '0',
	  `Tribe`         varchar(100) NOT NULL default '0',
	  `Genus`         varchar(100) NOT NULL default '0',

	  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($rv_taxonomy);
	unless ($sth->execute()) { print "\n\t Couldn't create table - does it exist already?\n\n"}

}

#***************************************************************************
# Subroutine:  load_rv_taxonomy_table
# Description: load screening database table 'BLAST_results'
#***************************************************************************
sub load_rv_taxonomy_table {

	my ($db_obj) = @_;
	
	my $dbh = $db_obj->{dbh}; 
	
	# Definition of the table
	my %taxonomy_fields = (
		name        => 'varchar',
		family      => 'varchar',
		subfamily   => 'varchar',
		tribe       => 'varchar',
		supertribe  => 'varchar',
		genus       => 'varchar',
	);
	my $taxonomy_table = MySQLtable->new('RV_taxonomy', $dbh, \%taxonomy_fields);
	$db_obj->{rv_taxonomy} = $taxonomy_table;
}

#***************************************************************************
# Subroutine:  show_title
# Description: show command line title blurb 
#***************************************************************************
sub show_title {

	$console->refresh();
	my $title       = 'add_rv_taxonomy_table.pl';
	my $version     = '1.0';
	my $description = 'Worked example script: add auxillary table to DIGS database';
	my $author      = 'Robert J. Gifford';
	my $contact	    = '<robert.gifford@glasgow.ac.uk>';
	$console->show_about_box($title, $version, $description, $author, $contact);
}

############################################################################
# EOF
############################################################################
