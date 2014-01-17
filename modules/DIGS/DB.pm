#!/usr/bin/perl -w
############################################################################
# Module:      DB.pm
# Description: The DB module for Pipeline databases
#              Contains routines for creating, and interacting with 
#              DBs in the Pipeline framework
# History:     June 2011: Created by Robert Gifford 
############################################################################
package DB;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;
use DBI;
############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::SeqIO;
use Base::FileIO;
use Base::DevTools;
use Base::Console;

# Components
use DIGS::GenomeControl;

# Database component modules
use Interface::MySQLtable;

############################################################################
# Globals
############################################################################

# Create base objects
my $seqio     = SeqIO->new();
my $fileio    = FileIO->new();
my $devtools  = DevTools->new();
my $console   = Console->new();
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Parameters
#***************************************************************************
sub new {

	my ($invocant, $parameters) = @_;

	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {
		
		# Paths
		db_name              => undef,
		
		# DB connection variables
		server                => $parameters->{server},
		username              => $parameters->{username},
		password              => $parameters->{password},
		
		# Screening database core tables
		status_table          => 0,
		blast_results_table   => 0,
		extracted_table       => 0,
	
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# LOAD DATABASES
############################################################################

#***************************************************************************
# Subroutine:  load_screening_db 
# Description: load database  
#***************************************************************************
sub load_screening_db {

	my ($self, $db_name) = @_;

	# Get connection variables from self
	my $server   = $self->{server};
	my $username = $self->{username};
	my $password = $self->{password};
	unless ($server and $username and $password and $db_name) { die; }

	# Set name
	$self->{db_name} = $db_name;
   	
	# Load tables from the screening result DB we've been given
	my $dbh = DBI->connect("dbi:mysql:$db_name:$server", $username, $password);
	unless ($dbh) {
        print "\n\n\t ### Creating '$db_name' screening database";
		$self->create_screening_db($db_name);
		$dbh = DBI->connect("dbi:mysql:$db_name:$server", $username, $password);
	}
	else {
		print "\n\n\t ### Loading '$db_name' screening database";
	}
	$self->{dbh} = $dbh;

	# Main Screening DB tables
	$self->load_blast_results_table($dbh);	
	$self->load_extracted_table($dbh);	
	$self->load_status_table($dbh);
}

############################################################################
# LOAD SCREENING DATABASE TABLES
############################################################################

#***************************************************************************
# Subroutine:  load_blast_results_table
# Description: load screening database table 'BLAST_results'
#***************************************************************************
sub load_blast_results_table {

	my ($self, $dbh) = @_;
	
	# Definition of the table
	my %blast_fields = (
		probe_name       => 'varchar',
		probe_gene       => 'varchar',
		probe_type       => 'varchar',
		organism         => 'varchar',
		chunk_name       => 'varchar',
		scaffold         => 'varchar',
		orientation      => 'varchar',
		bit_score        => 'float',
		identity         => 'varchar',
		e_value_num      => 'float',
		e_value_exp      => 'int',
		subject_start    => 'int',
		subject_end      => 'int',
		query_start      => 'int',
		query_end        => 'int',
		align_len        => 'int',
		mismatches       => 'int',
		gap_openings     => 'int',
	);
	my $blast_table = MySQLtable->new('BLAST_results', $dbh, \%blast_fields);
	$self->{blast_results_table} = $blast_table;
}

#***************************************************************************
# Subroutine:  load_extracted_table
# Description: load screening database table 'Extracted'
#***************************************************************************
sub load_extracted_table {

	my ($self, $dbh) = @_;

	# Definition of the table
	my %extract_fields = (
		blast_id         => 'blast_id',
		probe_name       => 'varchar',
		probe_gene       => 'varchar',
		probe_type       => 'varchar',
		assigned_to      => 'varchar',
		assigned_to_gene => 'varchar',
		probe_type       => 'varchar',
		organism         => 'varchar',
		chunk_name       => 'varchar',
		scaffold         => 'varchar',
		extract_start    => 'varchar',
		extract_end      => 'varchar',
		orientation      => 'varchar',
		bit_score        => 'float',
		identity         => 'varchar',
		e_value_num      => 'float',
		e_value_exp      => 'int',
	  	subject_start    => 'int',
	  	subject_end      => 'int',
		realex_start     => 'int',
        realex_end       => 'int',	
		query_start      => 'int',
	  	query_end        => 'int',
		align_len        => 'int',
		gap_openings     => 'int',
		mismatches       => 'int',
		sequence_length  => 'int',
		sequence         => 'text',
	);
	my $extract_table = MySQLtable->new('Extracted', $dbh, \%extract_fields);
	$self->{extracted_table} = $extract_table;
}

#***************************************************************************
# Subroutine:  load_status_table
# Description: load screening database table 'Status'
#***************************************************************************
sub load_status_table {

	my ($self, $dbh) = @_;

	# Definition of the table
	my %status_fields = (
		probe_id       => 'varchar',
		probe_name     => 'varchar',
		probe_gene     => 'varchar',
		organism       => 'varchar',
		version        => 'varchar',
		chunk_name     => 'varchar',
	);
	my $status_table = MySQLtable->new('Status', $dbh, \%status_fields);
	$self->{status_table} = $status_table;
}

############################################################################
# CREATE SCREENING DATABASE
############################################################################

#***************************************************************************
# Subroutine:  create_screening_db
# Description: create MySQL table objects for all screening DB tables
#***************************************************************************
sub create_screening_db {

	my ($self, $db_name) = @_;

	# Get connection variables from self
	my $server   = $self->{server};
	my $username = $self->{username};
	my $password = $self->{password};

	# CREATE THE DB 
    my $drh = DBI->install_driver("mysql");
    my $rc = $drh->func("createdb", $db_name, $server, $username, $password, 'admin');
    my $dbh = DBI->connect("dbi:mysql:$db_name:$server", $username, $password);
	unless ($dbh) {
		die "\n\t # Couldn't create screening database '$db_name'\n\n";
	}	

	$self->create_blast_results_table($dbh);
	$self->create_extracted_table($dbh);
	$self->create_status_table($dbh);
}

############################################################################
# CREATE SCREENING DATABASE TABLES
############################################################################

#***************************************************************************
# Subroutine:  create_blast_results_table
# Description: create MySQL 'BLAST_results' table
#***************************************************************************
sub create_blast_results_table {

	my ($self, $dbh) = @_;

	# BLAST results table 
	my $blast_results = "CREATE TABLE `BLAST_results` (
	  `Record_ID`     int(11) NOT NULL auto_increment,
	  `Probe_name`    varchar(100) NOT NULL default '0',
	  `Probe_gene`    varchar(100) NOT NULL default '0',
	  `Probe_type`    varchar(100) NOT NULL default '0',
	  `Organism`      varchar(100) NOT NULL default '0',
	  `Chunk_name`    varchar(100) NOT NULL default '0',
	  `Scaffold`      varchar(100) default 'NULL',
	  `Orientation`   varchar(100) NOT NULL default '0',
	  `Bit_score`     float NOT NULL default '0',
	  `Identity`      varchar(100) NOT NULL default '',
	  `e_value_num`   float  NOT NULL default '0',
	  `e_value_exp`   int(11)  NOT NULL default '0',
	  `Subject_start` int(11) NOT NULL default '0',
	  `Subject_end`   int(11) NOT NULL default '0',
	  `Query_start`   int(11) NOT NULL default '0',
	  `Query_end`     int(11) NOT NULL default '0',
	  `Align_len`     int(11) NOT NULL default '0',
	  `Gap_openings`  varchar(100) NOT NULL default '',
	  `Mismatches`    int(11) NOT NULL default '0',
	  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($blast_results);
	unless ($sth->execute()) { print "\n\t$blast_results\n\n\n"; exit;}

}

#***************************************************************************
# Subroutine:  create_extracted_table
# Description: create MySQL 'Extracted' table
#***************************************************************************
sub create_extracted_table {

	my ($self, $dbh) = @_;

	# Extracted sequences table 
	my $extracted = "CREATE TABLE `Extracted` (
	  `Record_ID`        int(11) NOT NULL auto_increment,
	  `BLAST_ID`         int(11) NOT NULL default '0',
	  `Probe_name`       varchar(100) NOT NULL default '0',
	  `Probe_gene`       varchar(100) NOT NULL default '0',
	  `Probe_type`       varchar(100) NOT NULL default '0',
	  `Organism`         varchar(100) NOT NULL default '0',
	  `Chunk_name`       varchar(100) default 'NULL',
	  `Scaffold`         varchar(100) default 'NULL',
	  `Extract_start`    int(11) NOT NULL default '0',
	  `Extract_end`      int(11) NOT NULL default '0',
	  `Sequence_length`  int(11) NOT NULL default '0',
	  `Sequence`         text NOT NULL,
	  `Assigned_to`      varchar(100) NOT NULL default '0',
	  `Assigned_to_gene` varchar(100) NOT NULL default '0',
	  `Orientation`      varchar(100) NOT NULL default '0',
	  `Bit_score`        float   NOT NULL default '0',
	  `Identity`         float   NOT NULL default '0',
	  `e_value_num`      float   NOT NULL default '0',
	  `e_value_exp`      int(11) NOT NULL default '0',
	  `Subject_start`    int(11) NOT NULL default '0',
	  `Subject_end`      int(11) NOT NULL default '0',
	  `Query_start`      int(11) NOT NULL default '0',
	  `Query_end`        int(11) NOT NULL default '0',
	  `RealEx_start`     int(11) NOT NULL default '0',
      `RealEx_end`       int(11) NOT NULL default '0',
	  `Align_len`        int(11) NOT NULL default '0',
	  `Gap_openings`     int(11) NOT NULL default '0',
	  `Mismatches`       int(11) NOT NULL default '0',
	  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($extracted);
	unless ($sth->execute()) { print "\n\t$extracted\n\n\n"; exit;}
}

#***************************************************************************
# Subroutine:  create_status_table
# Description: create MySQL 'Status' table
#***************************************************************************
sub create_status_table {

	my ($self, $dbh) = @_;

	# Status table (which BLAST queries have been executed)
	my $status = "CREATE TABLE `Status` (
	  `Record_ID`  int(11) NOT NULL auto_increment,
	  `Probe_ID`   varchar(100) NOT NULL default '',
	  `Probe_name` varchar(100) NOT NULL default '',
	  `Probe_gene` varchar(100) NOT NULL default '',
	  `Organism`   varchar(100) NOT NULL default '0',
	  `Version`    varchar(100) NOT NULL default '0',
	  `Chunk_name` varchar(100) NOT NULL default '0',
	  `Timestamp`  timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1";
	my $sth = $dbh->prepare($status);
	unless ($sth->execute()) { print "\n\t$status\n\n\n"; exit;}
}

############################################################################
# UPDATING DATABASES AND TABLES
############################################################################

#***************************************************************************
# Subroutine:  index_previously_executed_queries 
# Description: 
#***************************************************************************
sub index_previously_executed_queries {
	
	my ($self, $done_ref) = @_;
	
	my $status_table = $self->{status_table};
	my @data;
	my @fields = qw [ organism version chunk_name probe_name probe_gene ];
	$status_table->select_rows(\@fields, \@data);
	foreach my $data_ref (@data) {
		my $target_name = $data_ref->{chunk_name};
		my $version    = $data_ref->{version};
		my $organism   = $data_ref->{organism};
		my $probe_name = $data_ref->{probe_name};
		my $probe_gene = $data_ref->{probe_gene};
		#$devtools->print_hash($data_ref); die; # DEBUG
		unless ( $organism and $target_name and $version and $probe_name and $probe_gene) { 
		print "\n\t $organism and $target_name and $version and $probe_name and $probe_gene"; 
			die; 
		};
		my @key = ( $organism,   $version, $target_name, $probe_name, $probe_gene );
		my $key = join ('_', @key);
		#print "\n\t $key\n\n "; die;
		$done_ref->{$key} = 1;		
	}
	# DEBUG
	#$devtools->print_hash($done_ref);
}

############################################################################
# DROPPING AND FLUSHING DBs
############################################################################

#***************************************************************************
# Subroutine:  drop_screening_db
# Description: drop a screening database from the server
#***************************************************************************
sub drop_screening_db {

	my ($self) = @_;

	# Get connection variables from self
	my $server   = $self->{server};
	my $username = $self->{username};
	my $password = $self->{password};

	my $dbh = DBI->connect("dbi:mysql:ScreeningDB:$server", $username, $password);
	
	# Get database name	
	my $db_name = $self->{db_name};
	
	# Ask to make sure
	my $question = "\n\n\t Are you sure you want to DROP the $db_name database?";
	my $answer1 = $console->ask_yes_no_question($question);
	if ($answer1 eq 'y') {
		print "\n\t- - - PAUSING 5 seconds to allow for cancel - - -\n";
		sleep 5;
		
		my $dbh = DBI->connect("dbi:mysql:$db_name:$server", $username, $password);
		unless ($dbh) {	die "\n\t # Couldn't connect to $db_name database\n\n"; }
		
		# Execute the query
		my $drop = "DROP database  $db_name ";
		my $sth = $dbh->prepare($drop);
    	unless ($sth->execute()) { print $drop; exit;}	
	
		my $dbh2 = DBI->connect("dbi:mysql:ScreeningDB:$server", $username, $password);
		unless ($dbh2) {	die "\n\t # Couldn't connect to ScreeningDB database\n\n"; }
		my $cleanup = "DELETE FROM DB_detail WHERE DB_name = '$db_name'";	
		my $sth2 = $dbh2->prepare($cleanup);
    	unless ($sth2->execute()) { 
			print "\n\n\t Error: DB drop failed\n\n\n";
		}
	}
}

#***************************************************************************
# Subroutine:  flush_screening_db 
# Description: flush screening results database tables
#***************************************************************************
sub flush_screening_db {

	my ($self) = @_;
	
	# Ask to make sure
	my $db_name = $self->{db_name};
	my $question = "\n\n\t Are you sure you want to flush data in the $db_name database?";
	my $answer1 = $console->ask_yes_no_question($question);
	if ($answer1 eq 'y') {
		print "\n\t- - - PAUSING 3 seconds to allow for cancel - - -\n";
		sleep 3;

		# get tables
		my $status_table         = $self->{status_table};
		my $blast_results_table  = $self->{blast_results_table};
		my $extracted_table      = $self->{extracted_table};
		
		# Flush result tables
		$status_table->flush();
		$status_table->reset_primary_keys();
		$blast_results_table->flush();
		$blast_results_table->reset_primary_keys();
		$extracted_table->flush();
		$extracted_table->reset_primary_keys();
	}
}

############################################################################
# EOF
############################################################################
