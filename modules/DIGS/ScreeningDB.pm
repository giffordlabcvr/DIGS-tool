#!/usr/bin/perl -w
############################################################################
# Module:      ScreeningDB.pm
# Description: The DB module for DIGS databases
#              Contains routines for creating, and interacting with 
#              DBs in the DIGS framework
# History:     June 2011: Created by Robert Gifford 
############################################################################
package ScreeningDB;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;
use DBI;
############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::FileIO;
use Base::Console;

# Components
use DIGS::TargetDB;

# Database component modules
use Interface::MySQLtable;

############################################################################
# Globals
############################################################################

# Create base objects
my $fileio    = FileIO->new();
my $console   = Console->new();
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Create a new ScreeningDB.pm 'object'
#***************************************************************************
sub new {

	my ($invocant, $parameters) = @_;

	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {
		
		# Paths
		db_name              => undef,
		
		# DB connection variables
		server                => $parameters->{mysql_server},
		username              => $parameters->{mysql_username},
		password              => $parameters->{mysql_password},
		
		# Screening database core tables
		searches_table        => 0,
		blast_results_table   => 0,
		extracted_table       => 0,
		blast_chains_table    => 0,	
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# LOAD DATABASE & DATABASE TABLES
############################################################################

#***************************************************************************
# Subroutine:  load_screening_db 
# Description: load a screening database  
#***************************************************************************
sub load_screening_db {

	my ($self, $db_name) = @_;
	
	unless ($db_name) { die "\n\t DB name is not defined\n\n\n"; }

	# Get connection variables from self
	my $server   = $self->{server};
	my $username = $self->{username};
	my $password = $self->{password};
	unless ($server)   { die; }
	unless ($username) { die; }
	unless ($password) { die; }
	unless ($db_name)  { die; }

	# Set name
	$self->{db_name} = $db_name;
   
	# Load tables from the screening result DB we've been given
	my $dbh = DBI->connect("dbi:mysql:$db_name:$server", $username, $password);
	unless ($dbh) { die "\n\t Failed to connect to database\n\n\n"; }
	$self->{dbh} = $dbh;

	# Main Screening DB tables
	print   "\t  Connecting to DB:  $db_name";
	$self->load_searches_table($dbh);	
	$self->load_blast_results_table($dbh);	
	$self->load_extracted_table($dbh);	
	$self->load_blast_chains_table($dbh);	

}

############################################################################
# LOAD SCREENING DATABASE TABLES
############################################################################

#***************************************************************************
# Subroutine:  load_searches_table
# Description: load screening database table 'Searches_performed'
#***************************************************************************
sub load_searches_table {

	my ($self, $dbh) = @_;

	# Definition of the table
	my %searches_fields = (
		probe_id       => 'varchar',
		probe_name     => 'varchar',
		probe_gene     => 'varchar',
		genome_id      => 'varchar',
		organism       => 'varchar',
		data_type      => 'varchar',
		version        => 'varchar',
		target_name    => 'varchar',
	);
	my $searches_table = MySQLtable->new('Searches_performed', $dbh, \%searches_fields);
	$self->{searches_table} = $searches_table;
}

#***************************************************************************
# Subroutine:  load_blast_results_table
# Description: load screening database table 'BLAST_results'
#***************************************************************************
sub load_blast_results_table {

	my ($self, $dbh) = @_;
	
	# Definition of the table
	my %blast_fields = (
		extract_id       => 'int',
		probe_name       => 'varchar',
		probe_gene       => 'varchar',
		probe_type       => 'varchar',
		organism         => 'varchar',
		data_type        => 'varchar',
		version          => 'varchar',
		target_name      => 'varchar',
		scaffold         => 'varchar',
		orientation      => 'varchar',
		subject_start    => 'int',
		subject_end      => 'int',
		query_start      => 'int',
		query_end        => 'int',
		hit_length       => 'int',
		bit_score        => 'float',
		identity         => 'varchar',
		e_value_num      => 'float',
		e_value_exp      => 'int',
	  	subject_start    => 'int',
	  	subject_end      => 'int',
		query_start      => 'int',
	  	query_end        => 'int',
		gap_openings     => 'int',
		mismatches       => 'int',
	);
	my $blast_table = MySQLtable->new('BLAST_results', $dbh, \%blast_fields);
	$self->{blast_results_table} = $blast_table;
}

#***************************************************************************
# Subroutine:  load_extracted_table
# Description: load screening database table 'Extracted_sequences'
#***************************************************************************
sub load_extracted_table {

	my ($self, $dbh) = @_;

	# Definition of the table
	my %extract_fields = (
		organism         => 'varchar',
		version          => 'varchar',
		data_type        => 'varchar',
		target_name      => 'varchar',
		probe_type       => 'varchar',
		assigned_name    => 'varchar',
		assigned_gene    => 'varchar',
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
		query_start      => 'int',
	  	query_end        => 'int',
		align_len        => 'int',
		gap_openings     => 'int',
		mismatches       => 'int',
		sequence_length  => 'int',
		sequence         => 'text',
	);
	my $extract_table = MySQLtable->new('Extracted_sequences', $dbh, \%extract_fields);
	$self->{extracted_table} = $extract_table;
}

#***************************************************************************
# Subroutine:  load_blast_chains_table
# Description: load screening database table 'Extracted_sequences'
#***************************************************************************
sub load_blast_chains_table {

	my ($self, $dbh) = @_;

	# Definition of the table
	my %extract_fields = (
		extract_id       => 'int',
		probe_name       => 'varchar',
		probe_gene       => 'varchar',
		probe_type       => 'varchar',
		organism         => 'varchar',
		data_type        => 'varchar',
		version          => 'varchar',
		target_name      => 'varchar',
		scaffold         => 'varchar',
		orientation      => 'varchar',
		subject_start    => 'int',
		subject_end      => 'int',
		query_start      => 'int',
		query_end        => 'int',
		hit_length       => 'int',
		bit_score        => 'float',
		identity         => 'varchar',
		e_value_num      => 'float',
		e_value_exp      => 'int',
	  	subject_start    => 'int',
	  	subject_end      => 'int',
		query_start      => 'int',
	  	query_end        => 'int',
		gap_openings     => 'int',
		mismatches       => 'int',

	);
	my $extract_table = MySQLtable->new('BLAST_chains', $dbh, \%extract_fields);
	$self->{blast_chains_table} = $extract_table;
}

############################################################################
# CREATE SCREENING DATABASE TABLES
############################################################################

#***************************************************************************
# Subroutine:  create_screening_db
# Description: Create the tables of a DIGS screening DB
#***************************************************************************
sub create_screening_db {

	my ($self, $db_name) = @_;

	# Get connection variables from self
	my $server   = $self->{server};
	my $username = $self->{username};
	my $password = $self->{password};

	# CREATE THE DB 
	print "\n\t ### Creating '$db_name' screening database\n";
    my $drh = DBI->install_driver("mysql");
    my $rc = $drh->func("createdb", $db_name, $server, $username, $password, 'admin');
    my $dbh = DBI->connect("dbi:mysql:$db_name:$server", $username, $password);
	unless ($dbh) {
		die "\n\t # Couldn't create screening database '$db_name'\n\n";
	}	
	$self->create_searches_table($dbh);
	$self->create_blast_results_table($dbh);
	$self->create_extracted_table($dbh);
	$self->create_blast_chains_table($dbh);
}

#***************************************************************************
# Subroutine:  create_searches_table
# Description: create MySQL 'Searches_performed' table
#***************************************************************************
sub create_searches_table {

	my ($self, $dbh) = @_;

	# Searches_performed table (which BLAST queries have been executed)
	my $searches = "CREATE TABLE `Searches_performed` (
	  `Record_ID`   int(11) NOT NULL auto_increment,
	  `Probe_ID`    varchar(100) NOT NULL default '',
	  `Probe_name`  varchar(100) NOT NULL default '',
	  `Probe_gene`  varchar(100) NOT NULL default '',
	  `Genome_ID`   varchar(100) NOT NULL default '',
	  `Organism`    varchar(100) NOT NULL default '',
	  `Data_type`   varchar(100) NOT NULL default '',
	  `Version`     varchar(100) NOT NULL default '',
	  `Target_name` varchar(100) NOT NULL default '',
	  `Timestamp`   timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1";
	my $sth = $dbh->prepare($searches);
	unless ($sth->execute()) { print "\n\t$searches\n\n\n"; exit;}
}

#***************************************************************************
# Subroutine:  create_blast_results_table
# Description: create MySQL 'BLAST_results' table
#***************************************************************************
sub create_blast_results_table {

	my ($self, $dbh) = @_;

	# BLAST results table 
	my $blast_results = "CREATE TABLE `BLAST_results` (
	  `Record_ID`     int(11) NOT NULL auto_increment,
	  `Extract_ID`    int(11) default '0',

	  `Probe_name`    varchar(100) NOT NULL default '0',
	  `Probe_gene`    varchar(100) NOT NULL default '0',
	  `Probe_type`    varchar(100) NOT NULL default '0',
	  `Organism`      varchar(100) NOT NULL default '0',
	  `Data_type`     varchar(100) NOT NULL default '0',
	  `Version`       varchar(100) NOT NULL default '0',
	  `Target_name`   varchar(100) NOT NULL default '0',

	  `Scaffold`      varchar(100) NOT NULL default '0',
	  `Orientation`   varchar(100) NOT NULL default '0',
	  `Subject_start` int(11) NOT NULL default '0',
	  `Subject_end`   int(11) NOT NULL default '0',
	  `Query_start`   int(11) NOT NULL default '0',
	  `Query_end`     int(11) NOT NULL default '0',
	  `Hit_length`    int(11) NOT NULL default '0',

	  `Bit_score`        float   NOT NULL default '0',
	  `Identity`         float   NOT NULL default '0',
	  `e_value_num`      float   NOT NULL default '0',
	  `e_value_exp`      int(11) NOT NULL default '0',
	  `Gap_openings`     int(11) NOT NULL default '0',
	  `Mismatches`       int(11) NOT NULL default '0',

	  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($blast_results);
	unless ($sth->execute()) { print "\n\t$blast_results\n\n\n"; exit;}
}

#***************************************************************************
# Subroutine:  create_extracted_table
# Description: create MySQL 'Extracted_sequences' table
#***************************************************************************
sub create_extracted_table {

	my ($self, $dbh) = @_;

	# Extracted_sequences sequences table 
	my $extracted = "CREATE TABLE `Extracted_sequences` (
	  `Record_ID`        int(11) NOT NULL auto_increment,

	  `Organism`         varchar(100) NOT NULL default '0',
	  `Data_type`        varchar(100) NOT NULL default '0',
	  `Version`          varchar(100) NOT NULL default '0',
	  `Target_name`      varchar(100) NOT NULL default '0',
	  `Probe_type`       varchar(100) NOT NULL default '0',
	  `Scaffold`         varchar(100) NOT NULL default '0',
	  `Extract_start`    int(11) NOT NULL default '0',
	  `Extract_end`      int(11) NOT NULL default '0',
	  `Sequence_length`  int(11) NOT NULL default '0',
	  `Sequence`         text NOT NULL,
	  `Assigned_name`    varchar(100) NOT NULL default '0',
	  `Assigned_gene`    varchar(100) NOT NULL default '0',
	  `Orientation`      varchar(100) NOT NULL default '0',
	  `Bit_score`        float   NOT NULL default '0',
	  `Identity`         float   NOT NULL default '0',
	  `e_value_num`      float   NOT NULL default '0',
	  `e_value_exp`      int(11) NOT NULL default '0',
	  `Subject_start`    int(11) NOT NULL default '0',
	  `Subject_end`      int(11) NOT NULL default '0',
	  `Query_start`      int(11) NOT NULL default '0',
	  `Query_end`        int(11) NOT NULL default '0',
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
# Subroutine:  create_blast_chains_table
# Description: create MySQL 'BLAST_results' table
#***************************************************************************
sub create_blast_chains_table {

	my ($self, $dbh) = @_;

	# BLAST results table 
	my $blast_chains = "CREATE TABLE `BLAST_chains` (
	  `Record_ID`     int(11) NOT NULL auto_increment,
	  `Extract_ID`    int(11) NOT NULL default '0',

	  `Probe_name`    varchar(100) NOT NULL default '0',
	  `Probe_gene`    varchar(100) NOT NULL default '0',
	  `Probe_type`    varchar(100) NOT NULL default '0',
	  `Organism`      varchar(100) NOT NULL default '0',
	  `Data_type`     varchar(100) NOT NULL default '0',
	  `Version`       varchar(100) NOT NULL default '0',
	  `Target_name`   varchar(100) NOT NULL default '0',

	  `Scaffold`      varchar(100) default 'NULL',
	  `Orientation`   varchar(100) NOT NULL default '0',
	  `Subject_start` int(11) NOT NULL default '0',
	  `Subject_end`   int(11) NOT NULL default '0',
	  `Query_start`   int(11) NOT NULL default '0',
	  `Query_end`     int(11) NOT NULL default '0',
	  `Hit_length`    int(11) NOT NULL default '0',

	  `Bit_score`        float   NOT NULL default '0',
	  `Identity`         float   NOT NULL default '0',
	  `e_value_num`      float   NOT NULL default '0',
	  `e_value_exp`      int(11) NOT NULL default '0',
	  `Gap_openings`     int(11) NOT NULL default '0',
	  `Mismatches`       int(11) NOT NULL default '0',

	  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($blast_chains);
	unless ($sth->execute()) { print "\n\t$blast_chains\n\n\n"; exit;}
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

	# Get database name	
	my $db_name = $self->{db_name};
	
	# Confirm
	my $question = "\n\n\t Are you sure you want to DROP the $db_name database?";
	my $answer1 = $console->ask_yes_no_question($question);
	if ($answer1 eq 'y') {
		my $dbh = DBI->connect("dbi:mysql:$db_name:$server", $username, $password);
		unless ($dbh) {	die "\n\t # Couldn't connect to $db_name database\n\n"; }
		# Execute the query
		my $drop = "DROP database  $db_name ";
		my $sth = $dbh->prepare($drop);
    	unless ($sth->execute()) { print $drop; exit;}	
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
		#sleep 3;

		# get tables
		my $blast_results_table  = $self->{blast_results_table};
		my $extracted_table      = $self->{extracted_table};
		my $searches_table       = $self->{searches_table};
		my $blast_chains_table   = $self->{blast_chains_table};
		
		# Flush result tables
		$blast_results_table->flush();
		$blast_results_table->reset_primary_keys();
		$extracted_table->flush();
		$extracted_table->reset_primary_keys();
		$searches_table->flush();
		$searches_table->reset_primary_keys();
		$blast_chains_table->flush();
		$blast_chains_table->reset_primary_keys();

	}
}

############################################################################
# EXTEND SCREENING DB FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  create_ancillary_table
# Description: create an ancillary table for the screening DB
#***************************************************************************
sub create_ancillary_table {

	my ($self, $table_name, $fields_array_ref, $fields_hash_ref) = @_;
	
	my $dbh = $self->{dbh};
	unless ($dbh) { die "\n\t Couldn't retrieve database handle \n\n"; }
	my $db_name = $self->{db_name};
	print "\n\t Creating ancillary table '$table_name' in $db_name screening database";

	# Remove illegal characters from table name
	unless ($table_name) { die "\n\t No table name defined \n\n"; }
	$table_name =~ s/\s+//g;

	# Make the create statement
    my $create_lead = "CREATE TABLE IF NOT EXISTS `$table_name` (
      `Record_ID` int(11) NOT NULL auto_increment, ";

	my $create_middle = '';
	my $i= 0;
	foreach my $field (@$fields_array_ref) {
		$i++;
		my $type = $fields_hash_ref->{$i};
		unless ($type) {
			$type = "varchar(100) NOT NULL default '0'";
        }
		elsif ($type eq 'varchar') {
			$type = "varchar(100) NOT NULL default '0'";
		}
		$create_middle .= " `$field` $type, ";
	}

	my $create_tail = "`Timestamp` timestamp NOT NULL default
      CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
      PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";

	my $statement = $create_lead . $create_middle . $create_tail;
	my $sth = $dbh->prepare($statement);
	unless ($sth->execute()) { print "\n\t$statement\n\n\n"; exit; }

	return $table_name;
}

#***************************************************************************
# Subroutine:  get_ancillary_table_names
# Description: get the names of all ancillary tables in a screening DB
#***************************************************************************
sub get_ancillary_table_names {

	my ($self, $anc_tables_ref) = @_;
	
	my $dbh = $self->{dbh};
	unless ($dbh) { die "\n\t Couldn't retrieve database handle \n\n"; }
	my $show = "SHOW TABLES ";
	my $sth = $dbh->prepare($show);
   	unless ($sth->execute()) { print $show; exit;}	
	
	my $i = 0;
	while (my $row = $sth->fetchrow_arrayref) {
		foreach my $item (@$row) {
			chomp $item;
			$i++;
			if ($item eq 'Searches_performed')      { next; }
			elsif ($item eq 'BLAST_results')        { next; }
			elsif ($item eq 'Extracted_sequences')  { next; }
			else {
				push (@$anc_tables_ref, $item)
			}
		}
	}
}

#***************************************************************************
# Subroutine:  drop_ancillary_table
# Description: drop an ancillary table from the screening database 
#***************************************************************************
sub drop_ancillary_table {

	my ($self, $table_name) = @_;

	# Execute the query
	my $dbh = $self->{dbh};
	unless ($dbh) { die "\n\t Couldn't retrieve database handle \n\n"; }
	my $drop = "DROP TABLE IF EXISTS $table_name ";
	my $sth = $dbh->prepare($drop);
   	unless ($sth->execute()) { print $drop; exit;}	

}
	
############################################################################
# EOF
############################################################################
