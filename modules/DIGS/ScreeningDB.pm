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
my $devtools  = DevTools->new();
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
		active_set_table      => 0,
		digs_results_table    => 0,
		blast_chains_table    => 0,	
		loci_table            => 0,	
		loci_chains_table     => 0,
		
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
	$self->load_active_set_table($dbh);	
	$self->load_digs_results_table($dbh, 'digs_results');	
	$self->load_blast_chains_table($dbh);	

}

#***************************************************************************
# Subroutine:  load_searches_table
# Description: load screening database table 'searches_performed'
#***************************************************************************
sub load_searches_table {

	my ($self, $dbh) = @_;

	# Definition of the table
	my %searches_fields = (
		probe_id         => 'varchar',
		probe_name       => 'varchar',
		probe_gene       => 'varchar',
		target_id        => 'varchar',
		organism  => 'varchar',
		target_datatype  => 'varchar',
		target_version   => 'varchar',
		target_name      => 'varchar',
	);
	my $searches_table = MySQLtable->new('searches_performed', $dbh, \%searches_fields);
	$self->{searches_table} = $searches_table;
}

#***************************************************************************
# Subroutine:  load_active_set_table
# Description: load screening database table 'active_set'
#***************************************************************************
sub load_active_set_table {

	my ($self, $dbh) = @_;
	
	# Definition of the table
	my %active_set_fields = (
	
		digs_result_id       => 'int',
		probe_name       => 'varchar',
		probe_gene       => 'varchar',
		probe_type       => 'varchar',
		organism  => 'varchar',
		target_datatype  => 'varchar',
		target_version   => 'varchar',
		target_name      => 'varchar',
		scaffold         => 'varchar',
		orientation      => 'varchar',
		subject_start    => 'int',
		subject_end      => 'int',
		query_start      => 'int',
		query_end        => 'int',
		align_len        => 'int',
		bitscore         => 'float',
		identity         => 'varchar',
		evalue_num       => 'float',
		evalue_exp       => 'int',
	  	subject_start    => 'int',
	  	subject_end      => 'int',
		query_start      => 'int',
	  	query_end        => 'int',
		gap_openings     => 'int',
		mismatches       => 'int',
	);
	my $blast_table = MySQLtable->new('active_set', $dbh, \%active_set_fields);
	$self->{active_set_table} = $blast_table;
}

#***************************************************************************
# Subroutine:  load_digs_results_table
# Description: load screening database table 'digs_results'
#***************************************************************************
sub load_digs_results_table {

	my ($self, $dbh, $name) = @_;

	# Definition of the table
	my %digs_fields = (
	
		organism         => 'varchar',
		target_version   => 'varchar',
		target_datatype  => 'varchar',
		target_name      => 'varchar',
		probe_type       => 'varchar',
		assigned_name    => 'varchar',
		assigned_gene    => 'varchar',
		scaffold         => 'varchar',
		extract_start    => 'varchar',
		extract_end      => 'varchar',
		orientation      => 'varchar',
		bitscore         => 'float',
		identity         => 'varchar',
		evalue_num       => 'float',
		evalue_exp       => 'int',
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
	my $digs_table = MySQLtable->new($name, $dbh, \%digs_fields);
	my $table_name = $name . '_table';
	$self->{$table_name} = $digs_table;
}

#***************************************************************************
# Subroutine:  load_blast_chains_table
# Description: load screening database table 'digs_results'
#***************************************************************************
sub load_blast_chains_table {

	my ($self, $dbh) = @_;

	# Definition of the table
	my %extract_fields = (
	
		digs_result_id   => 'int',
		probe_name       => 'varchar',
		probe_gene       => 'varchar',
		probe_type       => 'varchar',
		organism         => 'varchar',
		target_datatype  => 'varchar',
		target_version   => 'varchar',
		target_name      => 'varchar',
		scaffold         => 'varchar',
		orientation      => 'varchar',
		subject_start    => 'int',
		subject_end      => 'int',
		query_start      => 'int',
		query_end        => 'int',
		align_len        => 'int',
		bitscore         => 'float',
		identity         => 'varchar',
		evalue_num       => 'float',
		evalue_exp       => 'int',
	  	subject_start    => 'int',
	  	subject_end      => 'int',
		query_start      => 'int',
	  	query_end        => 'int',
		gap_openings     => 'int',
		mismatches       => 'int',

	);
	my $extract_table = MySQLtable->new('blast_chains', $dbh, \%extract_fields);
	$self->{blast_chains_table} = $extract_table;
}

#***************************************************************************
# Subroutine:  load_loci_table
# Description: load screening database table 'loci'
#***************************************************************************
sub load_loci_table {

    my ($self, $dbh) = @_;

    # Definition of the loci table
    my %loci_fields = (

        organism         => 'varchar',
        target_version   => 'varchar',
        target_datatype  => 'varchar',
        target_name      => 'varchar',
        scaffold         => 'varchar',
        orientation      => 'varchar',
        assigned_name    => 'varchar',
        extract_start    => 'int',
        extract_end      => 'int',
        locus_structure  => 'text',

    );   
    my $loci_table = MySQLtable->new('loci', $dbh, \%loci_fields);
    $self->{loci_table} = $loci_table;

}

#***************************************************************************
# Subroutine:  load_loci_chains_table
# Description: load screening database table 'loci_chains'
#***************************************************************************
sub load_loci_chains_table {

    my ($self, $dbh) = @_;

    # Definition of the loci table
    my %loci_fields = (
        locus_id       => 'int',
        digs_result_id => 'int',
    );   
    my $loci_table = MySQLtable->new('loci_chains', $dbh, \%loci_fields);
    $self->{loci_chains_table} = $loci_table;
}

#***************************************************************************
# Subroutine:  load_nomenclature_table
# Description: load screening database table 'nomenclature'
#***************************************************************************
sub load_nomenclature_table {

	my ($self, $dbh) = @_;

	# Definition of the table
	my %nomenclature_fields = (	
        track_name       => 'varchar',
		assigned_name    => 'varchar',
		scaffold         => 'varchar',
		assigned_gene    => 'varchar',
		extract_start    => 'varchar',
		extract_end      => 'varchar',
		orientation      => 'varchar',
		sequence_length  => 'int',
		namespace_id     => 'varchar',	
	);
	
	my $nomenclature_table = MySQLtable->new('nomenclature', $dbh, \%nomenclature_fields);
	$self->{nomenclature_table} = $nomenclature_table;
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
	$self->create_active_set_table($dbh);
	$self->create_digs_results_table($dbh);
	$self->create_blast_chains_table($dbh);
}

#***************************************************************************
# Subroutine:  create_searches_table
# Description: create MySQL 'searches_performed' table
#***************************************************************************
sub create_searches_table {

	my ($self, $dbh) = @_;

	# searches_performed table (which BLAST queries have been executed)
	my $searches = "CREATE TABLE `searches_performed` (
	  `record_ID`        int(11) NOT NULL auto_increment,
	  `probe_ID`         varchar(100) NOT NULL default '',
	  `probe_name`       varchar(100) NOT NULL default '',
	  `probe_gene`       varchar(100) NOT NULL default '',
	  `target_id`        varchar(100) NOT NULL default '',
	  `organism`         varchar(100) NOT NULL default '0',
	  `target_datatype`  varchar(100) NOT NULL default '0',
	  `target_version`   varchar(100) NOT NULL default '0',
	  `target_name`      varchar(100) NOT NULL default '0',
	  `timestamp`   timestamp NOT NULL default CURRENT_TIMESTAMP 
	                on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`record_id`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1";
	my $sth = $dbh->prepare($searches);
	unless ($sth->execute()) { print "\n\t$searches\n\n\n"; exit;}
}

#***************************************************************************
# Subroutine:  create_active_set_table
# Description: create MySQL 'active_set' table
#***************************************************************************
sub create_active_set_table {

	my ($self, $dbh) = @_;

	# Active result set table 
	my $active_set = "CREATE TABLE `active_set` (
	  `record_ID`       int(11) NOT NULL auto_increment,
	  `digs_result_id`      int(11) default '0',

	  `probe_name`       varchar(100) NOT NULL default '0',
	  `probe_gene`       varchar(100) NOT NULL default '0',
	  `probe_type`       varchar(100) NOT NULL default '0',
	  `organism`         varchar(100) NOT NULL default '0',
	  `target_datatype`  varchar(100) NOT NULL default '0',
	  `target_version`   varchar(100) NOT NULL default '0',
	  `target_name`      varchar(100) NOT NULL default '0',

	  `scaffold`         varchar(100) NOT NULL default '0',
	  `orientation`      varchar(100) NOT NULL default '0',
	  `subject_start`    int(11) NOT NULL default '0',
	  `subject_end`      int(11) NOT NULL default '0',
	  `query_start`      int(11) NOT NULL default '0',
	  `query_end`        int(11) NOT NULL default '0',
	  `align_len`       int(11) NOT NULL default '0',

	  `bitscore`         float   NOT NULL default '0',
	  `identity`         float   NOT NULL default '0',
	  `evalue_num`       float   NOT NULL default '0',
	  `evalue_exp`       int(11) NOT NULL default '0',
	  `gap_openings`     int(11) NOT NULL default '0',
	  `mismatches`       int(11) NOT NULL default '0',

	  `timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`record_id`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($active_set);
	unless ($sth->execute()) { print "\n\t$active_set\n\n\n"; exit;}
}

#***************************************************************************
# Subroutine:  create_digs_results_table
# Description: create MySQL 'digs_results' table
#***************************************************************************
sub create_digs_results_table {

	my ($self, $dbh) = @_;

	# digs_results table 
	my $digs_results = "CREATE TABLE `digs_results` (
	  `record_ID`        int(11) NOT NULL auto_increment,
	  `organism`         varchar(100) NOT NULL default '0',
	  `target_datatype`  varchar(100) NOT NULL default '0',
	  `target_version`   varchar(100) NOT NULL default '0',
	  `target_name`      varchar(100) NOT NULL default '0',
	  `probe_type`       varchar(100) NOT NULL default '0',
	  `scaffold`         varchar(100) NOT NULL default '0',
	  `extract_start`    int(11) NOT NULL default '0',
	  `extract_end`      int(11) NOT NULL default '0',
	  `sequence_length`  int(11) NOT NULL default '0',
	  `sequence`         text NOT NULL,
	  `assigned_name`    varchar(100) NOT NULL default '0',
	  `assigned_gene`    varchar(100) NOT NULL default '0',
	  `orientation`      varchar(100) NOT NULL default '0',
	  `bitscore`         float   NOT NULL default '0',
	  `identity`         float   NOT NULL default '0',
	  `evalue_num`       float   NOT NULL default '0',
	  `evalue_exp`       int(11) NOT NULL default '0',
	  `subject_start`    int(11) NOT NULL default '0',
	  `subject_end`      int(11) NOT NULL default '0',
	  `query_start`      int(11) NOT NULL default '0',
	  `query_end`        int(11) NOT NULL default '0',
	  `align_len`        int(11) NOT NULL default '0',
	  `gap_openings`     int(11) NOT NULL default '0',
	  `mismatches`       int(11) NOT NULL default '0',

	  `timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`record_id`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($digs_results);
	unless ($sth->execute()) { print "\n\t$digs_results\n\n\n"; exit;}
}

#***************************************************************************
# Subroutine:  create_blast_chains_table
# Description: create MySQL 'active_set' table
#***************************************************************************
sub create_blast_chains_table {

	my ($self, $dbh) = @_;

	# Active result set table 
	my $blast_chains = "CREATE TABLE `blast_chains` (
	  `record_ID`       int(11) NOT NULL auto_increment,
	  `digs_result_id`      int(11) NOT NULL default '0',

	  `probe_name`      varchar(100) NOT NULL default '0',
	  `probe_gene`      varchar(100) NOT NULL default '0',
	  `probe_type`      varchar(100) NOT NULL default '0',
	  `organism`        varchar(100) NOT NULL default '0',
	  `target_datatype` varchar(100) NOT NULL default '0',
	  `target_version`  varchar(100) NOT NULL default '0',
	  `target_name`     varchar(100) NOT NULL default '0',

	  `scaffold`        varchar(100) default 'NULL',
	  `orientation`     varchar(100) NOT NULL default '0',
	  `subject_start`   int(11) NOT NULL default '0',
	  `subject_end`     int(11) NOT NULL default '0',
	  `query_start`     int(11) NOT NULL default '0',
	  `query_end`       int(11) NOT NULL default '0',
	  `align_len`       int(11) NOT NULL default '0',

	  `bitscore`        float   NOT NULL default '0',
	  `identity`        float   NOT NULL default '0',
	  `evalue_num`      float   NOT NULL default '0',
	  `evalue_exp`      int(11) NOT NULL default '0',
	  `gap_openings`    int(11) NOT NULL default '0',
	  `mismatches`      int(11) NOT NULL default '0',

	  `timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`record_id`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($blast_chains);
	unless ($sth->execute()) { print "\n\t$blast_chains\n\n\n"; exit;}
}

#***************************************************************************
# Subroutine:  create_loci_table
# Description: create MySQL 'loci' table
#***************************************************************************
sub create_loci_table {

    my ($self, $dbh) = @_;

    # consolidated loci table 
    my $loci = "CREATE TABLE `loci` (
    
        `record_id`         int(11) NOT NULL auto_increment,
        `organism`          varchar(100) NOT NULL default '0',
        `target_datatype`   varchar(100) NOT NULL default '0',
        `target_version`    varchar(100) NOT NULL default '0',
        `target_name`       varchar(100) NOT NULL default '0',

        `scaffold`          varchar(100) default 'NULL',
        `orientation`       varchar(100) NOT NULL default '0',

        `assigned_name`     varchar(100) NOT NULL default '0',
        `extract_start`     int(11) NOT NULL default '0',
        `extract_end`       int(11) NOT NULL default '0',
        `locus_structure`   text NOT NULL,
      
        `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
      PRIMARY KEY  (`record_id`)
    ) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
    my $sth = $dbh->prepare($loci);
    unless ($sth->execute()) { print "\n\t$loci\n\n\n"; exit;}
}

#***************************************************************************
# Subroutine:  create_loci_chains_table
# Description: create MySQL 'loci_chains' table
#***************************************************************************
sub create_loci_chains_table {

    my ($self, $dbh) = @_;

    # consolidated loci table 
    my $loci = "CREATE TABLE `loci_chains`  (
    
        `record_id`       int(11) NOT NULL auto_increment,
        `locus_id`        int(11) NOT NULL default '0',
        `digs_result_id`  int(11) NOT NULL default '0',
     
        `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
      PRIMARY KEY  (`record_id`)
    ) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
    my $sth = $dbh->prepare($loci);
    unless ($sth->execute()) { print "\n\t$loci\n\n\n"; exit;}
}

#***************************************************************************
# Subroutine:  create_nomenclature_table
# Description: create MySQL 'nomenclature' table
#***************************************************************************
sub create_nomenclature_table {

	my ($self, $dbh) = @_;

	#  Nomenclature table 
	my $nomenclature = "CREATE TABLE `nomenclature` (
	  `record_id`        int(11) NOT NULL auto_increment,
	  `track_name`       varchar(100) NOT NULL default '0',
	  `assigned_name`    varchar(100) NOT NULL default '0',
	  `scaffold`         varchar(100) NOT NULL default '0',
	  `extract_start`    int(11) NOT NULL default '0',
	  `extract_end`      int(11) NOT NULL default '0',
	  `sequence_length`  int(11) NOT NULL default '0',
	  `assigned_gene`    varchar(100) NOT NULL default '0',
	  `orientation`      varchar(100) NOT NULL default '0',
	  `namespace_id`     varchar(100) NOT NULL default '0',
	  `timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`record_id`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($nomenclature);
	unless ($sth->execute()) { print "\n\t$nomenclature\n\n\n"; exit;}
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
	unless ($db_name) { die; }
	
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
		my $active_set_table  = $self->{active_set_table};
		my $digs_results_table      = $self->{digs_results_table};
		my $searches_table       = $self->{searches_table};
		my $blast_chains_table   = $self->{blast_chains_table};
		
		# Flush result tables
		$active_set_table->flush();
		$active_set_table->reset_primary_keys();
		$digs_results_table->flush();
		$digs_results_table->reset_primary_keys();
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

	my $create_tail = "`timestamp` timestamp NOT NULL default
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
	
	# Get DB object
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
			if ($item eq 'searches_performed')      { next; }
			elsif ($item eq 'active_set')        { next; }
			elsif ($item eq 'digs_results')  { next; }
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

#***************************************************************************
# Subroutine:  backup_digs_results_table 
# Description: 
#***************************************************************************
sub backup_digs_results_table {

    my ($self) = @_;

	my $copy_name;
	my $is_unique = undef;
	my $i = 0;
	do {
		$i++;
		$copy_name = 'digs_results_' . $i;
		my $exists = $self->does_table_exist($copy_name);
		unless ($exists) { $is_unique = 'true'; }
	} until ($is_unique);
	
	print "\n\t Copying DIGS results to $copy_name";
	my $copy_sql_1 = "CREATE TABLE $copy_name LIKE digs_results;";
	my $copy_sql_2 = "INSERT $copy_name SELECT * FROM digs_results";
    my $dbh = $self->{dbh};
	my $sth = $dbh->prepare($copy_sql_1);
   	unless ($sth->execute()) { print $copy_sql_1; exit;}	
   	#$sth = $dbh->prepare($copy_sql_2);
   	#unless ($sth->execute()) { print $copy_sql_2; exit;}
   	
   	return $copy_name;

}

#***************************************************************************
# Subroutine:  translate_schema
# Description: 
#***************************************************************************
sub translate_schema {

	my ($self) = @_;

	# Get connection variables from self
	my $server   = $self->{server};
	my $username = $self->{username};
	my $password = $self->{password};
	my $db_name  = $self->{db_name};

	# Get database handle
	my $dbh = DBI->connect("dbi:mysql:$db_name:$server", $username, $password);
	unless ($dbh) {	die "\n\t # Couldn't connect to $db_name database\n\n"; }

	# Create new tables if they don't already exist
	my $searches_exists = $self->does_table_exist('searches_performed');
	unless ($searches_exists) {
		$self->create_searches_table($dbh);
	}
	my $active_set_exists = $self->does_table_exist('active_set');
	unless ($active_set_exists) {
		$self->create_active_set_table($dbh);
	}
	my $digs_results_exists = $self->does_table_exist('digs_results');
	unless ($digs_results_exists) {
		$self->create_digs_results_table($dbh);
	}
	my $blast_chains_exists = $self->does_table_exist('blast_chains');
	unless ($blast_chains_exists) {
		$self->create_blast_chains_table($dbh);
	}

	# Translate 'Extracted' to 'digs_results'
	$self->load_extracted_table($dbh);
	my $extracted_table    = $self->{extracted_table};
	my $digs_results_table = $self->{digs_results_table};
	$digs_results_table->flush();
	my @extracted_rows;
	my @extracted_fields = qw [ organism version data_type target_name probe_type 
	                            assigned_name assigned_gene scaffold 
	                            extract_start extract_end orientation 
	                            bit_score e_value_num e_value_exp align_len 
	                            gap_openings mismatches sequence_length sequence ];
	$extracted_table->select_rows(\@extracted_fields, \@extracted_rows);
	foreach my $row_ref (@extracted_rows) {	
		$row_ref->{target_version}  = $row_ref->{version};
		$row_ref->{target_datatype} = $row_ref->{data_type};
		$row_ref->{bitscore}        = $row_ref->{bit_score};
		$row_ref->{evalue_num}      = $row_ref->{e_value_num};
		$row_ref->{evalue_exp}      = $row_ref->{e_value_exp};
		$digs_results_table->insert_row($row_ref);
	}

	# Translate 'Status' to 'searches_performed'
	$self->load_status_table($dbh);
	my $status_table   = $self->{status_table};
	my $searches_table = $self->{searches_table};
	$searches_table->flush();
	my @status_rows;
	my @status_fields = qw [ probe_id probe_name probe_gene 
	                         genome_id organism data_type version target_name ];
	$status_table->select_rows(\@status_fields, \@status_rows);
	foreach my $row_ref (@status_rows) {		
		$row_ref->{target_id}       = $row_ref->{target_id};
		$row_ref->{target_version}  = $row_ref->{version};
		$row_ref->{target_datatype} = $row_ref->{data_type};
		$searches_table->insert_row($row_ref);
	}

}

#***************************************************************************
# Subroutine:  load_status_table
# Description: load screening database table 'Status'
# Note: Deprecated table - this function required for schema translate 
#***************************************************************************
sub load_status_table {

	my ($self, $dbh) = @_;

	# Definition of the table
	my %status_fields = (
		probe_id       => 'varchar',
		probe_name     => 'varchar',
		probe_gene     => 'varchar',
		target_id      => 'varchar',
		organism       => 'varchar',
		data_type      => 'varchar',
		version        => 'varchar',
		target_name    => 'varchar',
	);
	my $status_table = MySQLtable->new('Status', $dbh, \%status_fields);
	$self->{status_table} = $status_table;
}

#***************************************************************************
# Subroutine:  load_extracted_table
# Description: load screening database table 'Extracted'
# Note: Deprecated table - this function required for schema translate 
#***************************************************************************
sub load_extracted_table {

	my ($self, $dbh) = @_;

	# Definition of the table
	my %extract_fields = (
		blast_id         => 'blast_id',
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
	my $extract_table = MySQLtable->new('Extracted', $dbh, \%extract_fields);
	$self->{extracted_table} = $extract_table;
}

#***************************************************************************
# Subroutine:  does_table_exist
# Description: check if a table exists in the screening DB
#***************************************************************************
sub does_table_exist {

    my ($self, $table_name) = @_;

	# Get connection variables from self
	my $server   = $self->{server};
	my $username = $self->{username};
	my $password = $self->{password};
	my $db_name  = $self->{db_name};

	unless ($server)   { die; }
	unless ($username) { die; }
	unless ($password) { die; }

	# Set name
	my $info_db_name = 'information_schema';
   
	# Load tables from the screening result DB we've been given
	my $dbh = DBI->connect("dbi:mysql:$info_db_name:$server", $username, $password);
	unless ($dbh) { die "\n\t Failed to connect to database\n\n\n"; }

	my $query = "SELECT * FROM information_schema.tables
                 WHERE table_schema = '$db_name' 
  			     AND table_name = '$table_name'
		         LIMIT 1;";
	my $sth = $dbh->prepare($query);
	unless ($sth->execute()) { print $query; exit; }
	my $row_count = 0;
	my $exists = undef;
	while (my $row = $sth->fetchrow_arrayref) {
		$exists = 'TRUE';		
	}

	return $exists;	
}

#***************************************************************************
# Subroutine:  does_db_exist
# Description: check if a screening database exists
#***************************************************************************
sub does_db_exist {

    my ($self, $db_name) = @_;

	unless ($db_name)  { die; }

	# Get connection variables from self
	my $server   = $self->{server};
	my $username = $self->{username};
	my $password = $self->{password};

	unless ($server)   { die; }
	unless ($username) { die; }
	unless ($password) { die; }

	# Set name
	my $info_db_name = 'information_schema';
   
	# Load tables from the screening result DB we've been given
	my $dbh = DBI->connect("dbi:mysql:$info_db_name:$server", $username, $password);
	unless ($dbh) { die "\n\t Failed to connect to database\n\n\n"; }
	my $query = "SELECT * FROM information_schema.tables
                 WHERE table_schema = '$db_name' 
		         LIMIT 1;";
	my $sth = $dbh->prepare($query);
	unless ($sth->execute()) { print $query; exit; }
	my $row_count = 0;
	my $exists = undef;
	while (my $row = $sth->fetchrow_arrayref) {
		$exists = 'TRUE';		
	}

	return $exists;	
}

############################################################################
# EOF
############################################################################
