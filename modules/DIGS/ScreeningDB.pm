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
# Description: Create a new ScreeningDB 'object'
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
		loci_table            => 0,	
		loci_chains_table     => 0,

		# Flags
		verbose                => $parameters->{verbose},
		force                  => $parameters->{force},
		
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# UPDATE SCREENING DB
############################################################################

#***************************************************************************
# Subroutine:  update_db
# Description: update the screening DB based on a completed round of DIGS
#***************************************************************************
sub update_db {

	my ($self, $to_delete_ref, $extracted_ref, $table_name) = @_;
		
	# Get parameters from self
	my $verbose             = $self->{verbose};
	my $digs_results_table  = $self->{$table_name}; 
	my $active_set_table    = $self->{active_set_table}; 
	# DEBUG $devtools->print_array($extracted_ref); die;

	# Delete superfluous data from the digs_results table
	my $deleted = '0';
	foreach my $id (@$to_delete_ref) {
		my $extracted_where = " WHERE record_id = $id ";	
		if ($verbose) { print "\n\t\t    - Deleting redundant locus '$id'"; }
		$digs_results_table->delete_rows($extracted_where);
		$deleted++;
	}

	# Iterate through the extracted sequences
	foreach my $locus_ref (@$extracted_ref) {
		
		# Insert the data to the digs_results table
		my $digs_result_id = $digs_results_table->insert_row($locus_ref);
		unless ($digs_result_id) { die; }
		if ($verbose) { print "\n\t\t\t # Created new result row '$digs_result_id'"; }
	}	
	
	# Flush the active set table
	$active_set_table->flush();

	# Return the number
	return $deleted;
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
	
	unless ($db_name) { 
		die "\n\t Undefined db_name in load_screening_db fxn\n\n\n";
	}

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
	$self->load_searches_table($dbh);	
	$self->load_active_set_table($dbh);	
	$self->load_digs_results_table($dbh, 'digs_results');	

	
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
		organism         => 'varchar',
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
        extract_start    => 'int',
        extract_end      => 'int',
		
		assigned_name    => 'varchar',
		assigned_gene    => 'varchar',
		subject_start    => 'int',
		subject_end      => 'int',
		query_start      => 'int',
		query_end        => 'int',
		align_len        => 'int',
		bitscore         => 'float',
		identity         => 'varchar',
		evalue_num       => 'float',
		evalue_exp       => 'int',
		gap_openings     => 'int',
        

		locus_structure  => 'text',
        sequence_length  => 'int',
        sequence         => 'text',
        
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
	print "\n\t ### Creating '$db_name' screening database\n\n";
    my $drh = DBI->install_driver("mysql");
    my $rc = $drh->func("createdb", $db_name, $server, $username, $password, 'admin');
    my $dbh = DBI->connect("dbi:mysql:$db_name:$server", $username, $password);
	unless ($dbh) {
		die "\n\t # Couldn't create screening database '$db_name'\n\n";
	}	
	$self->create_searches_table($dbh);
	$self->create_active_set_table($dbh);
	$self->create_digs_results_table($dbh);
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

	  `scaffold`         varchar(200) NOT NULL default '0',
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
	  `scaffold`         varchar(200) NOT NULL default '0',
	  `extract_start`    int(11) NOT NULL default '0',
	  `extract_end`      int(11) NOT NULL default '0',
	  `sequence_length`  int(11) NOT NULL default '0',
	  `sequence`         text NOT NULL,
	  `assigned_name`    varchar(100) NOT NULL default '0',
	  `assigned_gene`    varchar(100) NOT NULL default '0',
	  `orientation`      varchar(100) NOT NULL default '0',
	  `bitscore`         int(11) NOT NULL default '0',
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
# Subroutine:  create_loci_table
# Description: create MySQL 'loci' table
#***************************************************************************
sub create_loci_table {

    my ($self, $dbh) = @_;

    # Consolidated loci table 
    my $loci = "CREATE TABLE `loci` (
    
        `record_id`         int(11) NOT NULL auto_increment,
        `organism`          varchar(100) NOT NULL default '0',
        `target_datatype`   varchar(100) NOT NULL default '0',
        `target_version`    varchar(100) NOT NULL default '0',
        `target_name`       varchar(100) NOT NULL default '0',

        `scaffold`          varchar(200) default 'NULL',
        `orientation`       varchar(100) NOT NULL default '0',
        `extract_start`     int(11) NOT NULL default '0',
        `extract_end`       int(11) NOT NULL default '0',

        `assigned_name`     varchar(100) NOT NULL default '0',
        `assigned_gene`     varchar(100) NOT NULL default '0',
	    
        `subject_start`     int(11) NOT NULL default '0',
	`subject_end`       int(11) NOT NULL default '0',
	`query_start`       int(11) NOT NULL default '0',
	`query_end`         int(11) NOT NULL default '0',
	`align_len`         int(11) NOT NULL default '0',

	`bitscore`          float   NOT NULL default '0',
	`identity`          float   NOT NULL default '0',
	`evalue_num`        float   NOT NULL default '0',
	`evalue_exp`        int(11) NOT NULL default '0',
	`gap_openings`    int(11) NOT NULL default '0',
	`mismatches`      int(11) NOT NULL default '0',

        `locus_structure`   text NOT NULL,
        `sequence_length`   int(11) NOT NULL default '0',
        `sequence`          longtext NOT NULL,
      
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

	# Get tables
	my $active_set_table   = $self->{active_set_table};
	my $digs_results_table = $self->{digs_results_table};
	my $searches_table     = $self->{searches_table};
		
	# Flush result tables
	$active_set_table->flush();
	$active_set_table->reset_primary_keys();
	$digs_results_table->flush();
	$digs_results_table->reset_primary_keys();
	$searches_table->flush();
	$searches_table->reset_primary_keys();

}

############################################################################
# EXTEND SCREENING DB FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  do_track_table_dialogue
# Description: create an annotation track table for the screening DB
#***************************************************************************
sub do_track_table_dialogue {

	my ($self, $data_ref, $fields_array_ref, $fields_hash_ref) = @_;
	
	my $dbh = $self->{dbh};

	# Show the options
	my @choices = qw [ 1 2 3 ];
	print "\n\t\t 1. Create new track table";
	print "\n\t\t 2. Append data to existing track table";
	print "\n\t\t 3. Flush existing track table and import fresh data\n";
	my $question = "\n\t Choose an option";
	my $option = $console->ask_simple_choice_question($question, \@choices);

	# Declare the variables & data structures we need
	my $table_name;
	my $table_to_use;
	my %extra_tables;
	my @extra_tables;
	if ($option eq '1') { # Get name of new table
		my $table_name_question = "\n\t What is the name of the new table?";
		$table_name = $console->ask_question($table_name_question);
		$table_to_use = $table_name;
	}
	# or choose one of the ancillary tables already in the DB
	else {

		# Get the ancillary tables in this DB
		$self->get_ancillary_table_names(\@extra_tables);
		
		my $table_num = 0;
		foreach my $extra_table_name (@extra_tables) {
			$table_num++;
			$extra_tables{$table_num} = $extra_table_name;
			print "\n\t\t Table $table_num: '$extra_table_name'";
		}
		my @table_choices = keys %extra_tables;

		my $question5 = "\n\n\t Apply to which of the above tables?";
		my $answer5   = $console->ask_simple_choice_question($question5, \@table_choices);
		$table_to_use = $extra_tables{$answer5};
		unless ($table_to_use) { die; }
	}

	if ($option eq 1) { # Create table 
		$self->create_nomenclature_tracks_table($dbh, $table_to_use);
	}

	# Load the table
	$self->load_tracks_table($dbh, $table_to_use);
	my $track_table = $self->{nomenclature_tracks_table};
	unless ($track_table) {
		$devtools->print_hash($self);
		die;
	}
	
	if ($option eq 3)   {  # Flush the table if requested
		$track_table->flush();
		$track_table->reset_primary_keys();
	}

	# Insert the data	
	print "\n\n\t #### IMPORTING to table '$table_name'";
	my $row_count = $self->import_data_to_track_table($data_ref, $fields_array_ref, $fields_hash_ref);
	
	return $table_to_use;
}

#***************************************************************************
# Subroutine:  do_ancillary_table_dialogue
# Description: create an ancillary table for the screening DB
#***************************************************************************
sub do_ancillary_table_dialogue {

	my ($self, $fields_array_ref, $fields_hash_ref) = @_;
	
	# Declare the variables & data structures we need
	my %extra_tables;
	my @extra_tables;

	# Show the options
	my @choices = qw [ 1 2 3 ];
	print "\n\t\t 1. Create new ancillary table";
	print "\n\t\t 2. Append data to existing ancillary table";
	print "\n\t\t 3. Flush existing ancillary table and import fresh data\n";
	my $question = "\n\t Choose an option";
	my $option = $console->ask_simple_choice_question($question, \@choices);

	my $table_name;
	my $table_to_use;
	if ($option eq '1') { # Get name of new table
		my $table_name_question = "\n\t What is the name of the new table?";
		$table_name = $console->ask_question($table_name_question);
	}
	# or choose one of the ancillary tables already in the DB
	else {

		# Get the ancillary tables in this DB
		$self->get_ancillary_table_names(\@extra_tables);
		
		my $table_num = 0;
		foreach my $table_name (@extra_tables) {
			$table_num++;
			$extra_tables{$table_num} = $table_name;
			print "\n\t\t Table $table_num: '$table_name'";
		}
		my @table_choices = sort by_number keys %extra_tables;

		my $question5 = "\n\n\t Apply to which of the above tables?";
		my $answer5   = $console->ask_simple_choice_question($question5, \@table_choices);
		$table_to_use = $extra_tables{$answer5};
		unless ($table_to_use) { die; }
	}

	if ($option eq 1) { # Create table 
		$table_to_use = $self->create_ancillary_table($table_name, $fields_array_ref, $fields_hash_ref);	
	}

	my $dbh = $self->{dbh};
	my $anc_table = MySQLtable->new($table_to_use, $dbh, $fields_hash_ref);
	#my $table = $table_to_use . '_table';
	$self->{$table_to_use} = $anc_table;

	if ($option eq 3)   {  # Flush the table if requested
		$anc_table->flush();
		$anc_table->reset_primary_keys();
	}
	
	return $table_to_use;
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
			if ($item eq 'searches_performed')    { next; }
			elsif ($item eq 'active_set')         { next; }
			elsif ($item eq 'digs_results')       { next; }
			else {
				push (@$anc_tables_ref, $item)
			}
		}
	}
}

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
	
	my $copy_sql_1 = "CREATE TABLE $copy_name LIKE digs_results;";
	my $copy_sql_2 = "INSERT $copy_name SELECT * FROM digs_results";
    my $dbh = $self->{dbh};
	my $sth = $dbh->prepare($copy_sql_1);
   	unless ($sth->execute()) { print $copy_sql_1; exit;}	
   	$sth = $dbh->prepare($copy_sql_2);
   	unless ($sth->execute()) { print $copy_sql_2; exit;}
   	
   	return $copy_name;

}

#***************************************************************************
# Subroutine:  create_empty_digs_results_table 
# Description: 
#***************************************************************************
sub create_empty_digs_results_table {

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
	
	my $copy_sql_1 = "CREATE TABLE $copy_name LIKE digs_results;";
    my $dbh = $self->{dbh};
	my $sth = $dbh->prepare($copy_sql_1);
   	unless ($sth->execute()) { print $copy_sql_1; exit;}	
   	
   	return $copy_name;

}


############################################################################
# TESTING FOR EXISTENCE OF TABLES & DATABASES
############################################################################

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
# SCREENING DATABASE - Database table retrieval and insert functions
############################################################################

#**************************************************************************
# Subroutine:  get_sorted_digs_results
# Description: get digs_results table rows sorted by scaffold & coordinates
#***************************************************************************
sub get_sorted_digs_results {

	my ($self, $data_ref, $where) = @_;;

	# Set statement to sort loci
	my $sort  = " ORDER BY organism, target_name, scaffold, extract_start ";
	if ($where) { $where .= $sort; }
	else        { $where  = $sort; }

	my $digs_results_table  = $self->{digs_results_table};
		
	# Set the fields to get values for
	my @fields = qw [ record_id organism 
	                  target_name target_version target_datatype
	                  assigned_name assigned_gene probe_type
	                  scaffold orientation
	                  bitscore gap_openings
	                  query_start query_end
	                  subject_start subject_end 
	                  mismatches align_len
                      evalue_num evalue_exp identity 
                      extract_start extract_end sequence_length ];
                      
	$digs_results_table->select_rows(\@fields, $data_ref, $where);
	
	# Set record_ID as 'digs_result_id' in all results
	foreach my $row_ref (@$data_ref) {
		$row_ref->{digs_result_id} = $row_ref->{record_id};
	}
}

#***************************************************************************
# Subroutine:  get sorted active set 
# Description: get active set rows, sorted by scaffold, in order of location
#***************************************************************************
sub get_sorted_active_set {
	
	my ($self, $data_ref, $where) = @_;;

	# Set statement to sort loci
	if ($where) { $where .= " ORDER BY scaffold, subject_start "; }
	else        { $where  = " ORDER BY scaffold, subject_start "; }

	# Get database tables
	my $active_set_table    = $self->{active_set_table};
	
	# Get sorted, combined extracted loci and new blast results	
	my @blast_fields = qw [ record_id digs_result_id  
	                        organism target_datatype target_version target_name
	                        probe_name probe_gene probe_type
	                        bitscore gap_openings
	                        query_start query_end 
	                        align_len mismatches 
                            evalue_num evalue_exp identity 
	                        scaffold orientation
	                        subject_start subject_end ];
	$active_set_table->select_rows(\@blast_fields, $data_ref, $where);
}

#***************************************************************************
# Subroutine:  insert_row_in_active_set_table
# Description: insert a BLAST result as a row into the active set table
#***************************************************************************
sub insert_row_in_active_set_table {

	my ($self, $query_ref, $hit_ref) = @_;

	# Get screening database table objects
	my $active_set_table = $self->{active_set_table};

	my $probe_id        = $query_ref->{probe_id};
	my $probe_name      = $query_ref->{probe_name};
	my $probe_gene      = $query_ref->{probe_gene};
	my $probe_type      = $query_ref->{probe_type};
	my $probe_path      = $query_ref->{probe_path};
	my $organism        = $query_ref->{organism};
	my $version         = $query_ref->{target_version};
	my $datatype        = $query_ref->{target_datatype};
	my $target_name     = $query_ref->{target_name};
	my $target_path     = $query_ref->{target_path};

	$hit_ref->{digs_result_id}  = 0;
	$hit_ref->{organism}        = $organism;
	$hit_ref->{target_version}  = $version;
	$hit_ref->{target_datatype} = $datatype;
	$hit_ref->{target_name}     = $target_name;
	$hit_ref->{probe_id}        = $probe_id;
	$hit_ref->{probe_name}      = $probe_name;
	$hit_ref->{probe_gene}      = $probe_gene;
	$hit_ref->{probe_type}      = $probe_type;
	$hit_ref->{subject_start}   = $hit_ref->{aln_start};  # Rename to match DB
	$hit_ref->{subject_end}     = $hit_ref->{aln_stop};   # Rename to match DB
	$hit_ref->{query_end}       = $hit_ref->{query_stop}; # Rename to match DB
	$active_set_table->insert_row($hit_ref);

}

#***************************************************************************
# Subroutine:  add_digs_results_to_active_set 
# Description: get all extracted loci for this target file (& probe)
#***************************************************************************
sub add_digs_results_to_active_set {
	
	my ($self, $data_ref) = @_;;

	# Get database tables
	my $active_set_table  = $self->{active_set_table};

	# Enter all relevant extracted loci into 'active_set' table 
	my $num_loci = scalar @$data_ref;
	foreach my $locus_ref (@$data_ref) {

		#print "\n\t\t # inserting extract ID $digs_result_id";
		#$devtools->print_hash($locus_ref);
		my $digs_result_id = $locus_ref->{record_id};
		
		# Translations
		$locus_ref->{digs_result_id}       = $digs_result_id;
		$locus_ref->{probe_name}       = $locus_ref->{assigned_name};
		$locus_ref->{probe_gene}       = $locus_ref->{assigned_gene};
		$locus_ref->{subject_start}    = $locus_ref->{extract_start};
		$locus_ref->{subject_end}      = $locus_ref->{extract_end};
		$locus_ref->{organism}  = $locus_ref->{organism};
		$active_set_table->insert_row($locus_ref);
	}
	
	return $num_loci;
}

############################################################################
# Importing data
############################################################################

#***************************************************************************
# Subroutine:  import_data_to_ancillary_table
# Description: load data to ancillary table
#***************************************************************************
sub import_data_to_ancillary_table {

	my ($self, $table_name, $data_ref, $fields_array_ref, $fields_hash_ref, $verbose) = @_;

	my $anc_table = $self->{$table_name};
	unless ($anc_table) { 
		#$devtools->print_hash($self);
		print "\n\t couldn't get table '$table_name'\n\n\n";
		die; 
	}
	#$devtools->print_hash($fields_hash_ref);
	#$devtools->print_array($data_ref); die;

	my $row_count = 0;
	foreach my $line (@$data_ref) { # Add data to the table
		$row_count++;
		chomp $line;
		my %insert;
		my @elements = split ("\t", $line);
		my $column_num = 0;
		foreach my $field (@$fields_array_ref) {
			my $value = $elements[$column_num];
			$column_num++;
			my $type  = $fields_hash_ref->{$column_num};
			if ($verbose) {
				print "\n\t Row count $row_count: importing value '$value' to field '$field'";
			}
			unless ($value) { 
				$value = 'NULL';
			}
			$insert{$field} = $value;
		}
		$anc_table->insert_row(\%insert);
	}
	
	return $row_count;
}

#***************************************************************************
# Subroutine:  import_data_to_track_table
# Description: load data to annotation track table
#***************************************************************************
sub import_data_to_track_table {

	my ($self, $data_ref, $fields_array_ref, $fields_hash_ref) = @_;

	my $tracks_table = $self->{nomenclature_tracks_table};
	unless ($tracks_table) { 
		$devtools->print_hash($self);
		print "\n\t couldn't get track table'\n\n\n";
		die; 
	}
	#$devtools->print_hash($fields_hash_ref);
	#$devtools->print_array($data_ref); die;

	my $row_count = 0;
	foreach my $line (@$data_ref) { # Add data to the table
		$row_count++;
		chomp $line;
		my %insert;
		my @elements = split ("\t", $line);
		my $column_num = 0;
		foreach my $field (@$fields_array_ref) {
			my $value = $elements[$column_num];
			$column_num++;
			my $type  = $fields_hash_ref->{$column_num};
			unless ($value) { 
				$value = 'NULL';
			}
			print "\n\t Row count $row_count: importing value '$value' to field '$field'";
			if ($field eq 'track_name'
			or  $field eq 'organism_code'
			or  $field eq 'locus_class'
			or  $field eq 'taxon'
			or  $field eq 'organism_code'
			or  $field eq 'gene'
			or  $field eq 'scaffold'
			or  $field eq 'start_position'
			or  $field eq 'end_position'
			or  $field eq 'orientation'
			or  $field eq 'namespace_id') {
				$insert{$field} = $value;
			}
		}
		$tracks_table->insert_row(\%insert);
	}
	
	return $row_count;
}

#***************************************************************************
# Subroutine:  import_data_to_searches_performed
# Description: import data to the 'searches_performed' table
#***************************************************************************
sub import_data_to_searches_performed {

	my ($self,  $data_path) = @_;

	# Get the table object
	my $searches_table = $self->{searches_table};

	# Get the data	
	my @data;
	$fileio->read_file(\@data, $data_path);

	# Iterate through the data and insert to the table
	foreach my $line (@data) {
		my %data;
		chomp $line;
		my @line = split("\t", $line);
		$data{probe_id}        = shift @line;
		$data{probe_name}      = shift @line;
		$data{probe_gene}      = shift @line;
		$data{target_id}       = shift @line;
		$data{organism}        = shift @line;
		$data{target_datatype} = shift @line;
		$data{target_version}  = shift @line;
		$data{target_name}     = shift @line;
		$searches_table->insert_row(\%data);	
	}
}

#***************************************************************************
# Subroutine:  import_data_to_digs_results
# Description: import data to the 'digs_results' table
#***************************************************************************
sub import_data_to_digs_results {

	my ($self, $data_path, $alt_table_name) = @_;

	# Get the table object
	my $digs_results_table;
	if ($alt_table_name) { 
		$digs_results_table = $self->{$alt_table_name}; 
	}
	else {
		$digs_results_table = $self->{digs_results_table};
	}
	
    # Read in the data from a tab delimited file
	my @data;
    my %fields;
    my @fields;
    $console->do_read_tabdelim_dialogue($data_path, \@data, \@fields, \%fields);
    #$devtools->print_array(\@data); die;
	#my $read = $fileio->read_file($data_path, \@data);
    #unless ($read) { die "\n\t\t Couldn't open file '$data_path'\n\n"; }
	
	# Iterate through the data and insert to the table
	my $rows = '0';
	foreach my $line (@data) {

		chomp $line;
		$line =~ s/"//g;
		my @line = split("\t", $line);
		my %data;
		$data{record_id}       = shift @line;
		$data{organism}        = shift @line;
		$data{target_datatype} = shift @line;
		$data{target_version}  = shift @line;
		$data{target_name}     = shift @line;
		$data{probe_type}      = shift @line;
		$data{scaffold}        = shift @line;
		$data{extract_start}   = shift @line;
		$data{extract_end}     = shift @line;
		$data{sequence_length} = shift @line;
		$data{sequence}        = shift @line;
		$data{assigned_name}   = shift @line;
		$data{assigned_gene}   = shift @line;
		$data{orientation}     = shift @line;
		$data{bitscore}        = shift @line;
		$data{identity}        = shift @line;
		$data{evalue_num}      = shift @line;
		$data{evalue_exp}      = shift @line;
		$data{subject_start}   = shift @line;
		$data{subject_end}     = shift @line;
		$data{query_start}     = shift @line;
		$data{query_end}       = shift @line;
		$data{align_len}       = shift @line;
		$data{gap_openings}    = shift @line;
		$data{mismatches}      = shift @line;
		$devtools->print_hash(\%data);
		$digs_results_table->insert_row(\%data);	
		$rows++;
	}	
	return $rows;
}

############################################################################
# UPDATE DB SCHEMA
############################################################################

#***************************************************************************
# Subroutine:  translate_schema
# Description: convert DIGS DB with previous schema into one with the new one
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

	# Translate 'Extracted' to 'digs_results'
	$self->load_extracted_table($dbh);
	my $extracted_table    = $self->{extracted_table};
	my $digs_results_table = $self->{digs_results_table};
	$digs_results_table->flush();
	my @extracted_rows;
	my @extracted_fields = qw [ organism version data_type target_name probe_type 
	                            assigned_name assigned_gene scaffold 
	                            extract_start extract_end orientation
                                subject_start subject_end query_start query_end 
	                            identity bit_score e_value_num e_value_exp align_len 
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
		$row_ref->{target_id}       = $row_ref->{genome_id};
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

############################################################################
# Basic
############################################################################

#***************************************************************************
# Subroutine:  by number
# Description: sort an array of integers by ascending numerical order 
#***************************************************************************
sub by_number { $a <=> $b }	

############################################################################
# EOF
############################################################################
