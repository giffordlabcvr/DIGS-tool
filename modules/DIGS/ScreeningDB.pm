#!/usr/bin/perl -w
############################################################################
# Module:      ScreeningDB.pm
# Description: The DB module for Pipeline databases
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
use Base::SeqIO;
use Base::FileIO;
use Base::DevTools;
use Base::Console;

# Components
use DIGS::TargetDB;

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
		server                => $parameters->{mysql_server},
		username              => $parameters->{mysql_username},
		password              => $parameters->{mysql_password},
		
		# Screening database core tables
		status_table          => 0,
		blast_results_table   => 0,
		extracted_table       => 0,
		loci_table            => 0,
		loci_link_table       => 0,
	
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# LOAD DATABASE & DATABASE TABLES
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
	print "\n\n\t ### Loading '$db_name' screening database";
	$self->load_status_table($dbh);	
	$self->load_blast_results_table($dbh);	
	$self->load_extracted_table($dbh);	
	$self->load_loci_table($dbh);
	$self->load_loci_link_table($dbh);

	# Check integrity of database
	my $extracted_count = $self->count_extracted_rows();
	my $blast_count     = $self->count_blast_rows();
	unless ($extracted_count eq $blast_count) { # Tables are out of sync
		print "\n\n\t Tables are out of sync";
		print "\n\t Extracted table:     $extracted_count rows";
		print "\n\t BLAST_results table: $blast_count rows";
		print "\n\t Rolling back screen...";
		$self->validate_db();
		sleep 3;
	}
	
	# If no results for most recent query, execute again just incase
	#$self->rollback_last_search();
}

############################################################################
# LOAD SCREENING DATABASE TABLES
############################################################################

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
		genome_id      => 'varchar',
		organism       => 'varchar',
		data_type      => 'varchar',
		version        => 'varchar',
		target_name    => 'varchar',
	);
	my $status_table = MySQLtable->new('Status', $dbh, \%status_fields);
	$self->{status_table} = $status_table;
}

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
		#realex_start     => 'int',
        #realex_end       => 'int',	
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
# Subroutine:  load_loci_table
# Description: load screening database table 'Loci'
#***************************************************************************
sub load_loci_table {

	my ($self, $dbh) = @_;

	# Definition of the loci table
	my %loci_fields = (
    	organism         => 'varchar',
        version          => 'varchar',
       	data_type        => 'varchar',
        target_name      => 'varchar',
        scaffold         => 'varchar',
        orientation      => 'varchar',
        assigned_name    => 'varchar',
        extract_start    => 'int',
        extract_end      => 'int',
        genome_structure => 'text',

	);
	my $loci_table = MySQLtable->new('Loci', $dbh, \%loci_fields);
	$self->{loci_table} = $loci_table;
}

#***************************************************************************
# Subroutine:  load_loci_link_table
# Description: load screening database table 'Loci'
#***************************************************************************
sub load_loci_link_table {

	my ($self, $dbh) = @_;

	# Definition of the loci table
	my %loci_fields = (
		locus_id       => 'int',
		extracted_id   => 'int',
	);
	my $loci_table = MySQLtable->new('Loci_link', $dbh, \%loci_fields);
	$self->{loci_link_table} = $loci_table;
}

############################################################################
# CREATE DATABASE TABLES
############################################################################

#***************************************************************************
# Subroutine:  create_loci_table
# Description: create MySQL 'Loci' table
#***************************************************************************
sub create_loci_table {

	my ($self, $dbh) = @_;

	# Loci table 
    my $loci = "CREATE TABLE `Loci` (
        `Record_ID`         int(11) NOT NULL auto_increment,
        `Organism`      varchar(100) NOT NULL default '0',
        `Data_type`     varchar(100) NOT NULL default '0',
        `Version`       varchar(100) NOT NULL default '0',
        `Target_name`   varchar(100) NOT NULL default '0',

        `Scaffold`      varchar(100) default 'NULL',
        `Orientation`   varchar(100) NOT NULL default '0',

        `Assigned_name`    varchar(100) NOT NULL default '0',
        `Extract_start`     int(11) NOT NULL default '0',
        `Extract_end`       int(11) NOT NULL default '0',
        `Genome_structure`  text NOT NULL,
	  
		`Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($loci);
	unless ($sth->execute()) { print "\n\t$loci\n\n\n"; exit;}
}

#***************************************************************************
# Subroutine:  create_loci_link_table
# Description: create MySQL 'Loci' table
#***************************************************************************
sub create_loci_link_table {

	my ($self, $dbh) = @_;

	# Loci table 
    my $loci = "CREATE TABLE `Loci_link` (
      `Record_ID`       int(11) NOT NULL auto_increment,
      `Locus_ID`        int(11) NOT NULL default '0',
      `Extracted_ID`    int(11) NOT NULL default '0',
      `Timestamp`       timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
      PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($loci);
	unless ($sth->execute()) { print "\n\t$loci\n\n\n"; exit;}
}

############################################################################
# CREATE SCREENING DATABASE TABLES
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
	print "\n\t ### Creating '$db_name' screening database\n";
    my $drh = DBI->install_driver("mysql");
    my $rc = $drh->func("createdb", $db_name, $server, $username, $password, 'admin');
    my $dbh = DBI->connect("dbi:mysql:$db_name:$server", $username, $password);
	unless ($dbh) {
		die "\n\t # Couldn't create screening database '$db_name'\n\n";
	}	

	$self->create_status_table($dbh);
	$self->create_blast_results_table($dbh);
	$self->create_extracted_table($dbh);
	$self->create_loci_table($dbh);
	$self->create_loci_link_table($dbh);
}

#***************************************************************************
# Subroutine:  create_status_table
# Description: create MySQL 'Status' table
#***************************************************************************
sub create_status_table {

	my ($self, $dbh) = @_;

	# Status table (which BLAST queries have been executed)
	my $status = "CREATE TABLE `Status` (
	  `Record_ID`   int(11) NOT NULL auto_increment,
	  `Probe_ID`    varchar(100) NOT NULL default '',
	  `Probe_name`  varchar(100) NOT NULL default '',
	  `Probe_gene`  varchar(100) NOT NULL default '',
	  `Genome_ID`   varchar(100) NOT NULL default '',
	  `Organism`    varchar(100) NOT NULL default '',
	  `Data_type`  varchar(100) NOT NULL default '',
	  `Version`     varchar(100) NOT NULL default '',
	  `Target_name` varchar(100) NOT NULL default '',
	  `Timestamp`   timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1";
	my $sth = $dbh->prepare($status);
	unless ($sth->execute()) { print "\n\t$status\n\n\n"; exit;}
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
# Subroutine:  create_loci_table
# Description: create MySQL 'Loci' table
#***************************************************************************
sub create_loci_table {

	my ($self, $dbh) = @_;

	# Loci table 
    my $loci = "CREATE TABLE `Loci` (
	  `Record_ID`         int(11) NOT NULL auto_increment,
	  `Assigned_to`       varchar(100) NOT NULL default '0',
	  `Assigned_notes`    text NOT NULL,
	  `Extract_start`     int(11) NOT NULL default '0',
	  `Extract_end`       int(11) NOT NULL default '0',
	  `Organism`          varchar(100) NOT NULL default '0',
	  `Orientation`       varchar(100) NOT NULL default '0',
	  `Chunk_name`        varchar(100) NOT NULL default '0',
	  `Scaffold`          varchar(100) default 'NULL',
	  `Genome_structure`  text NOT NULL,
	  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($loci);
	unless ($sth->execute()) { print "\n\t$loci\n\n\n"; exit;}
}

#***************************************************************************
# Subroutine:  create_loci_link_table
# Description: create MySQL 'Loci' table
#***************************************************************************
sub create_loci_link_table {

	my ($self, $dbh) = @_;

	# Loci table 
    my $loci = "CREATE TABLE `Loci_link` (
      `Record_ID`       int(11) NOT NULL auto_increment,
      `Locus_ID`        int(11) NOT NULL default '0',
      `Extracted_ID`    int(11) NOT NULL default '0',
      `Timestamp`       timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
      PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($loci);
	unless ($sth->execute()) { print "\n\t$loci\n\n\n"; exit;}
}

############################################################################
# UPDATING DATABASES AND TABLES
############################################################################

#***************************************************************************
# Subroutine:  index_previously_executed_queries 
# Description: index BLAST searches that have previously been executed
#***************************************************************************
sub index_previously_executed_queries {
	
	my ($self, $done_ref) = @_;
	
	my $status_table = $self->{status_table};
	unless ($status_table) { die "\n\t Status table not loaded\n\n"; }
	my @data;
	my @fields = qw [ record_id organism data_type version target_name
                      probe_name probe_gene ];
	my $where = " ORDER BY Record_ID ";
	$status_table->select_rows(\@fields, \@data, $where);
	
	# Index the executed searches
	foreach my $data_ref (@data) {
		my $organism    = $data_ref->{organism};
		my $data_type   = $data_ref->{data_type};
		my $version     = $data_ref->{version};
		my $target_name = $data_ref->{target_name};
		my $probe_name  = $data_ref->{probe_name};
		my $probe_gene  = $data_ref->{probe_gene};
		#$devtools->print_hash($data_ref); die; # DEBUG
		unless ( $organism and $data_type and $version and $target_name 
             and $probe_name and $probe_gene) { 
			die; 
		};
		my @genome = ( $organism , $data_type, $version );
		my $genome_id = join ('|', @genome);
		my $probe_id  = $probe_name . '_' .  $probe_gene;
		my @key = ( $genome_id, $target_name, $probe_id );
		my $key = join ('|', @key);
		#print "\n\t $key\n\n "; die;
		$done_ref->{$key} = 1;		
	}
	# DEBUG
	#$devtools->print_hash($done_ref); die;
}

############################################################################
# RETRIEVING DATA
############################################################################

#***************************************************************************
# Subroutine:  summarise_db
# Description: summarise a screening database 
#***************************************************************************
sub summarise_db {

	my ($self) = @_;
	
	my $db_name = $self->{db_name};
	
	# Summarise status table
	my $executed = $self->summarise_status_table();	
	$self->{status_table_count} = $executed;
	if ($executed) {

		print "\n\n\t ### Summarizing '$db_name' screening database";
		
		# Summarise BLAST_results  table
		$self->summarise_BLAST_results_table();
	
		# Summarise Extracted  table
		$self->summarise_extracted_table();
		print "\n\n";
	}
}

#***************************************************************************
# Subroutine:  summarise_status_table
# Description: summarise data in the Status table
#***************************************************************************
sub summarise_status_table {

	my ($self, $done_ref) = @_;
	
	my $status_table = $self->{status_table};
	my @data;
	my @fields = qw [ organism probe_name probe_gene ];
	$status_table->select_rows(\@fields, \@data);
	my $executed;
	my %organism;
	my %probe_id;
	foreach my $data_ref (@data) {

		my $organism   = $data_ref->{organism};
		my $probe_name = $data_ref->{probe_name};
		my $probe_gene = $data_ref->{probe_gene};
		my $probe_id   = $probe_name . '_' . $probe_gene;
	
		# Counts
		$executed++;
		# By organism
		if ($organism{$organism}) {
			$organism{$organism}++;
		}
		else {
			$organism{$organism} = 1;
		}
		# By probe ID
		if ($probe_id{$probe_id}) {
			$probe_id{$probe_id}++;
		}
		else {
			$probe_id{$probe_id} = 1;
		}
	}

	# Show data
	unless ($executed) { return 0; }
	print "\n\n\t # A total of $executed queries have previously been executed";
	print "\n\t #  ---";
	my @organisms = sort keys %organism; 
	foreach my $organism (@organisms) {
		my $by_organism = $organism{$organism};
		print "\n\t #  $by_organism searches of $organism databases";
	}
	print "\n\t #  ---";
	my @probe_ids = sort keys %probe_id; 
	foreach my $probe_id (@probe_ids) {
		my $by_probe_id = $probe_id{$probe_id};
		print "\n\t #  $by_probe_id searches using $probe_id";
	}
	sleep 1;
	return $executed;
}

#***************************************************************************
# Subroutine:  summarise BLAST_results table
# Description: summarise data in the BLAST_results table
#***************************************************************************
sub summarise_BLAST_results_table {

	my ($self, $done_ref) = @_;
	
	my $blast_table = $self->{blast_results_table};
	my @data;
	my @fields = qw [ organism target_name  ];
	push (@fields,  "count(*) AS 'number'");
	my $where = "GROUP BY  Organism, Target_name
                 ORDER BY  Organism, count(*) DESC";
	$blast_table->select_rows(\@fields, \@data, $where);
	
	my $blast = 0;
	foreach my $data_ref (@data) {
		my $number = $data_ref->{number};
		$blast = $blast + $number;
	}
	print "\n\n\t # The BLAST results table contains a total of '$blast' rows";
	print "\n\t #  ---";
	foreach my $data_ref (@data) {
		# get the data	
		my $organism    = $data_ref->{organism};
		my $chunk_name  = $data_ref->{target_name};
		my $number      = $data_ref->{number};
		print "\n\t #  $number hits in:";
		print "    $organism genome, target file $chunk_name\t";
	}
	sleep 1;
}

#***************************************************************************
# Subroutine:  summarise Extracted table
# Description: summarise data in the Extracted table
#***************************************************************************
sub summarise_extracted_table {

	my ($self, $done_ref) = @_;
	
	my $extracted_table = $self->{extracted_table};
	my @data;
	my @fields = qw [ organism assigned_name assigned_gene ];
	push (@fields,  "count(*) AS 'number'");
	my $where = "GROUP BY  Organism, Assigned_name, Assigned_gene 
                 ORDER BY  Organism, count(*) DESC";
	$extracted_table->select_rows(\@fields, \@data, $where);
	
	my $extracted = 0;
	foreach my $data_ref (@data) {
		my $number = $data_ref->{number};
		$extracted = $extracted + $number;
	}
	print "\n\n\t # The extracted table contains a total of '$extracted' rows";
	print "\n\t #  ---";
	foreach my $data_ref (@data) {
		# get the data	
		my $organism         = $data_ref->{organism};
		my $assigned_name    = $data_ref->{assigned_name};
		my $assigned_gene    = $data_ref->{assigned_gene};
		my $number           = $data_ref->{number};
		print "\n\t #  $number matches to:\t";
		print "$assigned_name, $assigned_gene";
		print "\t in $organism";
	}
	sleep 1;
}

#***************************************************************************
# Subroutine:  retrieve_sequences
# Description: get sequences and create FASTA
#***************************************************************************
sub retrieve_sequences {

	my ($self, $data_ref, $select_ref, $where) = @_;
	
	my $extracted_table = $self->{extracted_table};
	unless ($extracted_table) { die; }
	my @data;
	$extracted_table->select_rows($select_ref, \@data, $where);
	foreach my $row (@data) {
		my $sequence = $row->{sequence};
		my @header;
		#$devtools->print_hash($row);
		foreach my $field (@$select_ref) {
			if ($field eq 'sequence') { next; }
			my $value = $row->{$field};
			push (@header, $value);
		}
		my $header   = join("_", @header);
		my $fasta    = ">$header\n$sequence\n\n";
		push (@$data_ref, $fasta);
	}
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
	
	# Ask to make sure
	my $question = "\n\n\t Are you sure you want to DROP the $db_name database?";
	my $answer1 = $console->ask_yes_no_question($question);
	if ($answer1 eq 'y') {
		#sleep 3;
		
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
		my $status_table         = $self->{status_table};
		
		# Flush result tables
		$blast_results_table->flush();
		$blast_results_table->reset_primary_keys();
		$extracted_table->flush();
		$extracted_table->reset_primary_keys();
		$status_table->flush();
		$status_table->reset_primary_keys();
	}
}

#***************************************************************************
# Subroutine:  count_blast_rows
# Description: count rows in the BLAST_results table 
#***************************************************************************
sub count_blast_rows {

	my ($self) = @_;

	my $blast_results_table = $self->{blast_results_table};
	my @fields = qw [ COUNT(*) ];
	my @data;
	$blast_results_table->select_rows(\@fields, \@data);
	my $data_ref = shift @data;
	my $count = $data_ref->{'COUNT(*)'};
	#print "\n\t BLAST results table: $count rows";
	return $count;	
}

#***************************************************************************
# Subroutine:  count_extracted_rows
# Description: count rows in the extracted table 
#***************************************************************************
sub count_extracted_rows {

	my ($self) = @_;

	my $extracted_table = $self->{extracted_table};
	my @fields = qw [ COUNT(*) ];
	my @data;
	$extracted_table->select_rows(\@fields, \@data);
	my $data_ref = shift @data;
	my $count = $data_ref->{'COUNT(*)'};
	#print "\n\t Extracted table: $count rows";
	return $count;	
}

############################################################################
# STANDARD SCREENING DB TABLE FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  index BLAST results by record id
# Description: Index loci in BLAST_results table by the 'record_id' field
#***************************************************************************
sub index_BLAST_results_by_record_id {
	
	my ($self, $data_ref, $where) = @_;

	# Get relevant variables and objects
	my $blast_table = $self->{blast_results_table}; 
	my @fields = qw [ record_id 
	                  organism data_type version target_name 
                      probe_name probe_gene ];
	my @record_ids;
	$blast_table->select_rows(\@fields, \@record_ids, $where);	
	foreach my $hit_ref (@record_ids) {
		my $record_id = $hit_ref->{record_id};
		if ($data_ref->{$record_id}) { die; } # BLAST ID should be unique
		$data_ref->{$record_id} = $hit_ref;	
	}
}

#***************************************************************************
# Subroutine:  index extracted loci by BLAST id
# Description: Index loci in Extracted table by the 'blast_id' field
#***************************************************************************
sub index_extracted_loci_by_blast_id {
	
	my ($self, $previously_extracted_ref, $where) = @_;

	# Get relevant variables and objects
	my $extracted_table = $self->{extracted_table}; 
	my @fields = qw [ blast_id ];
	my @blast_ids;
	$extracted_table->select_rows(\@fields, \@blast_ids, $where);	
	foreach my $hit_ref (@blast_ids) {
		my $blast_id = $hit_ref->{blast_id};
		if ($previously_extracted_ref->{$blast_id}) { die; } # BLAST ID should be unique
		$previously_extracted_ref->{$blast_id} = 1;	
	}
}

#***************************************************************************
# Subroutine:  get_blast_hits_to_extract
# Description: Get BLAST hits to extract
#***************************************************************************
sub get_blast_hits_to_extract {
	
	my ($self, $query_ref, $hits_ref) = @_;

	# Get data structures and variables from self
	my $blast_results_table = $self->{blast_results_table};

	# Get parameters for this query
	my $probe_name  = $query_ref->{probe_name};
	my $probe_gene  = $query_ref->{probe_gene};
	my $target_name = $query_ref->{target_name}; # target file name

	# Get all BLAST results from table (ordered by sequential targets)
	my $where = " WHERE Target_name = '$target_name'
                  AND probe_name = '$probe_name' 
                  AND probe_gene = '$probe_gene'
	              ORDER BY scaffold, subject_start";

	my @fields = qw [ record_id 
	               organism data_type version 
                   probe_name probe_gene probe_type
				   orientation scaffold target_name
                   subject_start subject_end 
		           query_start query_end ];
	$blast_results_table->select_rows(\@fields, $hits_ref, $where); 
}

############################################################################
# VALIDATION FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  validate_db
# Description: validate screening database 
#***************************************************************************
sub validate_db {

	my ($self) = @_;
	
	# Get tables	
	my $extracted_table = $self->{extracted_table};
	my $blast_table     = $self->{blast_results_table};
	my $status_table    = $self->{status_table};

	# Index data in the BLAST results and Extracted tables	
	my %blast_results;
	$self->index_BLAST_results_by_record_id(\%blast_results);
	my %extracted;
	$self->index_extracted_loci_by_blast_id(\%extracted);

	# Check for BLAST_results table rows that lack counterparts in the Extracted table
	my @missing;
	my @ids = sort by_number keys %blast_results;
	foreach my $record_id (@ids) {	
		unless ($extracted{$record_id}) {
			my $data_ref = $blast_results{$record_id};
			unless ($data_ref) { die; }
			push (@missing, $data_ref); 
		}
	}
	#$devtools->print_array(\@missing); die;

	# Check for Extracted table rows that lack counterparts in the BLAST_results table
	my @blast_ids = sort by_number keys %extracted;
	foreach my $blast_id (@blast_ids) {	
		unless ($blast_results{$blast_id}) {
			my $where = " WHERE BLAST_ID = $blast_id ";
			$extracted_table->delete_rows($where);
		}
	}
	
	# Rollback any searches where no sequence was captured
	my %rollback_searches;
	foreach my $missing_ref (@missing) {
		# Store the search information
		my $record_id   = $missing_ref->{record_id};
		my $organism    = $missing_ref->{organism};
		my $data_type   = $missing_ref->{data_type};
		my $version     = $missing_ref->{version};
		my $target_name = $missing_ref->{target_name};
		my $probe_name  = $missing_ref->{probe_name};
		my $probe_gene  = $missing_ref->{probe_gene};
		my @target = ( $organism , $data_type, $version, $target_name );
		my $target_id = join ('|', @target);
		my @key = ( $target_id, $probe_name, $probe_gene );
		my $key = join ('|', @key);
		$rollback_searches{$key} = $missing_ref;
	}
	#$devtools->print_hash(\%rollback_searches); die;

	my @keys = keys %rollback_searches;
	foreach my $key (@keys) {
		print "\n\n\t Cleaning up incomplete search...";
		print "\n\t\t KEY: '$key'\n";
		my $missing_ref = $rollback_searches{$key};
		$self->rollback($missing_ref);
	}
	
	my $extracted_count = $self->count_extracted_rows();
	my $blast_count     = $self->count_blast_rows();
	unless ($extracted_count eq $blast_count) { die; }
	print "\n\t Tables re-synced";
	print "\n\t Extracted table:     $extracted_count rows";
	print "\n\t BLAST_results table: $blast_count rows";
	sleep 3;

}

#***************************************************************************
# Subroutine:  safe_rollback
# Description: 
#***************************************************************************
sub rollback_last_search {

	my ($self) = @_;
	
	# Get tables	
	my $status_table    = $self->{status_table};
	my @data;
	my @fields = qw [ record_id 
	                  organism data_type version target_name
                      probe_name probe_gene ];
	my $where  = " ORDER BY Record_ID ";
	$status_table->select_rows(\@fields, \@data, $where);
	my $last_search = pop @data;
	if ($last_search) {	
		$self->rollback($last_search);
	}
}

#***************************************************************************
# Subroutine:  rollback
# Description: 
#***************************************************************************
sub rollback {

	my ($self, $missing_ref) = @_;

	# Get tables	
	my $extracted_table = $self->{extracted_table};
	my $blast_table     = $self->{blast_results_table};
	my $status_table    = $self->{status_table};

	my $record_id   = $missing_ref->{record_id};
	my $organism    = $missing_ref->{organism};
	my $data_type   = $missing_ref->{data_type};
	my $version     = $missing_ref->{version};
	my $target_name = $missing_ref->{target_name};
	my $probe_name  = $missing_ref->{probe_name};
	my $probe_gene  = $missing_ref->{probe_gene};
	unless ($organism and $data_type and $version and $target_name) { die; }
	unless ($record_id and $probe_name and $probe_gene) { die; }

	my @genome = ( $organism , $data_type, $version );
	my $genome_id = join ('|', @genome);
	my @probe = ( $probe_name, $probe_gene );
	my $probe_id = join ('_', @probe); 
	  
	# Delete rows from the Status table for this search
	my $where1 = " WHERE probe_id    = '$probe_id'
					 AND genome_id   = '$genome_id'
					 AND target_name = '$target_name'";
	$status_table->delete_rows($where1);
	#print "\n\t '$where1' ";
	
	# Delete rows from the BLAST_results and Extracted tables for this search
	my @data;
	my @fields = qw [ record_id ];
	my $where2 =   " WHERE organism  = '$organism'
					 AND data_type   = '$data_type'
					 AND version     = '$version'
					 AND target_name = '$target_name'
					 AND probe_name  = '$probe_name'
					 AND probe_gene  = '$probe_gene'";
	#print "\n\t '$where2' ";
	$blast_table->select_rows(\@fields, \@data, $where2);
	foreach my $data_ref (@data) {
		my $record_id = $data_ref->{record_id};
		#print "\n\t Record ID = '$record_id'";
		my $id_where1 = " WHERE record_id = $record_id ";
		$blast_table->delete_rows($id_where1);
		my $id_where2 =   " WHERE blast_id  = $record_id ";
		$extracted_table->delete_rows($id_where2);
	}
}

#***************************************************************************
# Subroutine:  by_number
# Description: by number - for use with perl 'sort'  (cryptic but works) 
#***************************************************************************
sub by_number { $a <=> $b }	

############################################################################
# EOF
############################################################################
