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
		server                => $parameters->{mysql_server},
		username              => $parameters->{mysql_username},
		password              => $parameters->{mysql_password},
		
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
	unless ($dbh) { die "\n\t Failed to connect to database\n\n\n"; }
	$self->{dbh} = $dbh;

	# Main Screening DB tables
	print "\n\n\t ### Loading '$db_name' screening database";
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
		data_type        => 'varchar',
		version          => 'varchar',
		target_name       => 'varchar',
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
		organism         => 'varchar',
		version          => 'varchar',
		data_type        => 'varchar',
		target_name      => 'varchar',
		assigned_name      => 'varchar',
		assigned_gene => 'varchar',
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
	print "\n\n\t ### Creating '$db_name' screening database";
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
	  `Data_type`     varchar(100) NOT NULL default '0',
	  `Version`       varchar(100) NOT NULL default '0',
	  `Target_name`   varchar(100) NOT NULL default '0',

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
	  `Organism`         varchar(100) NOT NULL default '0',
	  `Data_type`        varchar(100) NOT NULL default '0',
	  `Version`          varchar(100) NOT NULL default '0',
	  `Target_name`      varchar(100) NOT NULL default '0',
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
# Description: 
#***************************************************************************
sub summarise_db {

	my ($self) = @_;
	
	my $db_name = $self->{db_name};
	print "\n\n\t ### Summarizing '$db_name' screening database";
	
	# Summarise status table
	my $executed = $self->summarise_status_table();	
	$self->{status_table_count} = $executed;
	if ($executed) {

		# Summarise BLAST_results  table
		$self->summarise_BLAST_results_table();
	
		# Summarise Extracted  table
		$self->summarise_extracted_table();
	}
}

#***************************************************************************
# Subroutine:  summarise_status_table
# Description: 
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
# Description: 
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
# Description: 
#***************************************************************************
sub summarise_extracted_table {

	my ($self, $done_ref) = @_;
	
	my $extracted_table = $self->{extracted_table};
	my @data;
	my @fields = qw [ organism assigned_name assigned_gene ];
	push (@fields,  "count(*) AS 'number'");
	my $where = "GROUP BY  Organism, Assigned_name 
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
		print "\n\t #  $number matches to:  ";
		print " $assigned_name, $assigned_gene";
		print " in $organism \t";
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
		print "\n\t- - - PAUSING 5 seconds to allow for cancel - - -\n";
		sleep 5;
		
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
		print "\n\t- - - PAUSING 3 seconds to allow for cancel - - -\n";
		sleep 3;

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

############################################################################
# EOF
############################################################################
