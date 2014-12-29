#!/usr/bin/perl -w
############################################################################
# Module:      SequenceDB.pm
# Description: The DB module for GLUE
# History:     June 2014: Created by Robert Gifford 
############################################################################
package SequenceDB;

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

#***************************************************************************
# Subroutine:  new
# Description: Parameters
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {
	
		# Paths and constants
		mode                  => $parameter_ref->{mode},
		process_id            => $parameter_ref->{process_id},
		output_type           => $parameter_ref->{output_type},
		
		# DB connection variables
		server                => $parameter_ref->{server},
		username              => $parameter_ref->{username},
		password              => $parameter_ref->{password},
		
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# Public member functions
############################################################################

#***************************************************************************
# Subroutine:  summarise_sequence_db
# Description: summarise a sequence DB
#***************************************************************************
sub summarise_sequence_db {

	my ($self) = @_;

	my $db_name = $self->{db_name};
	
	# Summarise sequence table
	#my $seqs = $self->summarise_sequence_table();	
	my $seqs;
	
	$self->summarise_location_table();	
	
	# Summarise genotype table
	if ($seqs) {
		$self->summarise_genotype_table();	
	}
}

#***************************************************************************
# Subroutine:  summarise Sequence table
# Description: 
#***************************************************************************
sub summarise_sequence_table {

	my ($self, $done_ref) = @_;
	
	my $sequence_table = $self->{sequence_table};
	unless ($sequence_table) { die;}

	my @data;
	my @fields = qw [ location ];
	push (@fields,  "count(*) AS 'number'");
	my $where = "GROUP BY  location 
                 ORDER BY  location ";
	$sequence_table->select_rows(\@fields, \@data, $where);
	my $count = 0;
	foreach my $data_ref (@data) {
		my $number = $data_ref->{number};
		$count = $count + $number;
	}
	print "\n\n\t # The Sequence table contains a total of '$count' rows";
	print "\n\t #  ---";
	foreach my $data_ref (@data) { # Get the data	
		my $country = $data_ref->{country};
		my $number  = $data_ref->{number};
		print "\n\t #  $number sequence from '$country'";
	}
	sleep 1;
	die;
}

#***************************************************************************
# Subroutine:  summarise location table
# Description: 
#***************************************************************************
sub summarise_location_table {

	my ($self, $done_ref) = @_;
	
	my $location_table = $self->{location_table};
	unless ($location_table) { die;}

	my @data;
	my @fields = qw [ location_id ];
	push (@fields,  "count(*) AS 'number'");
	my $where = "WHERE location_type = 'country'
                 GROUP BY  location_id 
                 ORDER BY  location_id ";
	$location_table->select_rows(\@fields, \@data, $where);
	my $count = 0;
	foreach my $data_ref (@data) {
		my $number = $data_ref->{number};
		$count = $count + $number;
	}
	print "\n\n\t # The Location table contains a total of '$count' rows";
	print "\n\t #  ---";
	foreach my $data_ref (@data) { # Get the data	
		my $country = $data_ref->{location_id};
		my $number  = $data_ref->{number};
		print "\n\t #  $number sequences from:\t '$country'";
	}
	sleep 1;
}

############################################################################
# LOADING SEQUENCE DATABASE
############################################################################

#***************************************************************************
# Subroutine:  load_sequence_db 
# Description: LOAD SEQUENCE DATABASE
#***************************************************************************
sub load_sequence_db {

	my ($self, $db_name) = @_;

	# Get connection variables from self
	my $server   = $self->{server};
	my $username = $self->{username};
	my $password = $self->{password};
	unless ($server and $username and $password and $db_name) {
		$devtools->print_hash($self); die; 
	}

	# Set name
	$self->{db_name} = $db_name;
   	
	# Load tables from the screening result DB we've been given
	my $dbh = DBI->connect("dbi:mysql:$db_name:$server", $username, $password);
	unless ($dbh) {
        print "\n\t ### Creating '$db_name' screening database\n ";
		$self->create_sequence_db($db_name);
	}
	else {
		print "\n\t ### Loading '$db_name' screening database\n ";
	}
	$self->{dbh} = $dbh;

	# Main Screening DB tables
	$self->load_sequence_table($dbh);	
	$self->load_genotype_table($dbh);	
	$self->load_location_table($dbh);	
	$self->load_mutation_table($dbh);	
	#$self->load_mvariant_table($dbh);	
}

#***************************************************************************
# Subroutine:  load_sequence_table
# Description: load sequence table
#***************************************************************************
sub load_sequence_table {

	my ($self, $dbh) = @_;
	
	# Definition of the table
	my %sequence_fields = (
		sequence_id          => 'varchar',
		sample_id            => 'varchar',
		sequence_date        => 'text',
		sequence_length      => 'int',
		start                => 'int',
		stop                 => 'int',
		isolation_date       => 'varchar',
		isolation_date_class => 'varchar',
		sequence             => 'text',
	);
	my $sequence_table = MySQLtable->new('Sequence', $dbh, \%sequence_fields);
	$self->{sequence_table} = $sequence_table;
}

#***************************************************************************
# Subroutine:  load_mutation_table
# Description: load mutation table
#***************************************************************************
sub load_mutation_table {

	my ($self, $dbh) = @_;
	
	# Definition of the table
	my %mutation_fields = (
		sequence_id      => 'varchar',
		gene_id          => 'varchar',
		level            => 'varchar',
		reference        => 'varchar',
		position         => 'int',
		mutation         => 'varchar',
	);
	my $mutation_table = MySQLtable->new('Mutation', $dbh, \%mutation_fields);
	$self->{mutation_table} = $mutation_table;
}

#***************************************************************************
# Subroutine:  load_genotype_table
# Description: load genotype table
#***************************************************************************
sub load_genotype_table {

	my ($self, $dbh) = @_;
	
	# Definition of the table
	my %genotype_fields = (
		sequence_id      => 'varchar',
		genotype         => 'varchar',
		genotype_method  => 'varchar',
		score            => 'varchar',
	);
	my $genotype_table = MySQLtable->new('Genotype', $dbh, \%genotype_fields);
	$self->{genotype_table} = $genotype_table;
}

#***************************************************************************
# Subroutine:  load_location_table
# Description: load genotype table
#***************************************************************************
sub load_location_table {

	my ($self, $dbh) = @_;
	
	# Definition of the table
	my %location_fields = (
		sequence_id      => 'varchar',
		location_id      => 'varchar',
		location_type    => 'varchar',
	);
	my $location_table = MySQLtable->new('Location', $dbh, \%location_fields);
	$self->{location_table} = $location_table;
}

############################################################################
# CREATE SEQUENCE DATABASE
############################################################################

#***************************************************************************
# Subroutine:  create_sequence_db
# Description: 
#***************************************************************************
sub create_sequence_db {

	my ($self, $db_name) = @_;

	# Get connection variables from self
	my $server   = $self->{server};
	my $username = $self->{username};
	my $password = $self->{password};
	unless ($server and $username and $password and $db_name) { die; }

	# CREATE THE DB 
	print "\n\n\t ### Creating '$db_name' screening database";
    my $drh = DBI->install_driver("mysql");
    my $rc = $drh->func("createdb", $db_name, $server, $username, $password, 'admin');
    my $dbh = DBI->connect("dbi:mysql:$db_name:$server", $username, $password);
	unless ($dbh) {
		die "\n\t # Couldn't create screening database '$db_name'\n\n";
	}	

	$self->create_sequence_table($dbh);
	$self->create_genotype_table($dbh);
	$self->create_location_table($dbh);
	$self->create_mutation_table($dbh);
	#$self->create_mvariant_table($dbh);
}

#***************************************************************************
# Subroutine:  create_sequence_table
# Description: create MySQL 'Sequence' table
#***************************************************************************
sub create_sequence_table {

	my ($self, $dbh) = @_;

	# Sequence table 
	my $sequence = "CREATE TABLE `Sequence` (

	  `Record_ID`     int(11) NOT NULL auto_increment,

	  `Sequence_ID`     varchar(100) NOT NULL default '0',
	  `Sample_ID`       varchar(100) NOT NULL default '0',
	  `Isolation_date` varchar(100) NOT NULL default '0',
	  `Isolation_date_class`  varchar(100) NOT NULL default '0',
	  `Start`           int(11) NOT NULL default '0',
	  `Stop`            int(11) NOT NULL default '0',
	  `Sequence_date`   varchar(100) NOT NULL default '0',
	  `Sequence_length` int(11) NOT NULL default '0',
	  `Sequence`        text NOT NULL,

	  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($sequence);
	unless ($sth->execute()) { print "\n\t$sequence\n\n\n"; exit;}

}

#***************************************************************************
# Subroutine:  create_genotype_table
# Description: create MySQL 'Sequence' table
#***************************************************************************
sub create_genotype_table {

	my ($self, $dbh) = @_;

	# Genotype table 
	my $genotype = "CREATE TABLE `Genotype` (

	  `Record_ID`     int(11) NOT NULL auto_increment,

	  `Sequence_ID`     varchar(100) NOT NULL default '0',
	  `Genotype`        varchar(100) NOT NULL default '0',
	  `Genotype_method` varchar(100) NOT NULL default '0',
	  `Score`           varchar(100) NOT NULL default '0',

	  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($genotype);
	unless ($sth->execute()) { print "\n\t$genotype\n\n\n"; exit;}
}

#***************************************************************************
# Subroutine:  create_location_table
# Description: create MySQL 'Sequence' table
#***************************************************************************
sub create_location_table {

	my ($self, $dbh) = @_;

	# Genotype table 
	my $location = "CREATE TABLE `Location` (

	  `Record_ID`     int(11) NOT NULL auto_increment,

	  `Sequence_ID`   varchar(100) NOT NULL default '0',
	  `Location_ID`   varchar(100) NOT NULL default '0',
	  `Location_type` varchar(100) NOT NULL default '0',

	  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($location);
	unless ($sth->execute()) { print "\n\t$location\n\n\n"; exit;}
}

#***************************************************************************
# Subroutine:  create_mutation_table
# Description: create MySQL 'Sequence' table
#***************************************************************************
sub create_mutation_table {

	my ($self, $dbh) = @_;

	# Mutation table 
	my $mutation = "CREATE TABLE `Mutation` (

	  `Record_ID`     int(11) NOT NULL auto_increment,

	  `Sequence_ID`   varchar(100) NOT NULL default '0',
	  `Gene_ID`       varchar(100) NOT NULL default '0',
	  `Level`         varchar(100) NOT NULL default '0',
	  `Reference`     varchar(100) NOT NULL default '0',
	  `Position`      int(11) NOT NULL default '0',
	  `Mutation`      varchar(100) NOT NULL default '0',

	  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($mutation);
	unless ($sth->execute()) { print "\n\t$mutation\n\n\n"; exit;}
}

#***************************************************************************
# Subroutine:  create_mvariant_table
# Description: create MySQL 'Sequence' table
#***************************************************************************
sub create_mvariant_table {

	my ($self, $dbh) = @_;

	# Minority variant table 
	my $mvariant = "CREATE TABLE `MVariant` (

	  `Record_ID`     int(11) NOT NULL auto_increment,

	  `Sequence_ID`   varchar(100) NOT NULL default '0',
	  `Gene_ID`       varchar(100) NOT NULL default '0',
	  `Position`      int(11) NOT NULL default '0',
	  `Fraction`      varchar(100) NOT NULL default '0',
	  `Read_depth`    int(11) NOT NULL default '0',
	  `Ref_NT`        varchar(100) NOT NULL default '0',
	  `Ref_AA`        varchar(100) NOT NULL default '0',
	  `Ref_codon`     varchar(100) NOT NULL default '0',
	  `Codon`         varchar(100) NOT NULL default '0',
	  `NT`            varchar(100) NOT NULL default '0',
	  `AA`            varchar(100) NOT NULL default '0',

	  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($mvariant);
	unless ($sth->execute()) { print "\n\t$mvariant\n\n\n"; exit;}
}

############################################################################
# DROPPING AND FLUSHING DBs
############################################################################

#***************************************************************************
# Subroutine:  drop_sequence_db
# Description: drop a screening database from the server
#***************************************************************************
sub drop_sequence_db {

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
# Subroutine:  flush_sequence_db 
# Description: flush screening results database tables
#***************************************************************************
sub flush_sequence_db {

	my ($self) = @_;
	
	# Ask to make sure
	my $db_name = $self->{db_name};
	my $question = "\n\n\t Are you sure you want to flush data in the $db_name database?";
	my $answer1 = $console->ask_yes_no_question($question);
	if ($answer1 eq 'y') {

		# get tables
		my $sequence_table  = $self->{sequence_table};
		my $genotype_table  = $self->{genotype_table};
		my $location_table  = $self->{location_table};
		my $mutation_table  = $self->{mutation_table};
		my $mvariant_table  = $self->{mvariant_table};
		
		# Flush result tables
		$sequence_table->flush();
		$sequence_table->reset_primary_keys();
		$genotype_table->flush();
		$genotype_table->reset_primary_keys();
		$mutation_table->flush();
		$mutation_table->reset_primary_keys();
		$location_table->flush();
		$location_table->reset_primary_keys();
		#$mvariant_table->flush();
		#$mvariant_table->reset_primary_keys();
	}
}

############################################################################
# EOF
############################################################################
