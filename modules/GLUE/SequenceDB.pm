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
	
		# Paths and constants	
		refseq_use_path       => $parameter_ref->{refseq_use_path},
		db_filter_path        => $parameter_ref->{db_filter_path},
		
		# Data structures
		db_filters            => $parameter_ref->{db_filters},
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
	my $seqs = $self->summarise_sequence_table();	
	
	# Summarise location table
	$self->summarise_location_table();	
	
	# Summarise genotype table
	if ($seqs) {
		$self->summarise_genotype_table();	
	}
	
	# Summarise mutation table
	if ($seqs) {
		$self->summarise_mutation_table();
	}	
}

#***************************************************************************
# Subroutine:  summarise mutation table
# Description: 
#***************************************************************************
sub summarise_mutation_table {

	my ($self) = @_;
	
	# Create typical list
	print "\n\n\t # Summarizing the mutation table";
	#my $question = "\n\n\t Create a typical list";
	my $question = "\n\n\t Enter a cutoff for typical:";
	my $upper  = 0;
	my $lower  = 1;
	my $cutoff = $console->ask_float_with_bounds_question($question, $upper, $lower);
	#print $cutoff;

	# Get position data (what is the coverage at all positions for each refseq and gene)
	my $position_table = $self->{position_table};
	unless ($position_table) { die;}
	my @pos_fields = qw [ gene_id position pos_count ];
	my @position_data;
	my %position_index;
	$position_table->select_rows(\@pos_fields, \@position_data);
	foreach my $data_ref (@position_data) {
		my $gene_id   = $data_ref->{gene_id};
		my $position  = $data_ref->{position};
		my $pos_count = $data_ref->{pos_count};
		#$devtools->print_hash($data_ref); die;
		my $position_key = $gene_id . ':' . $position;
		$position_index{$position_key} = $pos_count;
	}

	# Get mutation data (what mutations are present and what in what number of sequences)
	my $mutation_table = $self->{mutation_table};
	unless ($mutation_table) { die;}
	my @data;
	my @fields = qw [ gene_id reference position mutation ];
	push (@fields,  "count(*) AS 'number'");
	my $where = " GROUP BY Gene_ID, Reference, Position, Mutation
                  ORDER BY Gene_ID, Position; ";
	$mutation_table->select_rows(\@fields, \@data, $where);
	my $count = 0;
	foreach my $data_ref (@data) {
		my $number = $data_ref->{number};
		$count = $count + $number;
	}
	print "\n\n\t # The Mutation table contains a total of '$count' rows";
	print "\n\t #  ---";
	#$devtools->print_hash(\%position_index); die;
	my @typical;
	foreach my $data_ref (@data) { # Get the data	
		my $gene_id = $data_ref->{gene_id};
		my $reference = $data_ref->{reference};
		my $position = $data_ref->{position};
		my $mutation = $data_ref->{mutation};
		my $number   = $data_ref->{number};
		my $position_key = $gene_id . ':' . $position;
		my $pos_count = $position_index{$position_key};
		my $frequency = $number / $pos_count;
		my $f_frequency = sprintf("%.3f", $frequency);
		print "\n\t #  $number of $pos_count ($f_frequency) mutation";
		print "$reference $position $mutation in $gene_id";
		my $line = "typical\tNULL\t$gene_id\t$reference\t$position\t$number\t$f_frequency\n";
		push (@typical, $line);
	}
	$fileio->write_file('typical_list.txt', \@typical);
	sleep 1;
}

#***************************************************************************
# Subroutine:  summarise sequence table
# Description: 
#***************************************************************************
sub summarise_sequence_table {

	my ($self, $done_ref) = @_;
	
	my $sequence_table = $self->{sequence_table};
	unless ($sequence_table) { die;}

	my @data;
	my @fields;
	push (@fields,  "count(*) AS 'number'");
	$sequence_table->select_rows(\@fields, \@data);
	#$devtools->print_array(\@data); die;
	my $data_ref = shift @data;
	my $number = $data_ref->{number};
	if ($number) {
		print "\n\n\t # The Sequence table contains a total of '$number' rows";
	}
	return $number;
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

#***************************************************************************
# Subroutine:  summarise genotype table
# Description: 
#***************************************************************************
sub summarise_genotype_table {

	my ($self, $done_ref) = @_;
	
	my $genotype_table = $self->{genotype_table};
	unless ($genotype_table) { die;}

	my @data;
	my @fields = qw [ genotype genotype_method ];
	push (@fields,  "count(*) AS 'number'");
	my $where = "GROUP BY Genotype, Genotype_method ORDER BY Genotype_method, COUNT(*) DESC";
	$genotype_table->select_rows(\@fields, \@data, $where);
	my $count = scalar @data;
	print "\n\n\t # The Location table contains a total of '$count' rows";
	print "\n\t #  ---";
	foreach my $data_ref (@data) { # Get the data	
		my $number  = $data_ref->{number};
		my $genotype= $data_ref->{genotype};
		my $genotype_method = $data_ref->{genotype_method};
		print "\n\t #  $number group '$genotype' by '$genotype_method' genotype";
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
	#unless ($dbh) {
    #    print "\n\t ### Creating '$db_name' screening database\n ";
	#	$self->create_sequence_db($db_name);
	#}
	#else {
	#	print "\n\t ### Loading '$db_name' screening database\n ";
	#}
	$self->{dbh} = $dbh;

	# Main Screening DB tables
	$self->load_sequence_table($dbh);	
	$self->load_filter_table($dbh);	
	$self->load_genotype_table($dbh);	
	$self->load_location_table($dbh);	
	$self->load_mutation_table($dbh);	
	$self->load_position_table($dbh);	
	#$self->load_mvariant_table($dbh);	

	# Load filters
	my %filters;
	$self->load_db_filters($db_name, \%filters);
	$self->{db_filters} = \%filters;
	#$devtools->print_hash(\%filters);die;

}

#***************************************************************************
# Subroutine:  load_db_filters
# Description: load_db_filters
#***************************************************************************
sub load_db_filters {

	my ($self, $db_name, $filter_ref) = @_;
	
	my $db_filter_path = $self->{db_filter_path};

	unless ($db_filter_path) {
		print "\n\t No filters found for DB '$db_name'";
		return;
	}

	# Parse the 'SCREENDB' block
	my @filters;
	my $read_success = $fileio->read_file($db_filter_path, \@filters);
	unless ($read_success) {
		die "\n\t No filters found for DB '$db_name'\n\n";
	}
	#$devtools->print_array(\@filters);
	
	my $start = 'BEGIN GENOTYPE';
	my $stop  = 'ENDBLOCK';
	my %translations;
	$fileio->read_standard_field_value_block(\@filters, $start, $stop, \%translations, 1);
	#$devtools->print_hash(\%translations);
	$filter_ref->{genotype} = \%translations;

	$start = 'BEGIN HOST';
	$stop  = 'ENDBLOCK';
	my %host_translations;
	$fileio->read_standard_field_value_block(\@filters, $start, $stop, \%host_translations, 1);
	#$devtools->print_hash(\%host_translations); die;
	$filter_ref->{host} = \%host_translations;

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
		isolate_id           => 'varchar',
		isolate_source       => 'varchar',
		host                 => 'varchar',
		sequence_date        => 'text',
		sequence_length      => 'int',
		start                => 'int',
		stop                 => 'int',
		isolation_date       => 'varchar',
		isolation_date_class => 'varchar',
		definition           => 'text',
		sequence             => 'text',
	);
	my $sequence_table = MySQLtable->new('Sequence', $dbh, \%sequence_fields);
	$self->{sequence_table} = $sequence_table;
}

#***************************************************************************
# Subroutine:  load_filter_table
# Description: load filter table
#***************************************************************************
sub load_filter_table {

	my ($self, $dbh) = @_;
	
	# Definition of the table
	my %filter_fields = (
		sequence_id          => 'varchar',
		filter_rule          => 'varchar',
		rule_category        => 'varchar',
	);
	my $filter_table = MySQLtable->new('Filter', $dbh, \%filter_fields);
	$self->{filter_table} = $filter_table;
}

#***************************************************************************
# Subroutine:  load_position_table
# Description: load position table
#***************************************************************************
sub load_position_table {

	my ($self, $dbh) = @_;
	
	# Definition of the table
	my %position_fields = (
		refseq_id        => 'varchar',
		gene_id          => 'varchar',
		position         => 'int',
		pos_count        => 'int',
	);
	my $position_table = MySQLtable->new('Position', $dbh, \%position_fields);
	$self->{position_table} = $position_table;
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
	$self->create_filter_table($dbh);
	$self->create_genotype_table($dbh);
	$self->create_location_table($dbh);
	$self->create_mutation_table($dbh);
	$self->create_position_table($dbh);
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
	  `Isolate_ID`      varchar(100) NOT NULL default '0',
	  `Isolate_source`  varchar(100) NOT NULL default '0',
	  `Isolation_date`  varchar(100) NOT NULL default '0',
	  `Isolation_date_class`  varchar(100) NOT NULL default '0',
	  `Host`            varchar(100) NOT NULL default '0',
	  `Start`           int(11) NOT NULL default '0',
	  `Stop`            int(11) NOT NULL default '0',
	  `Sequence_date`   varchar(100) NOT NULL default '0',
	  `Sequence_length` int(11) NOT NULL default '0',
	  `Definition`      text NOT NULL,
	  `Sequence`        text NOT NULL,

	  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($sequence);
	unless ($sth->execute()) { print "\n\t$sequence\n\n\n"; exit;}

}


#***************************************************************************
# Subroutine:  create_filter_table
# Description: create MySQL 'Sequence' table
#***************************************************************************
sub create_filter_table {

	my ($self, $dbh) = @_;

	# Sequence table 
	my $filter = "CREATE TABLE `Filter` (
	  `Record_ID`     int(11) NOT NULL auto_increment,
	  `Sequence_ID`     varchar(100) NOT NULL default '0',
	  `Filter_rule`     varchar(100) NOT NULL default '0',
	  `Rule_category`     varchar(100) NOT NULL default '0',
	  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($filter);
	unless ($sth->execute()) { print "\n\t$filter\n\n\n"; exit;}

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
# Subroutine:  create_position_table
# Description: create sequence database  'Position' table
#***************************************************************************
sub create_position_table {

	my ($self, $dbh) = @_;

	# Position table 
	my $position = "CREATE TABLE `Position` (

	  `Record_ID`     int(11) NOT NULL auto_increment,

	  `Refseq_ID`     varchar(100) NOT NULL default '0',
	  `Gene_ID`       varchar(100) NOT NULL default '0',
	  `Position`      int(11) NOT NULL default '0',
	  `Pos_count`     int(11) NOT NULL default '0',

	  `Timestamp` timestamp NOT NULL default CURRENT_TIMESTAMP on update CURRENT_TIMESTAMP,
	  PRIMARY KEY  (`Record_ID`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
	my $sth = $dbh->prepare($position);
	unless ($sth->execute()) { print "\n\t$position\n\n\n"; exit;}
}

#***************************************************************************
# Subroutine:  create_mutation_table
# Description: create sequence database 'Mutation' table
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
# Description: create sequence database 'Sequence' table
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
		my $filter_table    = $self->{filter_table};
		my $genotype_table  = $self->{genotype_table};
		my $location_table  = $self->{location_table};
		my $mutation_table  = $self->{mutation_table};
		my $position_table  = $self->{position_table};
		#my $mvariant_table  = $self->{mvariant_table};
		
		# Flush result tables
		$sequence_table->flush();
		$sequence_table->reset_primary_keys();
		$filter_table->flush();
		$filter_table->reset_primary_keys();
		$genotype_table->flush();
		$genotype_table->reset_primary_keys();
		$mutation_table->flush();
		$mutation_table->reset_primary_keys();
		$position_table->flush();
		$position_table->reset_primary_keys();
		$location_table->flush();
		$location_table->reset_primary_keys();
		#$mvariant_table->flush();
		#$mvariant_table->reset_primary_keys();
	}
}

############################################################################
# LOADING EXCLUDE SET
############################################################################

#***************************************************************************
# Subroutine:  load_exclude_set 
# Description: load exclude set
#***************************************************************************
sub load_exclude_set {

	my ($self, $file) = @_;

	# Read the file
	unless ($file) { die; }
	my @exclude;
	$fileio->read_file($file, \@exclude);
	
	# Get the filter rule for excluding these
	my $filter = 'seqlen_short';	
	
	# Get further information
	my $category = '500nt';	

	# Index the filter table
	my $filter_table    = $self->{filter_table};
	my @filter;
	my @fields = qw [ sequence_id ];
	$filter_table->select_rows(\@fields, \@filter);
	my %exclude_ids;
	foreach my $data_ref (@filter) {
		my $filter_id = $data_ref->{sequence_id};
		$exclude_ids{$filter_id} = 1;	
	}

	# Upload info
	foreach my $line (@exclude) {

		my @line = split("\t", $line);
		my $id = shift @line;
		if ($exclude_ids{$id}) { next; }
		my %data;
		$data{sequence_id} = $id;
		$data{filter_rule} = $filter;
		$data{rule_category} = $category;
		$filter_table->insert_row(\%data);
	}
}

#***************************************************************************
# Subroutine:  generate_typical_list_from_db 
# Description: 
#***************************************************************************
sub generate_typical_list_from_db {

	my ($self, $refseq, $threshold) = @_;
	
	# HACK TO DO: fix this, no hardcoding
	$threshold    = 0.01;
	my $listname  = 'typical';
	my $version   = '1.0';
	my $color     = 'black';

	# Get tables
	my $position_table = $self->{position_table};
	my $mutation_table = $self->{mutation_table};

	# Get the reference amino acids for each RefSeq ORF
	my $genes_ref = $refseq->{hash_genes};

	# Get all the position information
	my @fields = qw [ gene_id position pos_count ];
	my @data;
	my $refseq_name = $refseq->{name};
	my $where = " WHERE Refseq_ID = '$refseq_name'";
	$position_table->select_rows(\@fields, \@data, $where);
	my %position_counts;
	foreach my $row_ref (@data) {
		my $gene_id   = $row_ref->{gene_id};
		my $position  = $row_ref->{position};
		my $pos_count = $row_ref->{pos_count};
		my $pos_key = $gene_id . ':' . $position;	
		$position_counts{$pos_key} = $pos_count;
	}

	# Get counts for all mutations
	my @mut_fields = qw [ gene_id position mutation ];
	my $count = 'COUNT(*)';
	push (@mut_fields, $count);
	my @mut_data;
	my $order = 'GROUP BY gene_id, position, mutation
                 ORDER BY gene_id, position, mutation';
	$mutation_table->select_distinct(\@mut_fields, \@mut_data, $order);

	my @list_output;
	push (@list_output, "Begin Mutlist;\n");
	push (@list_output, "name=$listname;\n");
	push (@list_output, "color=$color;\n");

	my $hash_genes_ref = $refseq->{hash_genes};
	foreach my $mut_data_ref (@mut_data) {
		my $gene_id  = $mut_data_ref->{gene_id};
		my $position = $mut_data_ref->{position};
		my $mutation = $mut_data_ref->{mutation};
		my $number  = $mut_data_ref->{'COUNT(*)'};
		my $pos_key = $gene_id . ':' . $position;
		my $mut_key = $gene_id . ':' . $position . ':' . $mutation;
		#$devtools->print_hash($mut_data_ref); die;

		my $gene_ref     = $hash_genes_ref->{$gene_id};
		my $indexed_aas  = $gene_ref->{indexed_aas};
		my $reference_aa = $indexed_aas->{$position}; 
		#$devtools->print_hash($gene_ref); die;
		
		my $pos_count = $position_counts{$pos_key};
		unless ($pos_count) { die; }
		my $mut_freq = $number / $pos_count;
		my $f_freq = sprintf("%.3f", $mut_freq);
		#print "\n\t Mutation $mut_key: frequency: $f_freq";		
		if ($f_freq >= $threshold) {
			my @row;
			push (@row, $listname);
			push (@row, $version);
			push (@row, $gene_id);
			push (@row, $reference_aa);
			push (@row, $position);
			push (@row, $mutation);
			push (@row, $f_freq);
			my $row = join ("\t", @row);
			push (@list_output, "$row\n");
			
		}	
	}
	
	push (@list_output, "Endblock;\n");
	my $outfile_name = $refseq_name . '_typical_list.txt'; 
	$fileio->write_file($outfile_name, \@list_output);
}

############################################################################
# EOF
############################################################################
