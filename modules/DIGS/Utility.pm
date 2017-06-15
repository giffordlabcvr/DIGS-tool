#!usr/bin/perl -w
############################################################################
# Module:      Utility.pm   
# Description: DIGS tool utility functions
# History:     December  2017: Created by Robert Gifford 
############################################################################
package Utility;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::FileIO;
use Base::Console;
use Base::DevTools;

# Program components
use DIGS::ScreenBuilder; # Functions to set up screen

############################################################################
# Globals
############################################################################

# Base objects
my $fileio    = FileIO->new();
my $console   = Console->new();
my $devtools  = DevTools->new();

# Program components
use DIGS::Initialise;    # Initialises the DIGS tool
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create new Utility 'object'
#***************************************************************************
sub new {

	my ($invocant, $digs_obj) = @_;
	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {

		# Global settings
		process_id             => $digs_obj->{process_id},
		program_version        => $digs_obj->{program_version},
		
		# Member classes 
		digs_obj               => $digs_obj,
		blast_obj              => $digs_obj->{blast_obj},

		# MySQL database connection parameters
		mysql_username         => $digs_obj->{mysql_username}, 
		mysql_password         => $digs_obj->{mysql_password},
		db_name                => '',   # Obtained from control file or user
		mysql_server           => '',   # Obtained from control file or user

	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# HELP INFO PAGE FOR UTILITY OPTIONS
############################################################################

#***************************************************************************
# Subroutine:  show_utility_help_page
# Description: show help page information for utility functions
#***************************************************************************
sub show_utility_help_page {

	my ($self) = @_;

	# Create utility help menu
	$console->refresh();
	my $program_version = $self->{program_version};

    my $HELP   = "\n\t ### DIGS version $program_version - utility functions help menu\n";

       $HELP  .= "\n\n\t ### Summarising target databases\n"; 	   
	   $HELP  .= "\n\t   -g=1   Summarise targets (brief summary, by species)";
	   $HELP  .= "\n\t   -g=2   Summarise targets (long, by individual target file)\n";

       $HELP  .= "\n\t ### Managing DIGS screening DBs\n"; 
	   $HELP  .= "\n\t   -d=1   Import tab-delimited data"; 
	   $HELP  .= "\n\t   -d=2   Flush core tables"; 
	   $HELP  .= "\n\t   -d=3   Drop tables";
	   $HELP  .= "\n\t   -d=4   Drop a screening DB"; 

	   $HELP  .= "\n\n"; 

	print $HELP;
}

############################################################################
# TOP LEVEL HANDLER
############################################################################

#***************************************************************************
# Subroutine:  run_utility_process
# Description: top-level handler for DIGS tool utility functions 
#***************************************************************************
sub run_utility_process {

	my ($self, $infile, $database, $genomes) = @_;

 	# Show title
	my $digs_obj = $self->{digs_obj};
	$digs_obj->show_title();  

	if ($genomes) {	# Run a target database summary process
		$self->run_target_db_summary_process($genomes);		
	}	
	else { # Run a screening database management process

		# Use a control file to connect to database
		if ($infile) {
			$self->parse_ctl_file_and_connect_to_db($infile);
		}
		else {
			$self->do_load_db_dialogue();		
		}
		if ($database) {
			$self->run_screening_db_management_process($database);
		}
	}
}

############################################################################
# SECOND LEVEL HANDLERS FOR UTILITY FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  run_target_db_summary_process
# Description: handler for functions summarising target databases
#***************************************************************************
sub run_target_db_summary_process {

	my ($self, $option) = @_;

	my $digs_obj = $self->{digs_obj};
	
	if ($option eq 1)    { # Summarise target genome directory (short)
		my $target_db_obj = TargetDB->new($digs_obj);
		$target_db_obj->summarise_targets_short();
	}
	elsif ($option eq 2) { # Summarise target genome directory (long)
		my $target_db_obj = TargetDB->new($digs_obj);
		$target_db_obj->summarise_targets_long();
	}
	else {
		print "\n\t  Unrecognized option '-g=$option'\n";
	}
}

#***************************************************************************
# Subroutine:  run_screening_db_management_process
# Description: handler for DIGS screening database utility fxns 
#***************************************************************************
sub run_screening_db_management_process {

	my ($self, $option, $ctl_file) = @_;

	my $digs_obj = $self->{digs_obj};

	# Hand off to functions 
	if ($option eq 1) { # Manage ancillary tables in a screening DB
		$self->import_data();
	}
	elsif ($option eq 2) { # Flush screening DB
		my $db = $digs_obj->{db};
		my $db_name = $db->{db_name};
		my $question = "\n\n\t  Are you sure you want to flush data in the $db_name database?";
		my $answer1 = $console->ask_yes_no_question($question); # Ask to make sure
		if ($answer1 eq 'y') { $db->flush_screening_db(); }
	}
	elsif ($option eq 3) {	# Drop tables
		my $db = $digs_obj->{db};
		my $db_name = $db->{db_name};
		$self->drop_db_table();
	}
	elsif ($option eq 4) { # Drop screening DB 
		my $db = $digs_obj->{db};
		$db->drop_screening_db();    
	}
	elsif ($option eq 5) {
		$self->extract_track_sequences($ctl_file);
	}
	else {
		print "\n\t  Unrecognized option '-d=$option'\n";
	}
}

############################################################################
# DATABASE MANAGEMENT FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  import_data
# Description: import data to DIGS screening DB tables
#***************************************************************************
sub import_data {

	my ($self) = @_;

	# Get database handle, die if we can't 
	my $digs_obj = $self->{digs_obj};

	my $db = $digs_obj->{db};
	unless ($db) { die; }
	my $dbh = $db->{dbh};
	unless ($dbh) { die "\n\t Couldn't retrieve database handle \n\n"; }

	# Get the file path
	print "\n\n\t #### WARNING: This function expects a tab-delimited data table with column headers!";
	my $question = "\n\n\t Please enter the path to the file with the table data and column headings\n\n\t";
	my $infile = $console->ask_question($question);
	unless ($infile) { die; }

	# Read in the data from a tab delimited file
	my @data;
	my %fields;
	my @fields;
	$console->do_read_tabdelim_dialogue($infile, \@data, \@fields, \%fields);
	#$devtools->print_array(\@data); die;

	# Get a reference to a table object for the ancillary table
	my $anc_table = $db->do_ancillary_table_dialogue(\@fields, \%fields);
	#$devtools->print_hash(\%fields); die;
	#$devtools->print_array(\@fields); die;

	# Insert the data	
	print "\n\n\t #### IMPORTING to table '$anc_table'";
	my $verbose = $self->{verbose}; # Get 'verbose' flag setting
	my $row_count = $db->import_data_to_ancillary_table($anc_table, \@data, \@fields, \%fields, $verbose);

}

#***************************************************************************
# Subroutine:  drop_db_table
# Description: drop a DIGS screening DB table
#***************************************************************************
sub drop_db_table {

	my ($self) = @_;

	# Get database handle, die if we can't 
	my $digs_obj = $self->{digs_obj};

	my $db = $digs_obj->{db};
	my $db_name = $db->{db_name};
	
	unless ($db) { die; }
	my $dbh = $db->{dbh};
	unless ($dbh) { die "\n\t Couldn't retrieve database handle \n\n"; }

	# Get the ancillary tables in this DB
	my @tables;
	my %tables;
	$db->get_ancillary_table_names(\@tables);
	my $num_choices = scalar @tables;
	if ($num_choices) {

		print "\n\n\t  # Drop ancillary tables in DIGS screening DB '$db_name'\n";	
		my $table_num = 0;
		foreach my $table_name (@tables) {
			$table_num++;
			$tables{$table_num} = $table_name;
			print "\n\t\t Table $table_num: '$table_name'";
		}
		my @table_choices = sort keys %tables;
		
		my $question = "\n\n\t Apply to which of the above tables?";
		my $answer   = $console->ask_list_question($question, $num_choices);
		my $table_to_drop = $tables{$answer};
		unless ($table_to_drop) { die; }	

		$db->drop_ancillary_table($table_to_drop);
	}
	else {
		print "\n\n\t  # There are no ancillary tables in DIGS screening DB '$db_name'\n";	
	}
	
}

############################################################################
# DATA FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  extract_track_sequences
# Description: extract FASTA nucs from a genome assembly using an input track 
#***************************************************************************
sub extract_track_sequences {
	
	my ($self, $ctl_file) = @_;

	# Get database handle, die if we can't 
	my $digs_obj = $self->{digs_obj};
	my $db = $digs_obj->{db};
	unless ($db)       { die; }
	unless ($ctl_file) { die; }
	my $dbh = $db->{dbh};
	unless ($dbh) { die "\n\t Couldn't retrieve database handle \n\n"; }

	# Get paths, objects, data structures and variables from self
	my $blast_obj = $digs_obj->{blast_obj};
	my $verbose   = $self->{verbose}; # Get 'verbose' flag setting

	# Try opening control file
	my @ctl_file;
	my $valid = $fileio->read_file($ctl_file, \@ctl_file);
	unless ($valid) {  # Exit if we can't open the file
		die "\n\t ### Couldn't open control file '$ctl_file'\n\n\n ";
	}

	# If control file looks OK, store the path and parse the file
	$self->{ctl_file} = $ctl_file;
	my $loader_obj = ScreenBuilder->new($digs_obj);
	$loader_obj->parse_control_file($ctl_file, $digs_obj);

	# Load/create the screening database
	my $db_name = $loader_obj->{db_name};
	unless ($db_name) { die "\n\t Error: no DB name defined \n\n\n"; }
	$digs_obj->initialise_screening_db($db_name);

	# Get all targets
	my $target_db_obj = TargetDB->new($digs_obj);
	my %targets;
	$target_db_obj->read_target_directory(\%targets);
	#$devtools->print_hash(\%targets); die;

	# Try to read the tab-delimited infile
	print "\n\n\t #### WARNING: This function expects a tab-delimited data table with column headers!";
	my $question1 = "\n\n\t Please enter the path to the file with the table data and column headings\n\n\t";
	#my $infile = $console->ask_question($question1);
	my $infile = 'loci.txt';
	unless ($infile) { die; }
	my @infile;
	$fileio->read_file($infile, \@infile);

	# Get the header row
	my $header_row = shift @infile;
	my @header_row = split ("\t", $header_row);		
	print "\n\n\t The following column headers (i.e. table fields) were obtained\n";
	my $i;
	my @fields;
	my %fields;
	foreach my $element (@header_row) {
		chomp $element;
		$i++;
		$element =~ s/\s+/_/g;
		if ($element eq '') { $element = 'EMPTY_COLUMN_' . $i; } 
		print "\n\t\t Column $i: '$element'";
		push (@fields, $element);
		$fields{$element} = "varchar";
	}

	# Prompt user - did we read the file correctly?
	my $question2 = "\n\n\t Is this correct?";
	#my $answer2 = $console->ask_yes_no_question($question2);
	#if ($answer2 eq 'n') { # Exit if theres a problem with the infile
	#	print "\n\t\t Aborted!\n\n\n"; exit;
	#}

	# Iterate through the tracks extracting
	#$devtools->print_hash(\%targets);
	my @fasta;	
	foreach my $line (@infile) {
		
		chomp $line; # remove newline
		my @line = split("\t", $line);
		my $record_id     = shift @line;
		my $organism      = shift @line;
		my $type          = shift @line;
		my $version       = shift @line;
		my $name          = shift @line;
        my $scaffold      = shift @line;
        my $orientation   = shift @line;
        my $assigned_name = shift @line;
        my $extract_start = shift @line;
        my $extract_end   = shift @line;
        my $gene          = shift @line;
		
		unless ($organism and $version and $type and $name) { die; }
		
		my $key = $organism . '|' . $type . '|' . $version;
		my $target_data = $targets{$key};
		unless ($target_data) {
			print "\n\t NO DATA FOR GENOME ID '$key'";
		}
		my $group = $target_data->{grouping};
		#print "\n\t GROUP FOR '$key': '$group'";
		
		my $genome_use_path = $digs_obj->{genome_use_path};
		unless ($genome_use_path) {	die; }
		my @target_path;
		push (@target_path, $genome_use_path);
		push (@target_path, $group);
		push (@target_path, $organism);
		push (@target_path, $type);
		push (@target_path, $version);
		push (@target_path, $name);
		my $target_path = join('/', @target_path);
	
		my %data;
		$data{start}       = $extract_start;
        $data{end}         = $extract_end;
        $data{scaffold}    = $scaffold; 
        $data{orientation} = $orientation;	
		#print "\n\t TARGET $target_path";

		my $sequence = $blast_obj->extract_sequence($target_path, \%data);
		unless ($sequence) {	
			print  "\n\t Sequence extraction failed for record '$record_id'";
			sleep 1;
		}
		else {
			my $header = $key . '|' . $scaffold . '|' . $assigned_name . "_$gene";
			#my $header = $name . "_$gene";
			#$header =~ s/\(/\./g;
			#$header =~ s/\)//g;
			print "\n\t\t Got sequence for $record_id: $header";
			my $digs_fasta = ">$header" . "\n$sequence\n";
			push (@fasta, $digs_fasta);;
		}
	}

	my $outfile = 'extracted.DIGS.fna';
	$fileio->write_file($outfile, \@fasta);

}

############################################################################
# INITIALISATION
############################################################################

#***************************************************************************
# Subroutine:  parse_ctl_file_and_connect_to_db
# Description: connect to a DIGS screening DB by parsing a DIGS control file
#***************************************************************************
sub parse_ctl_file_and_connect_to_db {

	my ($self, $infile) = @_;

	my $digs_obj = $self->{digs_obj};
	
	# Try opening control file
	my @ctl_file;
	my $valid = $fileio->read_file($infile, \@ctl_file);
	unless ($valid) {  # Exit if we can't open the file
		die "\n\t ### Couldn't open control file '$infile'\n\n\n ";
	}
	
	# If control file looks OK, store the path and parse the file
	$self->{ctl_file} = $infile;
	my $loader_obj = ScreenBuilder->new($digs_obj);
	$loader_obj->parse_control_file($infile, $digs_obj);

	# Store the ScreenBuilder object (used later)
	$self->{loader_obj} = $loader_obj;

	# Load/create the screening database
	my $db_name = $loader_obj->{db_name};
	unless ($db_name) { die "\n\t Error: no DB name defined \n\n\n"; }
	
	my $initialise_obj = Initialise->new($digs_obj);
	$initialise_obj->initialise_screening_db($digs_obj, $db_name);
}

#***************************************************************************
# Subroutine:  do_load_db_dialogue
# Description: connect to a DIGS screening DB
#***************************************************************************
sub do_load_db_dialogue {

	my ($self, $infile) = @_;

	my $digs_obj = $self->{digs_obj};
	
	# Load/create the screening database
	my $question = "\t  Enter the name of a DIGS screening database";
	my $db_name = $console->ask_question($question);
	unless ($db_name) { die "\n\t Error: no DB name defined \n\n\n"; }
	$digs_obj->{mysql_server}   = 'localhost';

	# Create the screening DB object
	my $db_obj = ScreeningDB->new($digs_obj);

	# Check if this screening DB exists, if not then create it
	my $db_exists = $db_obj->does_db_exist($db_name);
	unless ($db_exists) {
		print "\n\t  Could not connect to screening DB '$db_name'";
		print "\n\t  Exiting.\n\n\n"; exit;	
	}
	
	# Load the database
	print   "\n\t  Connecting to DB:  $db_name";
	$db_obj->load_screening_db($db_name);	
	$digs_obj->{db} = $db_obj; # Store the database object reference 

}

############################################################################
# EOF
############################################################################
