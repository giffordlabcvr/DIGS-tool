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
		genome_use_path        => $digs_obj->{genome_use_path},
		db_name                => '',   # Obtained from control file or user
		mysql_server           => '',   # Obtained from control file or user

	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# TOP LEVEL HANDLER
############################################################################

#***************************************************************************
# Subroutine:  run_utility_process
# Description: top-level handler for DIGS tool utility functions 
#***************************************************************************
sub run_utility_process {

	my ($self, $infile, $database, $genomes, $utility) = @_;

 	# Show title
	my $digs_obj = $self->{digs_obj};
	$digs_obj->show_title();  

	if ($genomes) {	# Run a target database summary process
		$self->run_target_db_summary_process($genomes);		
	}
	elsif ($database) { # Run a target database summary process

		if ($database eq 6) {
			$self->extract_track_sequences($infile);
		}
		else {

			# Use a control file to connect to database
			if ($infile) {
				$self->parse_ctl_file_and_connect_to_db($infile);
			}
			else {
				print "\n\t  No infile was supplied: enter parameters below\n\n";
				$self->do_load_db_dialogue();		
			}
			if ($database) {
				$self->run_screening_db_management_process($database);
			}
		}
	}
	else { 
		die;
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
		$self->append_to_digs_results($ctl_file);
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

#***************************************************************************
# Subroutine:  append_to_digs_results
# Description: append to digs_result table
#***************************************************************************
sub append_to_digs_results {

	my ($self) = @_;

	# Get database handle, die if we can't 
	my $digs_obj = $self->{digs_obj};

	my $db = $digs_obj->{db};
	unless ($db) { die; }

	# Get the file path
	print "\n\n\t #### WARNING: This function expects a tab-delimited data table with column headers!";
	my $question = "\n\n\t Please enter the path to the file with the table data and column headings\n\n\t";
	my $infile = $console->ask_question($question);
	unless ($infile) { die; }

	# Insert the data	
	print "\n\n\t #### IMPORTING '$infile' to digs_results table";
	$db->import_data_to_digs_results($infile)

}

############################################################################
# DATA FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  extract_track_sequences
# Description: extract FASTA nucs from a genome assembly using an input track 
#***************************************************************************
sub extract_track_sequences {
	
	my ($self, $infile) = @_;

	# Index genomes by key ( organism | type | version )
	my $digs_obj = $self->{digs_obj};
    my $genome_use_path  = $digs_obj->{genome_use_path};
	unless ($genome_use_path) { die "\n\t Path to genomes is not set\n\n\n"; }
	my $target_db_obj = TargetDB->new($digs_obj);
	my %server_data;
	$target_db_obj->read_target_directory(\%server_data);
	#$devtools->print_hash(\%server_data);
	#die;

	# Read infile
	my @infile;
	$fileio->read_file($infile, \@infile);

	# Get the header row
	my $header_row = shift @infile;
	my %fields;
	my %params;
    $self->check_input_and_set_flank_params($header_row, \%params, \%fields);

    # Iterate through the tracks extracting
	my @fasta;	
	foreach my $line (@infile) {
		
		chomp $line; # remove newline
		#print "\n\t ## LINE: $line";
        
        # Index the line properties, using field_names as keys
        my %indexed;
        $self->create_indexed_line($line, $header_row, \%fields, \%indexed);
         
        # Get the correct path to the indexed target file
		my $target_path = $self->create_target_path(\%server_data, \%indexed);
		$params{target_path} = $target_path;
		
        # Now do the extraction
		if ($target_path) {
        	my $digs_fasta =  $self->extract_seq_and_flanks(\%params, \%indexed);
			push (@fasta, $digs_fasta);
		}
	}

	my $outfile = 'extracted.DIGS.fna';
	$fileio->write_file($outfile, \@fasta);

}

	
#***************************************************************************
# Subroutine:  check_input_and_set_flank_params
# Description: 
#***************************************************************************
sub check_input_and_set_flank_params {

	my ($self, $header_row, $params_ref, $fields_ref) = @_; 

	my @header_row = split ("\t", $header_row);		
	print "\n\n\t The following column headers (i.e. table fields) were obtained\n";
	my $i = '0';
	foreach my $element (@header_row) {
		chomp $element;
		$element =~ s/\s+/_/g;
		print "\n\t\t Column $i: '$element'";
		$fields_ref->{$element} = $i;
		$i++;
	}
	print "\n\n\t\t CHECK COLUMN HEADINGS LOOK CORRECT!\n"; sleep 1; 
    
    my $question3 = "\n\t  Set 3' flank";
    my $question5 = "\n\t  Set 5' flank";
	my $flank3  = $console->ask_int_with_bounds_question($question3, '0', 100000);
	my $flank5  = $console->ask_int_with_bounds_question($question5, '0', 100000);
    $params_ref->{flank3} = $flank3;
    $params_ref->{flank5} = $flank5;

}
	
#***************************************************************************
# Subroutine:  extract_seq_and_flanks
# Description: 
#***************************************************************************
sub extract_seq_and_flanks {

	my ($self, $params_ref, $indexed_ref) = @_;

	my $blast_obj      = $self->{blast_obj};
	my $target_path    = $params_ref->{'target_path'};
	my $flank3         = $params_ref->{'flank3'};
	my $flank5         = $params_ref->{'flank5'};
	my $organism       = $indexed_ref->{'organism'};
	my $scaffold       = $indexed_ref->{'scaffold'};
	my $orientation    = $indexed_ref->{'orientation'};
	my $start          = $indexed_ref->{'extract_start'};
	my $end            = $indexed_ref->{'extract_end'};	
	my $record_id      = $indexed_ref->{'record_ID'};	
	my $virus_genus    = $indexed_ref->{'virus_genus'};	
  
    my $truncated3;
    my $truncated5;
    my $extract_start;
    my $extract_end;
    if ($orientation eq '+') {
       my $adjust_start = $start - $flank5;
       if ($adjust_start < 1) {
	      print  "\n\t\t\t 5' TRUNCATED!!!";
          $extract_start = 1;
          $truncated5 = abs($adjust_start);
	   }
       else {
          $extract_start = $adjust_start;
       }
       $extract_end = $end + $flank3;
	}
    elsif ($orientation eq '-') {
       my $adjust_start = $start - $flank3;
       if ($adjust_start < 1) {
	      print  "\n\t\t\t 5' TRUNCATED!!!";
          $extract_start = 1;
          $truncated3 = abs($adjust_start);
	   }
       else {
          $extract_start = $adjust_start;
       }
       $extract_end = $end + $flank5;
	}
 
	# Adjust the start coordinates according to NYT  
	my %data;
	$data{scaffold}    = $scaffold; 
	$data{orientation} = $orientation; # $orientation;	
	$data{start}       = $extract_start;
	$data{end}         = $extract_end;
	my $digs_fasta = '';

    my $expected_length = ($extract_end - $extract_start) + 1;
	my $sequence = $blast_obj->extract_sequence($target_path, \%data);
    my $extracted_length = length $sequence;
    if ($expected_length > $extracted_length) {
       
	   print  "\n\t\t\t 3' TRUNCATED!!!";
	   #print  "\n\t\t\t Expected length: '$expected_length', extracted length: '$extracted_length'";
       my $length_missing = ($expected_length - $extracted_length) + 1;
       $truncated5 = abs($length_missing);
	}

    unless ($sequence) {	
		print  "\n\t Sequence extraction failed for record";
		sleep 1;
	}
	else {
		my $header = $organism . ',' . $scaffold . ',' . $extract_start . ',' . $extracted_length;
		if ($truncated5) {
           $header .= "-T5($truncated5)";
        }
		if ($truncated3) {
           $header .= "-T5(0)-T3($truncated3)" . $extracted_length;
        }
		print "\n\t\t Got sequence for $header";
		$digs_fasta = ">$header" . "\n$sequence\n";
	}
	return $digs_fasta;

}

#***************************************************************************
# Subroutine:  create_indexed_line
# Description: 
#***************************************************************************
sub create_indexed_line {

	my ($self, $property_row, $field_names_row, $field_name_index, $indexed_ref) = @_;

	my @field_names_row = split ("\t", $field_names_row);		
	my $field_i = '0';
	my @line = split("\t", $property_row);
	foreach my $field (@field_names_row) {
		chomp $field;
		$field =~ s/\s+/_/g;
		my $property_i = $field_name_index->{$field};
		my $property = $line[$property_i];
		$indexed_ref->{$field} = $property;
		$field_i++;
		#print "\n\t ### Index: $field_i | $property_i: FIELD '$field' = '$property'";
	}
}

#***************************************************************************
# Subroutine:  create_target_path
# Description: 
#***************************************************************************
sub create_target_path {

	my ($self, $server_data, $indexed_ref) = @_;

	my $digs_obj = $self->{digs_obj};
    my $genome_use_path  = $digs_obj->{genome_use_path};
	unless ($genome_use_path) { die "\n\t Path to genomes is not set\n\n\n"; }
	my $target_path = '';
	
    #$devtools->print_hash($genome_data);
	my $organism       = $indexed_ref->{'organism'};
	my $type           = $indexed_ref->{'target_datatype'};
	my $version        = $indexed_ref->{'target_version'};
	my $target_name    = $indexed_ref->{'target_name'};
	my $scaffold       = $indexed_ref->{'scaffold'};
	my $orientation    = $indexed_ref->{'orientation'};

    # Get the top level categorisation
	my $key = $organism . '|' . $type . '|' . $version; # Create the key
	my $genome_data = $server_data->{$key};
	unless ($genome_data) { 
		print "\n\t\t\t Skipping row '$key': missing link data";
		sleep 1;
		return $target_path;
	}
	my $group = $genome_data->{grouping};
	unless ($group) { 
		print "\n\t\t\t Skipping row - missing genome path data (group)";
		return $target_path;
	}
	#print "\n\t TARGET $organism - grouping $group";
    my @target_path;
	push (@target_path, $genome_use_path);
	push (@target_path, $group);
	push (@target_path, $organism);
	push (@target_path, $type);
	push (@target_path, $version);
	push (@target_path, $target_name);
	$target_path = join('/', @target_path);

	return $target_path;

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
