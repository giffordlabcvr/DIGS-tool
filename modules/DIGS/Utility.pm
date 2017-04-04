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
	my $program_version = $self->{program_version};

    my $HELP   = "\n\t ### DIGS version $program_version - utility functions help menu\n";

       $HELP  .= "\n\t ### Managing DIGS screening DBs"; 
	   $HELP  .= "\n\t -d=1   Upload data to tables"; 
	   $HELP  .= "\n\t -d=2   Flush core tables"; 
	   $HELP  .= "\n\t -d=3   Drop DB tables";
	   $HELP  .= "\n\t -d=4   Drop a DIGS screening DB"; 
	   $HELP  .= "\n\t -d=5   Extract data\n";
	   #$HELP  .= "\n\t -d=2   Upload data to the digs_results table";

       $HELP  .= "\n\t ### Creating standard locus nomenclature"; 	   
	   $HELP  .= "\n\t -n=1   Create standard locus IDs\n"; 

       $HELP  .= "\n\t ### Summarizing target databases"; 	   
	   $HELP  .= "\n\t -g=1   Summarise targets (brief summary, by species)";
	   $HELP  .= "\n\t -g=2   Summarise targets (long, by individual target file)\n";

       $HELP  .= "\n\t ### Development and validation tools"; 	   
	   $HELP  .= "\n\t -u=1   Translate DB schema"; 
	   $HELP  .= "\n\t -u=2   Show BLAST chains (merged BLAST hits)"; 
	   $HELP  .= "\n\t -u=3   Show locus chains (merged digs_results table rows)"; 
	   $HELP  .= "\n\t -u=4   Show nomenclature chains (merged annotations)\n"; 

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

	my ($self, $infile, $database, $genomes, $utility) = @_;

 	# Show title
	my $digs_obj = $self->{digs_obj};
	$digs_obj->show_title();  

	if ($genomes) {
		$self->run_target_utility_process($genomes);		
	}
	elsif ($utility) {
		$self->run_data_utility_process($utility, $infile);	
	}
	else {

		# Use a control file to connect to database
		if ($infile) {
			$self->parse_ctl_file_and_connect_to_db($infile);
		}
		else {
			$self->do_load_db_dialogue();		
		}
		if ($database) {
			$self->run_screening_db_utility_process($database);
		}
		elsif ($utility) {
			$self->run_dev_validation_process($utility);		
		}
	}
}

############################################################################
# SECOND LEVEL HANDLERS FOR UTILITY FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  run_screening_db_utility_process
# Description: handler for DIGS screening database utility fxns 
#***************************************************************************
sub run_screening_db_utility_process {

	my ($self, $option, $ctl_file) = @_;

	my $digs_obj = $self->{digs_obj};

	# Hand off to functions 
	if ($option eq 1) { # Manage ancillary tables in a screening DB
		$self->upload_data();
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

#***************************************************************************
# Subroutine:  run_target_utility_process
# Description: handler for functions summarising target databases
#***************************************************************************
sub run_target_utility_process {

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
# Subroutine:  run_dev_validation_process
# Description: handler for DIGS development & validation functions
#***************************************************************************
sub run_dev_validation_process {

	my ($self, $option) = @_;

	my $digs_obj = $self->{digs_obj};
		
	if ($option eq 1) { # DB schema translation
		my $db_obj = $digs_obj->{db};
		$db_obj->translate_schema();
	}
	elsif ($option eq 2) { # Show the BLAST chains for each extracted locus
		$self->show_blast_chains();
	}
	elsif ($option eq 3) { # Show the digs_results breakdown for each consolidated locus
		$self->show_locus_chains();
	}
	elsif ($option eq 4) { # Show the alternative annotations breakdown for each ID-allocated locus
		my $nomenclature_obj = Nomenclature->new($self); 
		$nomenclature_obj->show_nomenclature_chains();
	}
	elsif ($option eq 5) {
		$self->fix_searches_performed_table();
	}
	elsif ($option eq 6) {
		$self->fix_searches_performed_table_2();
	}
	else {
		print "\n\t  Unrecognized option '-u=$option'\n";
	}
}

############################################################################
# DATABASE MANAGEMENT FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  upload_data
# Description: upload data to DIGS screening DB tables
#***************************************************************************
sub upload_data {

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
	#my $infile = $self->ask_question($question);
	my $infile = '../local/human/translations.txt';
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
	unless ($db) { die; }
	my $dbh = $db->{dbh};
	unless ($dbh) { die "\n\t Couldn't retrieve database handle \n\n"; }

	# Get the ancillary tables in this DB
	my @tables;
	my %tables;
	$db->get_ancillary_table_names(\@tables);

	print "\n\n\t  # Drop ancillary tables in DIGS screening DB\n";	
	my $table_num = 0;
	foreach my $table_name (@tables) {
		$table_num++;
		$tables{$table_num} = $table_name;
		print "\n\t\t Table $table_num: '$table_name'";
	}
	my @table_choices = keys %tables;
	my $num_choices = scalar @table_choices;
	
	my $question = "\n\n\t Apply to which of the above tables?";
	my $answer   = $console->ask_list_question($question, $num_choices);
	my $table_to_drop = $tables{$answer};
	unless ($table_to_drop) { die; }	
	$db->drop_ancillary_table($table_to_drop);

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
# DEVELOPMENT & VALIDATION FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  show_blast_chains
# Description: Show BLAST chains for all extracted sequences
#***************************************************************************
sub show_blast_chains {
	
	my ($self) = @_;

	# Get relevant variables and objects
	my $digs_obj = $self->{digs_obj};
	my $db = $digs_obj->{db};
	unless ($db) { die; } # Sanity checking
	my $digs_results_table = $db->{digs_results_table}; 
	my $blast_chains_table = $db->{blast_chains_table};
	my $extract_where = " ORDER BY record_id ";
	my @extracted_ids;
	my @fields = qw [ record_id assigned_name assigned_gene ];
	$digs_results_table->select_rows(\@fields, \@extracted_ids, $extract_where);	 

	# Iterate through the digs result rows	
	foreach my $hit_ref (@extracted_ids) {
		my $digs_result_id = $hit_ref->{record_id};
		my @chain;
		my $assigned_name = $hit_ref->{assigned_name};
		my $assigned_gene = $hit_ref->{assigned_gene};
		my @chain_fields = qw [ record_id probe_name probe_gene 
		                        organism target_name 
		                        scaffold subject_start subject_end
		                        bitscore identity align_len ];
		my $blast_where  = " WHERE digs_result_id = $digs_result_id ";
		   $blast_where .= " ORDER BY subject_start";
		$blast_chains_table->select_rows(\@chain_fields, \@chain, $blast_where);	
		print "\n\t ### BLAST hit chain for extracted locus $digs_result_id";
		print " ($assigned_name, $assigned_gene):";
		foreach my $hit_ref (@chain) {
			my $blast_id    = $hit_ref->{record_id};
			my $probe_name  = $hit_ref->{probe_name};
			my $probe_gene  = $hit_ref->{probe_gene};
			my $organism    = $hit_ref->{organism};
			my $scaffold    = $hit_ref->{scaffold};
			my $start       = $hit_ref->{subject_start};
			my $end         = $hit_ref->{subject_end};
			my $bitscore    = $hit_ref->{bitscore};
			my $identity    = $hit_ref->{identity};
			my $f_identity  = sprintf("%.2f", $identity);
			my $align_len   = $hit_ref->{align_len};
			print "\n\t\t $blast_id:\t Score: $bitscore, \%$f_identity identity ";
			print "across $align_len aa ($start-$end) to:\t $probe_name ($probe_gene) ";
		}
	}
}

#***************************************************************************
# Subroutine:  show_locus_chains
# Description: Show composition of consolidated loci
#***************************************************************************
sub show_locus_chains {
	
	my ($self) = @_;

	# Get relevant variables and objects
	my $digs_obj = $self->{digs_obj};
	my $db_ref = $digs_obj->{db};
	unless ($db_ref) { die; } # Sanity checking
	my $dbh = $db_ref->{dbh};

	my $loci_exists = $db_ref->does_table_exist('loci');
	unless ($loci_exists) {
		print "\n\t  The locus tables don't seem to exist (have you ran consolidate?\n\n";
		exit;
	}

	$db_ref->load_loci_table($dbh);
	$db_ref->load_loci_chains_table($dbh);
	#$devtools->print_hash($db); die;

	my $digs_results_table = $db_ref->{digs_results_table}; 
	my $loci_table         = $db_ref->{loci_table};
	my $loci_chains_table  = $db_ref->{loci_chains_table};
	unless ($digs_results_table and $loci_chains_table) {
		print "\n\t # Locus tables not found - run consolidate first\n\n\n";
		exit;
	}

	# Get all loci
	my $loci_where = " ORDER BY record_id ";
	my @loci;
	my @fields = qw [ record_id locus_structure ];
	$loci_table->select_rows(\@fields, \@loci, $loci_where);	 
	
	# Iterate through loci
	foreach my $locus_ref (@loci) {

		my $locus_id = $locus_ref->{record_id};
		print "\n\t ### Chain $locus_id ";	
		my $chain_where = " WHERE locus_id = $locus_id ";
		my @results;
		my @fields = qw [ record_id locus_id digs_result_id ];
		$loci_chains_table->select_rows(\@fields, \@results, $chain_where);	 
		
		foreach my $result_ref (@results) {

			my @digs_results;
			my $digs_result_id = $result_ref->{digs_result_id};
			my @result_fields = qw [ assigned_name assigned_gene 
			                         organism target_name 
			                         scaffold extract_start extract_end
			                         bitscore identity align_len ];
			my $where  = " WHERE record_id = $digs_result_id ";
			$digs_results_table->select_rows(\@result_fields, \@digs_results, $where);			

			foreach my $result_ref (@digs_results) {
		
				my $assigned_name = $result_ref->{assigned_name};
				my $assigned_gene = $result_ref->{assigned_gene};
				my $organism      = $result_ref->{organism};	
				my $scaffold      = $result_ref->{scaffold};
				my $start         = $result_ref->{extract_start};
				my $end           = $result_ref->{extract_end};
				my $bitscore      = $result_ref->{bitscore};
				my $identity      = $result_ref->{identity};
				my $align_len     = $result_ref->{align_len};
				my $f_identity    = sprintf("%.2f", $identity);
				print "\n\t\t $digs_result_id:\t Score: $bitscore, \%$f_identity identity ";
				print "across $align_len aa ($start-$end) to:\t $assigned_name ($assigned_gene) ";
	
			}
		}
	}
}

#===========================================================================
# TEMPORARY
#===========================================================================

#***************************************************************************
# Subroutine:  fix_searches_performed_table
# Description: 
#***************************************************************************
sub fix_searches_performed_table {

	my ($self) = @_;

	# Get relevant variables and objects
	my $digs_obj = $self->{digs_obj};

	my $db = $digs_obj->{db};
	unless ($db) { die; }
	my $dbh = $db->{dbh};
	unless ($dbh) { die "\n\t Couldn't retrieve database handle \n\n"; }
	my $searches_table = $db->{searches_table}; 
	my $where1 = " ORDER BY record_id ";
	my @searches;
	my @fields = qw [ record_id organism target_datatype target_version ];
	$searches_table->select_rows(\@fields, \@searches, $where1);	 

	# Iterate through the digs result rows	
	foreach my $row_ref (@searches) {
	
		my $record_id = $row_ref->{record_id};
		my $organism  = $row_ref->{organism};
		my $datatype  = $row_ref->{target_datatype};
		my $version   = $row_ref->{target_version};
		my @id = ( $organism, $datatype, $version );
		my $target_id = join ('|', @id);
		my $where2 = " WHERE record_id = $record_id ";
		my %data;
		$data{target_id} = $target_id;
		$searches_table->update(\%data, $where2);
		print "\n\t  UPDATED target_id field for $record_id to '$target_id'";
	}
}

#***************************************************************************
# Subroutine:  fix_searches_performed_table_2
# Description: 
#***************************************************************************
sub fix_searches_performed_table_2 {

	my ($self) = @_;

	# Get relevant variables and objects
	my $digs_obj = $self->{digs_obj};

	my $db = $digs_obj->{db};
	unless ($db) { die; }
	my $dbh = $db->{dbh};
	unless ($dbh) { die "\n\t Couldn't retrieve database handle \n\n"; }
	my $searches_table = $db->{searches_table}; 
	my $where1 = " ORDER BY record_id ";
	my @searches;
	my @fields = qw [ record_id organism target_datatype target_version target_id ];
	$searches_table->select_rows(\@fields, \@searches, $where1);	 

	# Iterate through the digs result rows	
	foreach my $row_ref (@searches) {
	
		my $record_id = $row_ref->{record_id};
		my $datatype  = $row_ref->{target_datatype};
		my $version   = $row_ref->{target_version};
		my $target_id = $row_ref->{target_id};
	
		my @id = split (/\|/, $target_id);
		my $organism  = shift @id;
			
		my $where2 = " WHERE record_id = $record_id ";
		my %data;
		$data{organism} = $organism;
		$searches_table->update(\%data, $where2);
		print "\n\t  UPDATED organism field for $record_id to '$organism'";
	}
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
	$digs_obj->initialise_screening_db($db_name);
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
