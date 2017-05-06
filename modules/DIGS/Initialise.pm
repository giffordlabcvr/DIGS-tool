#!usr/bin/perl -w
############################################################################
# Module:      Initialise.pm
# Description: Functions for initialising DIGS
# History:     April  2017: Created by Robert Gifford 
############################################################################
package Initialise;

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

############################################################################
# Globals
############################################################################

# Base objects
my $fileio    = FileIO->new();
my $console   = Console->new();
my $devtools  = DevTools->new();

# Maximum range for defragment
my $maximum   = 100000000;
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create new DIGS 'object'
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {

		# Global settings
		process_id             => $parameter_ref->{process_id},
		program_version        => $parameter_ref->{program_version},


	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# Main initialisation functions
############################################################################

#***************************************************************************
# Subroutine:  initialise 
# Description: do general set-up, then hand off to option-specific set-up fxns
#***************************************************************************
sub initialise {

	my ($self, $digs_obj, $option, $ctl_file) = @_;

	# DO GENERAL SET-UP 
	my $db_name = $self->do_general_setup($digs_obj, $option, $ctl_file);

	# LOAD/CREATE THE DATABASE
	$self->initialise_screening_db($digs_obj, $db_name);

	# SET-UP FOR DIGS SCREENING
	if ($option eq 2) { 	
		my $valid = $self->setup_for_a_digs_run($digs_obj);
		unless ($valid) { return 0; }
	}
	
	# SET-UP FOR REASSIGN
	if ($option eq 3) { 
		$self->setup_for_reassign($digs_obj);
	}
	
	# DO SET-UP NEEDED FOR BOTH DEFRAGMENT & CONSOLIDATE
	if ($option eq 4 or $option eq 5) { 
		$self->setup_for_defrag_or_consolidate($digs_obj, $option);
	}
	
	return 1;
}

#***************************************************************************
# Subroutine:  do_general_setup
# Description: do general set up for DIGS tool processes
#***************************************************************************
sub do_general_setup {

	my ($self, $digs_obj, $option, $ctl_file) = @_;

	# Try opening control file
	my @ctl_file;
	my $valid = $fileio->read_file($ctl_file, \@ctl_file);
	unless ($valid) {  # Exit if we can't open the file
		die "\n\t ### Couldn't open control file '$ctl_file'\n\n\n ";
	}

	# Parse the control file
	$digs_obj->{ctl_file} = $ctl_file;
	my $loader_obj = ScreenBuilder->new($digs_obj);
	$loader_obj->parse_control_file($ctl_file, $digs_obj, $option);
	
	# Create the output directories if running a screen or re-assigning results table
	if ($option eq 2 or $option eq 3) { # Need output directory for these options

		$self->create_output_directories($digs_obj);

		# Store the ScreenBuilder object (TODO check this is used later)
		$loader_obj->{report_dir} = $digs_obj->{report_dir};
		$loader_obj->{tmp_path}   = $digs_obj->{tmp_path};
		$digs_obj->{loader_obj}   = $loader_obj;
	}

	my $db_name = $loader_obj->{db_name};
	unless ($db_name) { die "\n\t Error: no DB name defined \n\n\n"; }
	return $db_name;
}

#***************************************************************************
# Subroutine:  create output directories
# Description: create a unique 'report' directory for this process
#***************************************************************************
sub create_output_directories {
	
	my ($self, $digs_obj) = @_;

	# Create a unique ID and report directory for this run
	my $process_id  = $digs_obj->{process_id};
	my $output_path = $digs_obj->{output_path};
	unless ($process_id)  { die; }
	unless ($output_path) { die; }
	
	# $devtools->print_hash($digs_obj); die;
	my $report_dir  = $output_path . 'result_set_' . $process_id;
	$fileio->create_directory($report_dir);
	$digs_obj->{report_dir}  = $report_dir . '/';
	print "\n\t Created report directory";
	print "\n\t Path: '$report_dir'";
	
	# Create the tmp directory inside the report directory
	my $tmp_path = $report_dir . '/tmp';
	$fileio->create_directory($tmp_path);
	$digs_obj->{tmp_path}   = $tmp_path;

	# Create log file
	my $log_file   = $report_dir . "/log.txt";
	$fileio->append_text_to_file($log_file, "DIGS process $process_id\n");
	$digs_obj->{log_file} = $log_file;

}

#***************************************************************************
# Subroutine:  initialise_screening_db
# Description: load a DIGS screening database (create if doesn't exist) 
#***************************************************************************
sub initialise_screening_db {

	my ($self, $digs_obj, $db_name) = @_;

	# Create the screening DB object
	my $db_obj = ScreeningDB->new($digs_obj);

	# Check if this screening DB exists, if not then create it
	my $db_exists = $db_obj->does_db_exist($db_name);
	unless ($db_exists) {
		$db_obj->create_screening_db($db_name);	
	}
	
	# Load the table handles into screening database object
	print   "\n\n\t  Connecting to DB:  $db_name";
	$db_obj->load_screening_db($db_name);	
	$digs_obj->{db} = $db_obj; # Store the database object reference 
}

#***************************************************************************
# Subroutine:  setup_for_a_digs_run
# Description: prepare database and DIGS query list prior to screening
#***************************************************************************
sub setup_for_a_digs_run {

	my ($self, $digs_obj) = @_;

	# Flush active set
	my $db         = $digs_obj->{db};
	my $loader_obj = $digs_obj->{loader_obj};
	unless ($loader_obj) { die; }  # Sanity checking
	unless ($db)         { die "\n\t Error: no DB defined \n\n\n"; }

	#print "\n\t  Flushing 'active_set' table\n";
	my $active_set_table = $db->{active_set_table};
	$active_set_table->flush();
	
	# Index previously executed searches
	my %done;
	$self->index_previously_executed_searches($digs_obj, \%done);
	$loader_obj->{previously_executed_searches} = \%done;

	# Finally, set up the screen
	my %queries;
	my $total_queries = $loader_obj->setup_screen($digs_obj, \%queries);
	unless ($total_queries)  { 
		print "\n\t  Exiting DIGS setup";	
		return 0
	}
		
	# Record queries 
	$digs_obj->{queries}         = \%queries;
	$digs_obj->{total_queries}   = $total_queries;
	$digs_obj->{defragment_mode} = 'defragment';
	
	return 1;
}

#***************************************************************************
# Subroutine:  index_previously_executed_searches 
# Description: index BLAST searches that have previously been executed
#***************************************************************************
sub index_previously_executed_searches {
	
	my ($self, $digs_obj, $done_ref) = @_;

	my $db = $digs_obj->{db};	
	my $searches_table = $db->{searches_table};
	unless ($searches_table) { die "\n\t Searches_performed table not loaded\n\n"; }
	my @data;
	my @fields = qw [ record_id
                      probe_name probe_gene 
	                  organism target_datatype target_version target_name ];
	my $where = " ORDER BY record_id ";
	$searches_table->select_rows(\@fields, \@data, $where);
	
	# Index the executed searches
	foreach my $data_ref (@data) {
		
		# Get the query parameters
		my $organism        = $data_ref->{organism};
		my $target_datatype = $data_ref->{target_datatype};
		my $version         = $data_ref->{target_version};
		my $target_name     = $data_ref->{target_name};
		my $probe_name      = $data_ref->{probe_name};
		my $probe_gene      = $data_ref->{probe_gene};
	
		# Sanity checking
		unless ( $organism )        { die; }
        unless ( $target_datatype ) { die; }
        unless ( $version )         { die; }
        unless ( $target_name )     { die; }
        unless ( $probe_name )      { die; }
        unless ( $probe_gene )      { die; }
		
		# Create the unique key for this search
		my @genome = ( $organism , $target_datatype, $version );
		my $target_id = join ('|', @genome);
		my $probe_id  = $probe_name . '_' .  $probe_gene;
		my @key = ( $target_id, $target_name, $probe_id );
		my $key = join ('|', @key);

		# Record the query, indexed by it's unique key
		$done_ref->{$key} = $data_ref;		
	}

	#$devtools->print_hash($done_ref); die; # DEBUG
}

#***************************************************************************
# Subroutine:  setup_for_reassign
# Description: do general set up for a reassign process
#***************************************************************************
sub setup_for_reassign {

	my ($self, $digs_obj) = @_;

	my $loader_obj = $digs_obj->{loader_obj};
	my $where = '';
	unless ($digs_obj->{force}) {
		# Option to enter a WHERE statement
		my $question = "\n\n\t  Enter a WHERE statement to limit reaasign (Optional)";
		$where = $console->ask_question($question);
	}
	# Get the assigned digs_results

	# Get database tables
	my @reassign_loci;
	my $db = $digs_obj->{db};
	my $digs_results_table  = $db->{digs_results_table};
		
	# Set the fields to get values for
	my @fields = qw [ record_id assigned_name assigned_gene 
	                  probe_type sequence ];
	$digs_results_table->select_rows(\@fields, \@reassign_loci, $where);
	$digs_obj->{reassign_loci} = \@reassign_loci;
	
	# Set up the reference library
	$loader_obj->setup_reference_libraries($digs_obj);

}

#***************************************************************************
# Subroutine:  setup_for_defrag_or_consolidate
# Description: do general set up for a defragment or consolidate process
#***************************************************************************
sub setup_for_defrag_or_consolidate {

	my ($self, $digs_obj, $option) = @_;

	my $loader_obj = $self->{loader_obj};

	if ($option eq 4 or $option eq 5) { 

		# Set target sequence files for screening
		my %targets;
		my $num_targets = $loader_obj->set_targets(\%targets);

		# Show error and exit if no targets found
		unless ($num_targets) {
			$loader_obj->show_no_targets_found_error();
		}

		# Get the 'group' field for each target file (part of the path to the file)
		my %target_groups;
		$loader_obj->set_target_groups(\%targets, \%target_groups);
		$digs_obj->{target_groups} = \%target_groups; 
	}

	# DO SET-UP NEEDED FOR DEFRAGMENT ONLY
	if ($option eq 4) { 
		$digs_obj->{defragment_mode} = 'defragment';	
		# Set up the reference library
		$loader_obj->setup_reference_libraries($self);
	}
	# DO SET-UP NEEDED FOR CONSOLIDATE ONLY
	elsif ($option eq 5) { 
		$self->set_up_consolidate_tables();
		$digs_obj->{defragment_mode} = 'consolidate';

		# Get contig lengths and capture in a table
		$self->calculate_contig_lengths($digs_obj);	

		# Get the parameters for consolidation
		my $range = $self->{consolidate_range};
		my $d_range = $self->{defragment_range};
		unless ($range) { 
			my $question1 = "\n\n\t # Set the range for consolidating digs results";
			$range = $console->ask_int_with_bounds_question($question1, $d_range, $maximum);		
		}
	
		# Set the parameters for consolidation
		my %consolidate_settings;
		$consolidate_settings{range} = $range;
		$consolidate_settings{start} = 'extract_start';
		$consolidate_settings{end}   = 'extract_end';
		$digs_obj->{consolidate_settings} = \%consolidate_settings;
	}
}

#***************************************************************************
# Subroutine:  set_up_consolidate_tables
# Description:
#***************************************************************************
sub set_up_consolidate_tables {

	my ($self) = @_;

 	# Create tables if they don't exist already
	my $db_ref = $self->{db};
	my $dbh = $db_ref->{dbh};
	my $loci_exists = $db_ref->does_table_exist('loci');
	unless ($loci_exists) {
		$db_ref->create_loci_table($dbh);
	}
	my $loci_chains_exists = $db_ref->does_table_exist('loci_chains');
	unless ($loci_chains_exists) {
		$db_ref->create_loci_chains_table($dbh);
	}
	my $contigs_exists = $db_ref->does_table_exist('contigs');
	unless ($contigs_exists) {
		$db_ref->create_contigs_table($dbh);
	}

	# Load tables
	$db_ref->load_loci_table($dbh);
	$db_ref->load_loci_chains_table($dbh);
	$db_ref->load_contigs_table($dbh);

	# Get table references and set up for this consolidation process
	my $loci_table        = $db_ref->{loci_table};
	my $loci_chains_table = $db_ref->{loci_chains_table};
	my $contigs_table     = $db_ref->{contigs_table};
	$loci_table->flush();
	$loci_chains_table->flush();
	$contigs_table->flush();
}

#***************************************************************************
# Subroutine:  calculate_contig_lengths
# Description: 
#***************************************************************************
sub calculate_contig_lengths {

	my ($self, $digs_obj) = @_;

	# Get database and tables
	my $db_ref = $self->{db};
	my $digs_results     = $db_ref->{digs_results_table};
	my $contigs_table    = $db_ref->{contigs_table};
	my $genome_use_path  = $digs_obj->{genome_use_path};
	my $target_group_ref = $digs_obj->{target_groups};
	
	my $question1 = "\n\n\t # Refresh contig length table";
	my $refresh = $console->ask_yes_no_question($question1);		

	if ($refresh eq 'y') {
	
		my @targets;
		my @fields = qw [ organism target_datatype target_version target_name ];
		$digs_results->select_distinct(\@fields, \@targets);
		foreach my $target_ref(@targets) {
			
			# Read the file and get the lengths of each contig	
			# Get the target details (and thus the target path)	
			#my $target_path = $self->get_target_file_path($target_ref);
			my $organism        = $target_ref->{organism};
			my $target_name     = $target_ref->{target_name};
			my $target_datatype = $target_ref->{target_datatype};
			my $target_version  = $target_ref->{target_version};
			my @genome = ( $organism , $target_datatype, $target_version, $target_name );
			my $target_id       = join ('|', @genome);
			my $target_group = $target_group_ref->{$target_id};
			unless ($target_group) { die; 
				print " \n\t No target group found for TARGET ID $target_id\n\n"; 
        		sleep 1;
				next;
			}
		
			# Construct the path to this target file
			my @path;
			push (@path, $genome_use_path);
			push (@path, $target_group);
			push (@path, $organism);
			push (@path, $target_datatype);
			push (@path, $target_version);
			push (@path, $target_name);
			my $target_path = join ('/', @path);
			my @contigs;
			$fileio->read_fasta($target_path, \@contigs, 'true');
			
			foreach my $contig_ref (@contigs) {
			
				my $header = $contig_ref->{header};
				my $length = $contig_ref->{seq_length};
				unless ($length and $header) {
					$devtools->print_hash($contig_ref); die;
				}
				
				$contig_ref->{organism}        = $organism;
				$contig_ref->{target_datatype} = $target_datatype;
				$contig_ref->{target_version}  = $target_version;
				$contig_ref->{target_name}     = $target_name;
				$contig_ref->{scaffold}        = $header;
				print "\n\t # $header: $length";
				
				$contigs_table->insert_row($contig_ref);
			}			
		}
	}
}

############################################################################
# EOF
############################################################################