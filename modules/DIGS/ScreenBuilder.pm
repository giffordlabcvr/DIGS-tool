#!/usr/bin/perl -w
############################################################################
# Module:      ScreenBuilder.pm
# Description: Module for settings up a DIGS screen
# History:     December 2012: Created by Robert Gifford 
############################################################################
package ScreenBuilder;

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

# Program components
use DIGS::ScreeningDB;

############################################################################
# Globals
############################################################################

# Create base objects
my $fileio    = FileIO->new();
my $console   = Console->new();

# BLAST defaults
my $default_tblastn_min = 50;
my $default_blastn_min  = 50;
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Create a new ScreenBuilder.pm 'object'
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;
	my $self = {
		# Member variables
		process_id           => $parameter_ref->{process_id},
		# Flags
		refresh_genomes      => $parameter_ref->{refresh_genomes},
		# Paths
		genome_use_path      => $parameter_ref->{genome_use_path},
		blast_bin_path       => $parameter_ref->{blast_bin_path},
	};
	bless ($self, $class);
	return $self;
}

############################################################################
# TOP-LEVEL HANDLER
############################################################################

#***************************************************************************
# Subroutine:  set_up_screen
# Description: Set up all the queries to execute, indexed by genome target
#***************************************************************************
sub set_up_screen  {
	
	my ($self, $pipeline_obj, $queries_ref) = @_;

	# Create the output directories
	$self->create_output_directories($pipeline_obj);
	
	# Set up the reference library for BLAST
	$self->setup_reference_library($pipeline_obj);

	# Import the probes
	my @probes;
	$self->setup_blast_probes(\@probes);

	# Set target sequence files for screening
	my %targets;
	$self->set_targets(\%targets);
	
	# Create the BLAST queries for this screen
	my $num_queries = $self->set_queries($pipeline_obj, \@probes, \%targets, $queries_ref);
}

############################################################################
# MAIN FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  parse control file
# Description: read an input file to get parameters for screening
#***************************************************************************
sub parse_control_file {

	my ($self, $ctl_file, $pipeline_obj, $override) = @_;
	
	# Read input file
	my @ctl_file;
	my $valid = $fileio->read_file($ctl_file, \@ctl_file);

	# Parse the 'SCREENDB' block
	$self->parse_screendb_block(\@ctl_file);

	# Parse the 'SCREENSETS' block
	unless ($override) {
		$self->parse_screensets_block(\@ctl_file);
	}

	# READ the 'TARGETS' block
	unless ($override) {
		$self->parse_targets_block(\@ctl_file);
	}
	
	# Set parameters in pipeline object
	$pipeline_obj->{db_name}                = $self->{db_name};
	$pipeline_obj->{mysql_server}           = $self->{mysql_server};
	$pipeline_obj->{mysql_username}         = $self->{mysql_username};
	$pipeline_obj->{mysql_password}         = $self->{mysql_password};
	
	$pipeline_obj->{server}                 = $self->{mysql_server};
	$pipeline_obj->{password}               = $self->{mysql_password};
	$pipeline_obj->{username}               = $self->{mysql_username};
	
	$pipeline_obj->{output_path}            = $self->{output_path};
	$pipeline_obj->{blast_orf_lib_path}     = $self->{blast_orf_lib_path};
	$pipeline_obj->{blast_utr_lib_path}     = $self->{blast_utr_lib_path};
	
	$pipeline_obj->{seq_length_minimum}     = $self->{seq_length_minimum};
	$pipeline_obj->{bit_score_min_tblastn}  = $self->{bit_score_min_tblastn};
	$pipeline_obj->{bit_score_min_blastn}   = $self->{bit_score_min_blastn};
	$pipeline_obj->{redundancy_mode}        = $self->{redundancy_mode};
	$pipeline_obj->{threadhit_probe_buffer} = $self->{threadhit_probe_buffer};
	$pipeline_obj->{threadhit_gap_buffer}   = $self->{threadhit_gap_buffer};
	$pipeline_obj->{threadhit_max_gap}      = $self->{threadhit_max_gap};
	$pipeline_obj->{target_paths}           = $self->{target_paths};
}

#***************************************************************************
# Subroutine:  setup_reference_library
# Description: handler for reference library set up
#***************************************************************************
sub setup_reference_library {
	
	my ($self, $pipeline_obj) = @_;

	# Get reference library params
	my $reference_type;
	my $ref_fasta;
	if ($self->{reference_aa_fasta}) {
 		$self->{reference_library_type} = 'aa';
		$ref_fasta = $self->{reference_aa_fasta};
		$reference_type = 'amino acid';
	}
	elsif ($self->{reference_na_fasta}) {
 		$self->{reference_library_type} = 'na';
		$ref_fasta = $self->{reference_na_fasta};
		$reference_type = 'nucleic acid';
	}
	unless ($ref_fasta)  {die; }
	
	my @ref_fasta;
	my $num_fasta;
	if ($ref_fasta) {

		my @fasta;
		$self->read_fasta($ref_fasta, \@fasta);
		$num_fasta = scalar @fasta;
		unless ($num_fasta) { die "\n\t  Reference library: $reference_type FASTA not found'\n\n\n"; }

		print "\n\t  Reference library: $num_fasta $reference_type sequences";
		my $i = 0;
		foreach my $seq_ref (@fasta) {
			$i++;
			my $header  = $seq_ref->{header};
			$header  =~ s/\s+/_/g;
			my %header_data;
			my $valid = $self->parse_fasta_header_data($header, \%header_data);
			if ($valid) {
				my $name      = $header_data{name};
				my $gene_name = $header_data{gene_name};
				my $seq    = $seq_ref->{sequence};
				my $fasta = ">$name" . "_$gene_name" . "\n$seq\n\n";
				push (@ref_fasta, $fasta);
			}
		}
	}

	# Format the library for BLAST
	if ($num_fasta) {
		if ($self->{reference_library_type} eq 'aa') {
			$self->create_blast_lib(\@ref_fasta, 'aa');
		}
		elsif ($self->{reference_library_type} eq 'na') {
			$self->create_blast_lib(\@ref_fasta, 'na');
		}
		else { die; }
	}

	# Set the paths to the BLAST-formatted libraries
	$pipeline_obj->{blast_utr_lib_path} = $self->{blast_utr_lib_path};
	$pipeline_obj->{blast_orf_lib_path} = $self->{blast_orf_lib_path};
	
}

#***************************************************************************
# Subroutine:  setup_blast_probes
# Description: Setting up probes for a screen
#***************************************************************************
sub setup_blast_probes {
	
	my ($self, $probes_ref) = @_;

	# Get parameters from self
	my $probe_type;
	my $query_fasta;
	if ($self->{query_aa_fasta}) {
 		$self->{probe_library_type} = 'aa';
		$query_fasta = $self->{query_aa_fasta};
		$probe_type = 'amino acid';
	}
	elsif ($self->{query_na_fasta}) {
 		$self->{probe_library_type} = 'na';
		$query_fasta = $self->{query_na_fasta};
		$probe_type = 'nucleic acid';
	}
	else { 
		die "\n\t No path to probes setting has been loaded, check control file\n\n\n";
	}

	# Read FASTA probe library
	my @fasta;
	$self->read_fasta($query_fasta, \@fasta);
	my $num_fasta = scalar @fasta;
	print "\n\t  Probes:            $num_fasta $probe_type sequences";
	my $i = 0;
	my $type = $self->{probe_library_type};
	foreach my $seq_ref (@fasta) {
		$i++;
		my $header  = $seq_ref->{header};
		my %header_data;
		my $valid = $self->parse_fasta_header_data($header, \%header_data);
		if ($valid) {
			my $name      = $header_data{name};
			my $gene_name = $header_data{gene_name};
			my $seq       = $seq_ref->{sequence};
			
			my %probe;
			$probe{probe_name}      = $name;
			$probe{probe_gene}      = $gene_name;
			$probe{probe_id}        = $name . "_$gene_name";
			$probe{sequence}        = $seq;
			
			if ($type eq 'aa') {
				$probe{probe_type}  = 'ORF';
				$probe{blast_alg}   = 'tblastn';
				$probe{bitscore_cutoff} = $self->{bit_score_min_tblastn};
			}
			if ($type eq 'na') {
				$probe{probe_type}  = 'UTR';
				$probe{blast_alg}   = 'blastn';
				$probe{bitscore_cutoff} = $self->{bit_score_min_blastn};
			}
			push(@$probes_ref, \%probe);	

		}
		else {
			print "\n\t Couldn't extract data from header '$header'"; 
		}
	}
	unless ($i) {
		die "\n\t No Probes were loaded - check path is correct\n\n\n"; 
	}
}

#***************************************************************************
# Subroutine:  set_targets
# Description: get information about the target sequence data files
#***************************************************************************
sub set_targets {
	
	my ($self, $targets_ref) = @_;

	# Initialise target sequence library
	my $genome_obj = TargetDB->new($self);
	if ($self->{refresh_genomes}) {
		$genome_obj->refresh_genomes($targets_ref);
	} 
	
	# Iterate through targets set target file paths
	my %paths;
	my %target_data;
	my $genome_use_path  = $self->{genome_use_path};
	my $target_paths_ref = $self->{target_paths};
	unless ($target_paths_ref) { die; } 
	foreach my $path (@$target_paths_ref) {
		
		my $full_path = $genome_use_path . "/$path";	
		my $exists = $fileio->check_directory_exists($full_path);
		my @leaves;
		if ($exists) {
			$fileio->read_directory_tree_leaves_simple($full_path, \@leaves);
		}
		else {
			$path =~ s/\/\//\//g;
			my @path = split(/\//, $path);
			my $file = pop @path;
			my %file;
			$file{file} = $file;
			$file{path} = $full_path;
			push (@leaves, \%file);
		}
		$self->read_genome_files(\@leaves, $targets_ref);		
	}
	my @keys = keys %$targets_ref;
	my $unique_targets = scalar @keys;
	unless ($unique_targets) {
		print "\n\n\t No target databases found";
		print " - check target paths are correctly specified in control file\n";
		print "\n\t  \$DIGS_GENOMES path is set to '$ENV{DIGS_GENOMES}'\n";
		print "\n\t  TARGETS block from control file has these paths:";
		foreach my $path (@$target_paths_ref) {
			print "\n\t\t $path";
		}
		print "\n\n";
		exit;
	}
	print "\n\t  Targets:           $unique_targets target files";
}

#***************************************************************************
# Subroutine:  set_queries
# Description: set up the individual BLAST searches 
#***************************************************************************
sub set_queries {
	
	my ($self, $pipeline_obj, $probes_ref, $targets_ref, $queries_ref) = @_;

	my $db = $pipeline_obj->{db};
	unless ($db) { die; }
	
	# Get data from self
	my $report_dir = $self->{report_dir};
	my $tmp_path   = $self->{tmp_path};
	unless ($report_dir and $tmp_path) { die; } 

	# Get the target database information
	my @target_names = sort keys %$targets_ref;
	my $num_targets = scalar @target_names;
	unless ($num_targets) { die "\n\t No target databases found\n\n\n";	}

	# Work out the current state with respect to searches performed
	my %done;
	$db->index_previously_executed_queries(\%done);

	# Get relevant member variables and objects
	my $path;
	my $num_probes = scalar @$probes_ref;
	unless ($num_probes) { die "\n\t no probes found\n\n\n"; }

	my $outstanding;
	foreach my $probe_ref (@$probes_ref) {
		my $blast_alg       = $probe_ref->{blast_alg};
		my $bitscore_cutoff = $probe_ref->{bitscore_cutoff};
		my $probe_type      = $probe_ref->{probe_type};
		my $probe_name      = $probe_ref->{probe_name};
		my $probe_gene      = $probe_ref->{probe_gene};
		my $sequence        = $probe_ref->{sequence};
		my $probe_len       = length $sequence;
		my $fasta = "\n>$probe_name\n$sequence";
		my $probe_id  = $probe_name . '_' . $probe_gene;
		my $query_seq_file = $report_dir . $probe_id;
		$fileio->write_text_to_file($query_seq_file, $fasta);
		$probe_ref->{probe_path}   = $query_seq_file;
		$probe_ref->{probe_length} = $probe_len;
		$probe_ref->{result_path}  = $tmp_path;

		# Iterate through targets
		foreach my $target_name (@target_names) {
			
			# Get target data
			my $target_ref   = $targets_ref->{$target_name};
			my $organism     = $target_ref->{organism};
			my $data_type    = $target_ref->{data_type};
			my $version      = $target_ref->{version};
			my $target_path  = $target_ref->{path};
			my $target_name  = $target_ref->{file};		
		
			# Sanity checking	
			unless ( $organism and  $version and $data_type and
                     $target_name and $probe_name and $probe_gene ) {
			 		die;
			}
			
			# Create a unique key for this genome
			my @genome = ( $organism , $data_type, $version );
			my $genome_id = join ('|', @genome);

			# Create a unique key for this query
			my @key = ( $genome_id, $target_name, $probe_id );
			my $key = join ('|', @key);

			if ($done{$key}) { 
				next; # Skip queries that have been issued
			} 

			# Else store the query
			$outstanding++;
			$probe_ref->{genome_id}   = $genome_id;		
			$probe_ref->{organism}    = $organism;		
			$probe_ref->{version}     = $version;
			$probe_ref->{data_type}   = $data_type;		
			$probe_ref->{target_name} = $target_name;		
			$probe_ref->{target_path} = $target_ref->{path};;		

			# Important - create a copy
			my %query = %$probe_ref;
			if ($queries_ref->{$probe_name}) {
				my $probe_query_ref = $queries_ref->{$probe_name};
				push(@$probe_query_ref, \%query);
			}
			else {
				my @probe_queries;
				push(@probe_queries, \%query);
				$queries_ref->{$probe_name} = \@probe_queries;
			}
		}
	}

	# Show number of queries loaded
	unless ($outstanding) { print "\n\n\t ### No outstanding searches were loaded\n"; }
	else { print "\n\t  Searches to run    $outstanding\n"; }
	return $outstanding;
}

############################################################################
# INTERNAL FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  parse_screendb_block
# Description: read a NEXUS sttyle block from a DIGS control file
#***************************************************************************
sub parse_screendb_block {

	my ($self, $file_ref) = @_;
	
	# Parse the 'SCREENDB' block
	my $start = 'BEGIN SCREENDB';
	my $stop  = 'ENDBLOCK';
	my $db_block = $fileio->read_standard_field_value_block($file_ref, $start, $stop, $self);
	unless ($db_block)  {
		die "\n\t Control file error: no 'SCREENDB' block found\n\n\n";
	}
	
	# Get the 'SCREENDB' block values and validate
	my $db_name  = $self->{db_name};
	my $server   = $self->{mysql_server};
	my $user     = $self->{mysql_username};
	my $password = $self->{mysql_password};

	unless ($db_name)  {
		die "\n\t Control file error: 'db_name' undefined in 'SCREENDB' block\n\n\n";
	}
	unless ($server)  {
		die "\n\t Control file error: 'mysql_server' undefined in 'SCREENDB' block\n\n\n";
	}
	unless ($user)  {
		die "\n\t Control file error: 'mysql_username' undefined in 'SCREENDB' block\n\n\n";
	}
	unless ($password)  {
		die "\n\t Control file error: 'mysql_password' undefined in 'SCREENDB' block\n\n\n";
	}
}

#***************************************************************************
# Subroutine:  parse_screensets_block
# Description: read an input file to get parameters for screening
#***************************************************************************
sub parse_screensets_block {

	my ($self, $file_ref) = @_;
	
	# Parse the 'SCREENSETS' block
	my $start = 'BEGIN SCREENSETS';
	my $stop  = 'ENDBLOCK';
	my $block = $fileio->read_standard_field_value_block($file_ref, $start, $stop, $self);
	unless ($block)  {
		die "\n\n\t Control file error: no 'SCREENSETS' block found\n\n\n";
	}

	# Get the 'SCREENSETS' block values and validate
	my $tblastn_min            = $self->{bit_score_min_tblastn};
	my $blastn_min             = $self->{bit_score_min_blastn};
	my $query_aa_fasta         = $self->{query_aa_fasta};
	my $reference_aa_fasta     = $self->{reference_aa_fasta};
	my $query_na_fasta         = $self->{query_na_fasta};
	my $reference_na_fasta     = $self->{reference_na_fasta};
	my $query_glue             = $self->{query_glue};
	my $reference_glue         = $self->{reference_glue};
	my $redundancy_mode        = $self->{redundancy_mode};
	my $threadhit_probe_buffer = $self->{threadhit_probe_buffer};
	my $threadhit_gap_buffer   = $self->{threadhit_gap_buffer};
	my $threadhit_max_gap      = $self->{threadhit_max_gap};
	my $output_path            = $self->{output_path};

	unless ($output_path) {
		print "\n\t Warning no output path defined, results folder will be created in current directory\n\n\n";
	}
	# Validation for a nucleic acid, FASTA-based screen
	if ($query_aa_fasta) { # If a set of protein probes has been specified
		# Check if BLAST bitscore or evalue minimum set
		unless ($blastn_min) { # Set to default minimum
			$self->{bit_score_min_tblastn} = $default_blastn_min;
		}
		unless ($reference_aa_fasta) { # Set to default minimum
		  die "\n\t Control file error: no AA reference library defined for AA query set\n\n\n";
		}
	}
	# Validation for a nucleic acid, FASTA-based screen
	if ($query_na_fasta) {
		unless ($tblastn_min) { # Set to default minimum
			$self->{bit_score_min_tblastn} = $default_tblastn_min;
		}
		unless ($reference_na_fasta) { # Set to default minimum
		  die "\n\t Control file error: no NT reference library defined for NT query set\n\n\n";
		}
	}
	unless ($query_aa_fasta or $query_na_fasta or $query_glue) {
		die "\n\t Control file error: no probe library defined\n\n\n";
	}
	unless ($redundancy_mode) {
		# Set extract mode to default (extract everything)
		$self->{redundancy_mode} = 1;
	}
	unless ($threadhit_probe_buffer) {
		die "\n\t Control file error: 'Screensets' block parameter 'threadhit_probe_buffer' is undefined. \n\n\n";
	}
	unless ($threadhit_gap_buffer) {
		die "\n\t Control file error: 'Screensets' block parameter 'threadhit_probe_buffer' is undefined. \n\n\n";
	}
	unless ($threadhit_max_gap) {
		die "\n\t Control file error: 'Screensets' block parameter 'threadhit_max_gap' is undefined. \n\n\n";
	}
}

#***************************************************************************
# Subroutine:  parse_targets_block
# Description: get paths to the target sequence databases for DIGS
#***************************************************************************
sub parse_targets_block {

	my ($self, $file_ref) = @_;

	# READ the 'TARGETS' block
	my @target_block;
	my $start = 'BEGIN TARGETS';
	my $stop  = 'ENDBLOCK';
	$fileio->extract_text_block($file_ref, \@target_block, $start, $stop);
	my $screenset_lines = scalar @target_block;
	unless ($screenset_lines)  {
		die "\n\n\t Control file error: nothing in 'TARGETS' block\n\n\n";
	}
	
	# Parse the target strings
	my $targets = 0;
	my @targets;
	foreach my $line (@target_block) {
		chomp $line;
		$line =~ s/\s+//g; # remove whitespace
		if ($line =~ /^\s*$/)   { next; } # discard blank line
		if ($line =~ /^\s*#/)   { next; } # discard comment line 
		push (@targets, $line);
		$targets++;
	}
	$self->{target_paths} = \@targets;
}

#***************************************************************************
# Subroutine:  create output directories
# Description: create a unique 'report' directory for this process
#***************************************************************************
sub create_output_directories {
	
	my ($self, $pipeline_obj) = @_;

	# Create a unique ID and report directory for this run
	my $process_id   = $self->{process_id};
	my $output_path = $self->{output_path};
	my $report_dir  = $output_path . 'result_set_' . $process_id;
	$fileio->create_unique_directory($report_dir);
	$self->{report_dir}  = $report_dir . '/';
	
	# Create print "\n\t Report dir $report_dir"; die;
	my $tmp_path = $report_dir . '/tmp';
	$fileio->create_unique_directory($tmp_path);
	$self->{tmp_path}   = $tmp_path . '/';
	$pipeline_obj->{tmp_path}   = $tmp_path;
	$pipeline_obj->{report_dir} = $report_dir;
}

#***************************************************************************
# Subroutine:  create_blast_lib
# Description: create protein sequence library for reciprocal BLAST
#***************************************************************************
sub create_blast_lib {
	
	my ($self, $lib_ref, $type) = @_;

	# Get params from self
	my $report_dir   = $self->{report_dir};
	unless ($report_dir) { die; }	
	
	# Copy file to the report directory
	my $lib_path = $report_dir . "/reference_lib.fas";
	$fileio->write_file($lib_path, $lib_ref); 
	
	# Set path to blast binary
	my $blast_program = 'makeblastdb';
	my $blast_bin_dir = $self->{blast_bin_path};
	my $bin_path;
	if ($blast_bin_dir) {
		 $bin_path = $self->{blast_bin_path} . $blast_program;
	}
	else { $bin_path = $blast_program; }

	# Execute command
	my $makedb_cmd;
	if ($type eq 'aa') {
		$makedb_cmd = "$bin_path -in $lib_path -dbtype prot > /dev/null";
		$self->{blast_orf_lib_path} = $lib_path; 
	}
	elsif ($type eq 'na') {
		$makedb_cmd = "$bin_path -in $lib_path -dbtype nucl> /dev/null";
		$self->{blast_utr_lib_path} = $lib_path;
	}
	my $result = system $makedb_cmd;
	if ($result) {
		print "\n\t Failed to format reference library for BLAST! \n\n";
		print "\n\t $makedb_cmd \n\n"; exit;	
	}
}

#***************************************************************************
# Subroutine:  parse_fasta_header_data
# Description: parse elements out of a structured FASTA header
#              (Header is split into two elements using underscore) 
#***************************************************************************
sub parse_fasta_header_data {
	
	my ($self, $header, $data_ref) = @_;

	my $name;
	my $gene_name;
	my $valid = 1;

	# Remove illegal characters from the header line: these include:
	# / : * ? " < > |   because we need to write files using header elements
	# '                 because quotes interfere with SQL statements
	$header =~ s/\|//g;
	$header =~ s/\///g;
	$header =~ s/\*//g;
	$header =~ s/\?//g;
	$header =~ s/://g;
	$header =~ s/"//g;
	$header =~ s/<//g;
	$header =~ s/>//g;
	$header =~ s/'//g;
	$header =~ s/\s+//g;

	# Retrieve data from the header line
	my @header = split (/_/, $header);
	$gene_name  = pop   @header;
	$name      = join('_', @header);
	
	unless ($name and $gene_name) { 
		print "\n\t FASTA HEADER FORMAT ERROR";
		print "\n\t HEADER = '$header'";
		print "\n\t Headers should include two elements separated by an underscore\n\n";
		die;
	}
	$data_ref->{name}     = $name;
	$data_ref->{gene_name} = $gene_name;
}

#***************************************************************************
# Subroutine:  read genome files
# Description: processes the top level (leaves) of the genome directory
#***************************************************************************
sub read_genome_files {
	
	my ($self, $leaves_ref, $targets_ref) = @_;

	foreach my $file_ref (@$leaves_ref) {

		my $file = $file_ref->{file};
		my $path = $file_ref->{path};
	
		my $file_type = $fileio->get_infile_type($file);
		if ($file_type eq 'fa' or $file_type eq 'fas' or $file_type eq 'fasta') {
			
			$path =~ s/\/\//\//g;
			my @path = split(/\//, $path);
			my $file     = pop @path;
			my $version  = pop @path;
			my $type     = pop @path;
			my $organism = pop @path;
			my $group    = pop @path;
			unless ($organism and $type and $version) { die; }
			my @target = ( $organism , $type, $version, $file );
			my $target_id = join ('|', @target);

			# Store using key
			my %data;
			$data{file}      = $file;
			$data{path}      = $path;
			$data{organism}  = $organism;
			$data{version}   = $version;
			$data{data_type} = $type;
			$data{group}     = $group;
			$targets_ref->{$target_id} = \%data;	
		}
	}
}

#***************************************************************************
# Subroutine:  read_fasta
# Description: read a fasta file into an array of hashes. 
# Arguments:   $file: the name of the file to read
#              $array_ref: reference to the hash array to copy to
#***************************************************************************
sub read_fasta {

	my ($self, $file, $array_ref, $identifier) = @_;
	
	unless ($identifier) { $identifier = 'SEQ_'; }

	# Read in the file or else return
	unless (open(INFILE, $file)) {
		print "\n\t  Cannot open file \"$file\"\n\n";
		return undef;
	}

	# Use process ID and time to create unique ID stem
	my $pid  = $$;
	my $time = time;
	my $alias_stem = $identifier;
	
	# Iterate through lines in the file
	my @raw_fasta = <INFILE>;
	close INFILE;
	my $header;
    my $sequence;
   	my $i = 0;
	foreach my $line (@raw_fasta) {
		
		#print "\n## $i";
		chomp $line;
		if    ($line =~ /^\s*$/)   { next; } # discard blank line
		elsif ($line =~ /^\s*#/)   { next; } # discard comment line 
		elsif ($line =~ /^BEGIN/)  { last; } # stop if we reach a data block
		elsif ($line =~ /^>/) {
			
			$line =~ s/^>//g;
			# new header, store any sequence held in the buffer
			if ($header and $sequence) {
				$i++;
				my $alias_id = $alias_stem . "_$i";
				$sequence = uc $sequence;
				my %seq_obj;
				$seq_obj{sequence}    = $sequence;
				$seq_obj{header}      = $header;
				$seq_obj{sequence_id} = $alias_id;
				push(@$array_ref, \%seq_obj);
			}
			# reset the variables 
			$line =~ s/^>//;
			$header = $line;
			$sequence = undef;
		}
		else {
			# keep line, add to sequence string
            $sequence .= $line;
     	}
    }
	
	# Before exit, store any sequence held in the buffer
	if ($header and $sequence) {
		$i++;
		my $alias_id = $alias_stem . "_$i";
		$sequence =~ s/\s+//g; # Remove whitespace
		$sequence = uc $sequence;
		my %seq_obj;
		$seq_obj{sequence}    = $sequence;
		$seq_obj{header}      = $header;
		$seq_obj{sequence_id} = $alias_id;
		push(@$array_ref, \%seq_obj);
	}
}

############################################################################
# EOF
############################################################################
