#!/usr/bin/perl -w
############################################################################
# Module:      ScreenBuild.pm
# Description: Set up a screening process 
# History:     December 2011: Created by Robert Gifford 
############################################################################
package ScreenBuild;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::FileIO;
use Base::SeqIO;
use Base::DevTools;
use Base::Console;
use Base::Sequence;

# Program components
use DIGS::DB;

############################################################################
# Globals
############################################################################

# Create base objects
my $fileio    = FileIO->new();
my $seqio     = SeqIO->new();
my $devtools  = DevTools->new();
my $console   = Console->new();

# BLASTn min
my $default_tblastn_min = 50;
my $default_blastn_min  = 50;
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Parameters
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {
		
		# Member variables
		process_id           => $parameter_ref->{process_id},
		output_path          => $parameter_ref->{output_path},
		
		# Paths
		genome_use_path        => $parameter_ref->{genome_use_path},
		blast_bin_path         => $parameter_ref->{blast_bin_path},
		
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# SCREENING SET UP
############################################################################

#***************************************************************************
# Subroutine:  set_up_screen
# Description: Set up all the queries to execute, indexed by genome target
#***************************************************************************
sub set_up_screen  {
	
	my ($self, $pipeline_obj, $queries_ref, $ctl_file) = @_;

	# Create the output directories
	$self->create_output_directories();
	
	# Set parameters for screening
	$self->parse_control_file($ctl_file);

	# Load screening database (includes some MacroLineage Tables
	$self->set_screening_db();

	# Show metrics for db
	my $db_obj = $self->{db}; 
	unless ($db_obj) { die; }
	$db_obj->summarise_db();
	
	# Set up the reference library
	if ($self->{reference_aa_fasta}) {
		$self->load_aa_fasta_reference_library();
	}
	if ($self->{reference_nt_fasta}) {
		$self->load_nt_fasta_reference_library();
	}

	# Set up the target sequences to screen
	print "\n\n\t ### Getting the target sequences to screen";
	my %targets;
	$self->set_targets(\%targets);
		
	# Set up the probes
	print "\n\n\t ### Setting up the sequence 'probes' for screening";
	my @probes;
	if ($self->{query_aa_fasta}) {
		$self->load_aa_fasta_probes(\@probes);
	}
	if ($self->{query_nt_fasta}) {
		$self->load_nt_fasta_probes(\@probes);
	}

	# Create the list of BLAST queries for screening
	print "\n\n\t ### Creating the BLAST queries\n";
	$self->set_queries(\@probes, \%targets, $queries_ref);
	my $num_queries = scalar keys %$queries_ref;
	unless ( $num_queries ) {
		print "\n\t ### No screening queries were loaded\n\n\n";
		return 0;
	}
	#$devtools->print_hash($queries_ref); die;	# DEBUG

	# transfer parameters from this object to the pipeline object
	$pipeline_obj->{db}                 = $self->{db};
	$pipeline_obj->{mysql_server}       = $self->{mysql_server};
	$pipeline_obj->{mysql_username}     = $self->{mysql_username};
	$pipeline_obj->{mysql_password}     = $self->{mysql_password};
	$pipeline_obj->{tmp_path}           = $self->{tmp_path};
	$pipeline_obj->{blast_orf_lib_path} = $self->{blast_orf_lib_path};
	$pipeline_obj->{blast_utr_lib_path} = $self->{blast_utr_lib_path};
	$pipeline_obj->{seq_length_minimum} = $self->{seq_length_minimum};
	#$pipeline_obj->{select_list}        = $self->{select_list};
	#$pipeline_obj->{where_statement}    = $self->{where_statement};
	return 1;
}

############################################################################
# SET REFERENCE LIBRARY FOR RECIPROCAL BLAST
############################################################################

#***************************************************************************
# Subroutine:  set_screening_db
# Description: 
#***************************************************************************
sub set_screening_db {

	my ($self) = @_;

	my $db_name = $self->{db_name};
	my $db_obj = DB->new($self);
	$db_obj->load_screening_db($db_name);	
	$self->{db} = $db_obj; # Store the database object reference 
}
	
############################################################################
# SET REFERENCE LIBRARY FOR RECIPROCAL BLAST
############################################################################

#***************************************************************************
# Subroutine:  load_aa_fasta_reference_library
# Description: 
#***************************************************************************
sub load_aa_fasta_reference_library {
	
	my ($self) = @_;

	# Get reference library params
	my $ref_aa_fasta = $self->{reference_aa_fasta};
	unless ($ref_aa_fasta)  {die; }
	
	# Format a reference AA library	
	print "\n\n\t ### Loading FASTA AA reference sequences";
	my @ref_aa_fasta;
	my $num_fasta;
	if ($ref_aa_fasta) {

		my @fasta;
		$seqio->read_fasta($ref_aa_fasta, \@fasta);
		#unless ($status) { die "\n\t Input error: couldn't open FASTA probe library\n\n"; }
		$num_fasta = scalar @fasta;
		unless ($num_fasta) {
			die "\n\t Reference library protein FASTA not found'\n\n\n";
		}
		print "\n\n\t   '$num_fasta' FASTA formatted sequences in reference library";
		my $i = 0;
		my $fail_count = 0;
		foreach my $seq_ref (@fasta) {
			$i++;
			my $header  = $seq_ref->{header};
			my %header_data;
			$self->parse_fasta_header_data($header, \%header_data, $fail_count);
			#$devtools->print_hash(\%header_data);
			my $name      = $header_data{name};
			my $gene_name = $header_data{gene_name};
			my $aa_seq    = $seq_ref->{sequence};
			my $fasta = ">$name" . "_$gene_name" . "\n$aa_seq\n\n";
			push (@ref_aa_fasta, $fasta);
		}
	}

	# Create the libraries
	if ($num_fasta) {
		$self->create_blast_aa_lib(\@ref_aa_fasta);
	}
}

#***************************************************************************
# Subroutine:  load_nt_fasta_reference_library
# Description: 
#***************************************************************************
sub load_nt_fasta_reference_library {
	
	my ($self) = @_;

	# Get reference library params
	my $ref_nt_fasta = $self->{reference_nt_fasta};
	
	# Format a reference NT library	
	print "\n\t ### Loading FASTA reference sequences\n\n";
	my $num_fasta;
	my @ref_nt_fasta;
	if ($ref_nt_fasta) {
		my @fasta;
		$seqio->read_fasta($ref_nt_fasta, \@fasta);
		#unless ($status) { die "\n\t Input error: couldn't open FASTA probe library\n\n"; }
		$num_fasta = scalar @fasta;
		unless ($num_fasta) {
			die "\n\t Reference library NT FASTA not found\n\n\n";
		}
		print "\n\n\t   '$num_fasta' FASTA formatted sequences to be used as probes";
		my $i = 0;
		my $fail_count = 0;
		foreach my $seq_ref (@fasta) {
			$i++;
			my $header  = $seq_ref->{header};
			my %header_data;
			my $mode = $self->{ref_fasta_header_mode};
			$self->parse_fasta_header_data($header, \%header_data, $fail_count, $mode);
			my $name     = $header_data{name};
			my $gene_name = $header_data{gene_name};
			my $nt_seq   = $seq_ref->{sequence};
			my $fasta = ">$name" . "_$gene_name" . "\n$nt_seq\n\n";
			push (@ref_nt_fasta, $fasta);
		}
	}
	
	# Create the libraries
	if ($num_fasta) {
		$self->create_blast_nt_lib(\@ref_nt_fasta);
	}
}


#***************************************************************************
# Subroutine:  create_blast_aa_lib
# Description: create protein sequence library for reciprocal BLAST
#***************************************************************************
sub create_blast_aa_lib {
	
	my ($self, $aa_lib_ref) = @_;

	# Get params from self
	my $report_dir   = $self->{report_dir};
	unless ($report_dir) { die; }	
	
	my $blast_program = 'makeblastdb';
	my $aa_lib_path = $report_dir . "/reference_lib_aa.fas";
	$fileio->write_file($aa_lib_path, $aa_lib_ref); 
	my $blast_bin_dir = $self->{blast_bin_path};
	my $bin_path;
	if ($blast_bin_dir) {
		 $bin_path = $self->{blast_bin_path} . $blast_program;
	}
	else {
		 $bin_path = $blast_program;
	}
	my $makedb_cmd = "$bin_path -in $aa_lib_path > /dev/null";
	#print "\n\t $makedb_cmd \n\n"; die;	
	system $makedb_cmd;
	$self->{blast_orf_lib_path} = $aa_lib_path; 
}

#***************************************************************************
# Subroutine:  create_blast_nt_lib
# Description: create nucleotide sequence library for reciprocal BLAST
#***************************************************************************
sub create_blast_nt_lib {

       my ($self, $nt_lib_ref) = @_;

       # Get params from self
       my $report_dir   = $self->{report_dir};
       unless ($report_dir) { die; }
       my $blast_program = 'makeblastdb';
       my $nt_lib_path = $report_dir . "/reference_lib_nt.fas";
       $fileio->write_file($nt_lib_path, $nt_lib_ref);

       my $blast_bin_dir = $self->{blast_bin_path};
       my $bin_path;
       if ($blast_bin_dir) {
			$bin_path = $self->{blast_bin_path} . $blast_program;
       }
       else {
			 $bin_path = $blast_program;
       }

       my $makedb_cmd = "$bin_path -in $nt_lib_path -dbtype nucl> /dev/null";
       #print "\n\t $makedb_cmd \n\n";
       system $makedb_cmd;
       $self->{blast_utr_lib_path} = $nt_lib_path;
}

############################################################################
# SET TARGET SEQUENCE FILES (i.e. files of contigs)
############################################################################

#***************************************************************************
# Subroutine:  set_targets
# Description: get information about the target sequence data files
#***************************************************************************
sub set_targets {
	
	my ($self, $targets_ref) = @_;

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
		    #print "\n\t Reading leaves for PATH $full_path";
			$fileio->read_directory_tree_leaves_simple($full_path, \@leaves);
			#$devtools->print_array(\@leaves); die;
		}
		else {
			#die "\n\t # Couldn't open directory '$path'\n\n\n";
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
	#$devtools->print_hash($targets_ref); die; # DEBUG
}


#***************************************************************************
# Subroutine:  read genome files
# Description: processes the top level (leaves) of the genome directory
#***************************************************************************
sub read_genome_files {
	
	my ($self, $leaves_ref, $targets_ref) = @_;

	#print "\n\t PATH '$path'";
	foreach my $file_ref (@$leaves_ref) {

		my $file = $file_ref->{file};
		my $path = $file_ref->{path};
	
		my $file_type = $fileio->get_infile_type($file);
		#if ($file_type eq 'fa' or $file_type eq 'fas' or $file_type eq 'fasta') {
		if ($file_type eq 'fa') {
			
			my %path_elements;
			$self->get_path_elements(\%path_elements, $path); 
			
			my $organism     = $path_elements{organism};
			my $data_type    = $path_elements{data_type};
			my $version      = $path_elements{version};
			unless ($organism and $data_type and $version) { die; }
			my @target = ( $organism , $data_type, $version, $file );
			my $target_id = join ('|', @target);
			#print "\n\t KEY $target_id"; die;

			# Store using key
			my %data;
			$data{file}      = $file;
			$data{path}      = $path;
			$data{organism}  = $organism;
			$data{version}   = $version;
			$data{data_type} = $data_type;
			$data{group}     = $path_elements{group};
			$targets_ref->{$target_id} = \%data;	
		}
	}


}

#***************************************************************************
# Subroutine:  get_path_elements
# Description: get the directory names within a genome path string 
#***************************************************************************
sub get_path_elements {
	
	my ($self, $elements_ref, $path) = @_;

	#print "\n\t PATH '$path'";
	$path =~ s/\/\//\//g;
	my @path = split(/\//, $path);
	my $file     = pop @path;
	my $version  = pop @path;
	my $type     = pop @path;
	my $organism = pop @path;
	my $group    = pop @path;
	$elements_ref->{organism}  = $organism;
	$elements_ref->{version}   = $version;
	$elements_ref->{data_type} = $type;
	$elements_ref->{group}     = $group;
	$elements_ref->{file}      = $file;
	#$devtools->print_array(\@path); die;
	#$devtools->print_hash($elements_ref); die;

}

############################################################################
# SET UP PROBES
############################################################################

#***************************************************************************
# Subroutine:  load_aa_fasta_probes
# Description: load probes from FASTA library 
#***************************************************************************
sub load_aa_fasta_probes {
	
	my ($self, $probes_ref) = @_;

	# Get parameters from self
	my $query_aa_fasta = $self->{query_aa_fasta};

	# Read FASTA nt probe library
	if ($query_aa_fasta) {
		my @fasta;
		$seqio->read_fasta($query_aa_fasta, \@fasta);
		#unless ($status) { die "\n\t Input error: couldn't open FASTA probe library\n\n"; }
		my $num_fasta = scalar @fasta;
		print "\n\n\t   '$num_fasta' FASTA formatted sequences will be used as probes";
		my $i = 0;
		my $fail_count = 0;
		foreach my $seq_ref (@fasta) {
			$i++;
			my $header  = $seq_ref->{header};
			my %header_data;
			$self->parse_fasta_header_data($header, \%header_data, $fail_count);
			my $name     = $header_data{name};
			my $gene_name = $header_data{gene_name};
			my $aa_seq   = $seq_ref->{sequence};
			#$devtools->print_hash(\%header_data);
			$self->add_aa_probe($probes_ref, $name, $gene_name, $aa_seq);
		}
	}
}

#***************************************************************************
# Subroutine:  load_nt_fasta_probes
# Description: load probes from FASTA library 
#***************************************************************************
sub load_nt_fasta_probes {

	my ($self, $probes_ref) = @_;

	# Get parameters from self
	my $query_nt_fasta = $self->{query_nt_fasta};

	# Read FASTA nt probe library
	if ($query_nt_fasta) {
		my @fasta;
		$seqio->read_fasta($query_nt_fasta, \@fasta);
		#unless ($status) { die "\n\t Input error: couldn't open FASTA probe library\n\n"; }
		my $num_fasta = scalar @fasta;
		print "\n\t '$num_fasta' FASTA formatted reference sequences will be used as probes";
		my $i = 0;
		my $fail_count = 0;
		foreach my $seq_ref (@fasta) {
			$i++;
			my $header  = $seq_ref->{header};
			my %header_data;
			$self->parse_fasta_header_data($header, \%header_data, $fail_count);
			my $name     = $header_data{name};
			my $gene_name = $header_data{gene_name};
			my $utr_seq   = $seq_ref->{sequence};
			#$devtools->print_hash(\%header_data);
			$self->add_na_probe($probes_ref, $name, $gene_name, $utr_seq);
		}
	}
}

#***************************************************************************
# Subroutine:  parse_fasta_header_data
# Description:
#***************************************************************************
sub parse_fasta_header_data {
	
	my ($self, $header, $data_ref, $fail_count) = @_;

	my $name;
	my $gene_name;
	my @header = split (/_/, $header);
	$gene_name  = pop   @header;
	$name      = join('_', @header);

	unless ($name) { 
		$fail_count++;
		$name = "unknown_$fail_count";
	}
	unless ($gene_name) { 
		$gene_name = "unknown";
	}
	
	$data_ref->{name}     = $name;
	$data_ref->{gene_name} = $gene_name;

	# DEBUG	
	#print "\n\t # HEADER $header";	
	#print "\n\t # NAME   $name";	
	#print "\n\t # ORF    $gene_name";	
}

#***************************************************************************
# Subroutine:  add_aa_probe
# Description: add an amino acid probe sequence to a BLAST query definition 
#***************************************************************************
sub add_aa_probe {
	
	my ($self, $probes_ref, $refseq_name, $aa_seq_name, $seq) = @_;

	my $bitscore_min = $self->{bit_score_min_tblastn};
	my %probe;
	$probe{blast_alg}       = 'tblastn';
	$probe{bitscore_cutoff} = $bitscore_min;
	$probe{probe_type}      = 'ORF';
	$probe{probe_name}      = $refseq_name;
	$probe{probe_gene}      = $aa_seq_name;
	$probe{probe_id}        = $refseq_name . "_$aa_seq_name";
	$probe{sequence}        = $seq;
	push(@$probes_ref, \%probe);	
}

#***************************************************************************
# Subroutine:  add_na_probe
# Description: add a nucleic acid probe sequence to a BLAST query definition 
#***************************************************************************
sub add_na_probe {
	
	my ($self, $probes_ref, $refseq_name, $na_seq_name, $seq) = @_;

	my $bitscore_min = $self->{bit_score_min_blastn};
	my %probe;
	$probe{blast_alg}       = 'blastn';
	$probe{bitscore_cutoff} = $bitscore_min;
	$probe{probe_type}      = 'UTR';
	$probe{probe_name}      = $refseq_name;
	$probe{probe_gene}      = $na_seq_name;
	$probe{probe_id}        = $refseq_name . "_$na_seq_name";
	$probe{sequence}        = $seq;
	push(@$probes_ref, \%probe);	

}

############################################################################
# SET UP QUERIES FOR SCREENING
############################################################################

#***************************************************************************
# Subroutine:  set_queries
# Description: set up the individual BLAST searches 
#***************************************************************************
sub set_queries {
	
	my ($self, $probes_ref, $targets_ref, $queries_ref) = @_;

	# Get data from self
	my $report_dir = $self->{report_dir};
	#$devtools->print_hash($targets_ref); #die; # DEBUG

	# Work out the current state with respect to searches performed
	my %done;
	my $db_obj = $self->{db};
	$db_obj->index_previously_executed_queries(\%done);
	#$devtools->print_hash(\%done); die; # DEBUG

	# Get relevant member variables and objects
	my $path;
	my $num_probes = scalar @$probes_ref;
	unless ($num_probes) {
		die "\n\t no probes found\n\n\n";
	}
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
		$probe_ref->{result_path}  = $self->{tmp_path};

		# Iterate through targets
		my @target_names = sort keys %$targets_ref;
		my $num_targets = scalar @target_names;
		unless ($num_targets) {
			die "\n\t no targets found\n\n\n";
		}
		foreach my $target_name (@target_names) {
			
			# Get target data
			my $target_ref   = $targets_ref->{$target_name};
			my $organism     = $target_ref->{organism};
			my $data_type    = $target_ref->{data_type};
			my $version      = $target_ref->{version};
			my $target_path  = $target_ref->{path};
			my $target_name  = $target_ref->{file};		
			#$devtools->print_hash($target_ref); # die; # DEBUG
			unless ( $organism and  $version and $data_type and
                     $target_name and $probe_name and $probe_gene ) {
			 		die;
			}
			my @genome = ( $organism , $data_type, $version );
			my $genome_id = join ('|', @genome);
			my @key = ( $genome_id, $target_name, $probe_id );
			my $key = join ('|', @key);
			if ($done{$key}) { 
				#print "\n\t ###### Skipping query: probe '$probe_id' vs '$target_name'";
				next; # Skip queries that have been issued
			} 
			#$devtools->print_hash(\%done); 
			#print "\n\t ###### KEY '$key'"; #die;	

			# Else store the query
			print "\n\t ###### Setting query: probe '$probe_id' vs '$target_name'";
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
}

############################################################################
# PARSING ETC
############################################################################

#***************************************************************************
# Subroutine:  create output directories
# Description: create a unique 'report' directory for this process
#***************************************************************************
sub create_output_directories {
	
	my ($self) = @_;

	# Get process ID to create unique output directory
	my $process_id   = $self->{process_id};
	print "\n\n\t ### Screening process ID is '$process_id'";
	
	# Create a unique ID and report directory for this run
	my $output_path = $self->{output_path};
	my $report_dir  = $output_path . $process_id;
	$fileio->create_unique_directory($report_dir);
	$self->{report_dir}  = $report_dir . '/';
	
	# Create print "\n\t Report dir $report_dir"; die;
	my $tmp_path = $report_dir . '/tmp';
	$fileio->create_unique_directory($tmp_path);
	$self->{tmp_path}   = $tmp_path . '/';
	#print "\n\t tmp path '$tmp_path'"; die;

}

#***************************************************************************
# Subroutine:  parse control file
# Description: read an input file to get parameters for screening
#***************************************************************************
sub parse_control_file {

	my ($self, $ctl_file) = @_;
	
	# Get parameters inherited from Pipeline.pm
	my $process_id    = $self->{process_id};
	my $genome_path   = $self->{genome_use_path};
	my $output_path   = $self->{output_path};
	unless ($genome_path and $process_id and $output_path) { die; }
	
	# Read input file
	my @ctl_file;
	my $valid = $fileio->read_file($ctl_file, \@ctl_file);

	# Parse the 'SCREENDB' block
	my $start = 'BEGIN SCREENDB';
	my $stop  = 'ENDBLOCK';
	my $db_block = $fileio->read_standard_field_value_block(\@ctl_file, $start, $stop, $self);
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

	# Parse the 'SCREENSETS' block
	$start = 'BEGIN SCREENSETS';
	$stop  = 'ENDBLOCK';
	my $block = $fileio->read_standard_field_value_block(\@ctl_file, $start, $stop, $self);
	unless ($block)  {
		die "\n\n\t Control file error: no 'SCREENSETS' block found\n\n\n";
	}

	# Get the 'SCREENSETS' block values and validate
	my $tblastn_min        = $self->{bit_score_min_tblastn};
	my $blastn_min         = $self->{bit_score_min_blastn};
	my $query_aa_fasta     = $self->{query_aa_fasta};
	my $query_nt_fasta     = $self->{query_nt_fasta};
	my $reference_aa_fasta = $self->{reference_aa_fasta};
	my $reference_nt_fasta = $self->{reference_nt_fasta};

	# Check the probe files and correspondence to parameters for BLAST
	if ($query_aa_fasta) { # If a set of protein probes has been specified
		# Check if BLAST bitscore or evalue minimum set
		unless ($blastn_min) { # Set to default minimum
			$self->{bit_score_min_tblastn} = $default_blastn_min;
		}
		unless ($reference_aa_fasta) { # Set to default minimum
		  die "\n\t Control file error: no AA reference library defined for AA query set\n\n\n";
		}
		# TODO Attempt to read the sequences
		# Validate reference and probe FASTA

	}
	if ($query_nt_fasta) {
		unless ($tblastn_min) { # Set to default minimum
			$self->{bit_score_min_tblastn} = $default_tblastn_min;
		}
		unless ($reference_nt_fasta) { # Set to default minimum
		  die "\n\t Control file error: no NT reference library defined for NT query set\n\n\n";
		}
		# TODO Attempt to read the sequences
		# Validate reference and probe FASTA
	}
	unless ($query_aa_fasta or $query_nt_fasta) {
		die "\n\t Control file error: no probe library defined\n\n\n";
	}

	# READ the 'TARGETS' block
	my @target_block;
	$start = 'BEGIN TARGETS';
	$stop  = 'ENDBLOCK';
	$fileio->extract_text_block(\@ctl_file, \@target_block, $start, $stop);
	my $screenset_lines = scalar @target_block;
	unless ($screenset_lines)  {
		die "\n\n\t Control file error: nothing in 'TARGETS' block\n\n\n";
	}
	
	# Parse the target strings
	my $targets = 0;
	my @targets;
	foreach my $line (@target_block) {
		if ($line =~ /^\s*#/)   { next; } # discard comment line 
		chomp $line;
		push (@targets, $line);
		$targets++;
	}
	$self->{target_paths} = \@targets;

	# READ the 'SCREENSQL' block
	$start = 'BEGIN SCREENSQL';
	$stop  = 'ENDBLOCK';
	my $sql_block = $fileio->read_sql_block(\@ctl_file, $start, $stop, $self);
	unless ($sql_block)  {
		die "\n\n\t Control file error: nothing in 'SCREENSQL' block\n\n\n";
	}
	my $select_list     = $self->{select_list};
	my $where_statement = $self->{where_statement};
	#$devtools->print_hash($self); die;
}

############################################################################
# EOF
############################################################################
