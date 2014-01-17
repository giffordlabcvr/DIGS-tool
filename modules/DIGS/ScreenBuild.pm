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
		
		# DB connection variables
		server                 => $parameter_ref->{server},
		username               => $parameter_ref->{username},
		password               => $parameter_ref->{password},
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# SCREENING SET UP
############################################################################

#***************************************************************************
# Subroutine:  set_up_screen
# Description: Set up all the queries to execute, indexed by genome chunk
#***************************************************************************
sub set_up_screen  {
	
	my ($self, $pipeline_obj, $queries_ref, $ctl_file) = @_;

	# Create the output directories
	$self->create_output_directories();
	
	# Set parameters for screening
	$self->parse_control_file($ctl_file);

	# Check we have everything we need
	$self->validate_screen_setup($self);
	
	# Load screening database (includes some MacroLineage Tables
	$self->set_screening_db();

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

	# transfer parameters from this object to the pipeline object
	$pipeline_obj->{tmp_path}           = $self->{tmp_path};
	$pipeline_obj->{blast_orf_lib_path} = $self->{blast_orf_lib_path};
	$pipeline_obj->{blast_utr_lib_path} = $self->{blast_utr_lib_path};
	$pipeline_obj->{seq_length_minimum} = $self->{seq_length_minimum};
	$pipeline_obj->{db}                 = $self->{db};

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
			die "\n\t Reference library protein FASTA not found\n\n\n";
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
			my $name     = $header_data{name};
			my $orf_name = $header_data{orf_name};
			my $aa_seq   = $seq_ref->{sequence};
			my $fasta = ">$name" . "_$orf_name" . "\n$aa_seq\n\n";
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
			my $orf_name = $header_data{orf_name};
			my $nt_seq   = $seq_ref->{sequence};
			my $fasta = ">$name" . "_$orf_name" . "\n$nt_seq\n\n";
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
	$fileio->write_output_file($aa_lib_path, $aa_lib_ref); 
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
       $fileio->write_output_file($nt_lib_path, $nt_lib_ref);

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
	my $target_paths_ref = $self->{target_paths};
	my @organisms = keys %$target_paths_ref;
	foreach my $organism (@organisms) {

		my $data_ref = $target_paths_ref->{$organism};
		my $path     = $data_ref->{path};
		unless ($path) { die; }
	
		#print "\n\t $organism: PATH $path";
		my $exists = $fileio->check_directory_exists($path);
		if ($exists) {
			
			my @leaves;
			$fileio->read_directory_tree_leaves_simple($path, \@leaves);
			
			#$devtools->print_array(\@leaves);
			foreach my $file_ref (@leaves) {
				my $file      = $file_ref->{file};
				#print "\n\t $file";
				my $file_type = $fileio->get_infile_type($file);
				if ($file_type eq 'fa') {
					$file_ref->{organism} = $organism;
					my $path = $file_ref->{path};
					my %path_elements;
					$self->get_path_elements(\%path_elements, $path); 

					my %data;
					$data{file}     = $file;
					$data{path}     = $path;
					$data{organism} = $organism;
					$data{version}  = $path_elements{version};
					my $key = $organism . '_' . $file;
					$targets_ref->{$key} = \%data;	
					#print "\n\t KEY $key";
				}
				elsif ($file_type eq 'fas') {
					die;
				}
			}
		}
	}
	#$devtools->print_hash($targets_ref); die; # DEBUG
}

#***************************************************************************
# Subroutine:  get_path_elements
# Description: get path elements
#***************************************************************************
sub get_path_elements {
	
	my ($self, $elements_ref, $path) = @_;

	#print "\n\t PATH '$path'";
	$path =~ s/\/\//\//g;
	my @path = split(/\//, $path);
	
	#$devtools->print_array(\@path); die;
	my $class    = $path[2];
	my $organism = $path[3];
	my $type     = $path[4];
	my $version  = $path[5];
	$elements_ref->{class}    = $class;
	$elements_ref->{organism} = $organism;
	$elements_ref->{type}     = $type;
	$elements_ref->{version}  = $version;
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
			my $orf_name = $header_data{orf_name};
			my $aa_seq   = $seq_ref->{sequence};
			#$devtools->print_hash(\%header_data);
			$self->add_aa_probe($probes_ref, $name, $orf_name, $aa_seq);
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
			my $orf_name = $header_data{orf_name};
			my $utr_seq   = $seq_ref->{sequence};
			$devtools->print_hash(\%header_data);
			$self->add_na_probe($probes_ref, $name, $orf_name, $utr_seq);
		}
	}
	die;
}

#***************************************************************************
# Subroutine:  parse_fasta_header_data
# Description:
#***************************************************************************
sub parse_fasta_header_data {
	
	my ($self, $header, $data_ref, $fail_count) = @_;

	my $name;
	my $orf_name;
	my @header = split (/_/, $header);
	$orf_name  = pop   @header;
	$name      = join('_', @header);

	unless ($name) { 
		$fail_count++;
		$name = "unknown_$fail_count";
	}
	unless ($orf_name) { 
		$orf_name = "unknown";
	}
	
	# DEBUG	
	#print "\n\t # HEADER $header";	
	#print "\n\t # NAME   $name";	
	#print "\n\t # ORF    $orf_name";	
	$data_ref->{name}     = $name;
	$data_ref->{orf_name} = $orf_name;

}

#***************************************************************************
# Subroutine:  add_aa_probe
# Description: 
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
# Description: 
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
		my $probe_id        = $probe_ref->{probe_id};
		my $sequence        = $probe_ref->{sequence};
		my $probe_len       = length $sequence;
		my $fasta = "\n>$probe_name\n$sequence";
		my $query_seq_file = $report_dir . $probe_id;
		$fileio->write_text_to_file($query_seq_file, $fasta);
		$probe_ref->{probe_path}   = $query_seq_file;
		$probe_ref->{probe_length} = $probe_len;
		$probe_ref->{result_path}  = $self->{tmp_path};
		print "\n\t probe $probe_id";		

		# Iterate through targets
		my @target_names = sort keys %$targets_ref;
		my $num_targets = scalar @target_names;
		unless ($num_targets) {
			die "\n\t no targets found\n\n\n";
		}
		foreach my $target_name (@target_names) {
			
			# Get target data
			my $target_ref   = $targets_ref->{$target_name};
			#$devtools->print_hash($target_ref); # die; # DEBUG
			my $organism     = $target_ref->{organism};
			my $target_path  = $target_ref->{path};
			my $target_name  = $target_ref->{file};		
			my $version      = $target_ref->{version};
			my @key = ( $organism, $version, $target_name, $probe_name, $probe_gene );
			my $key = join ('_', @key);
			if ($done{$key}) { 
				#print "\n\t ###### Skipping query: probe '$probe_id' vs '$target_name'";
				next; # Skip queries that have been issued
			} 
			
			# Else store the query
			print "\n\t ###### Setting query: probe '$probe_id' vs '$target_name'";
			$probe_ref->{organism}    = $organism;		
			$probe_ref->{target_name} = $target_name;		
			$probe_ref->{chunk_name}  = $target_name;	# TODO: remove (old field)
			$probe_ref->{target_path} = $target_ref->{path};;		
			$probe_ref->{version}     = $version;
			#$probe_ref->{source_type} = $target_ref->{source_type};

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
	
	# Read input file
	my @ctl_file;
	my $valid = $fileio->read_input_file($ctl_file, \@ctl_file);
	
	# Parse the PARAMS block
	my $start = 'BEGIN PARAMS';
	my $stop  = 'ENDBLOCK';
	$fileio->read_standard_field_value_block(\@ctl_file, $start, $stop, $self);
	#$devtools->print_hash($self);

	# Read the targets block
	#print "\n\t ### Reading target genomes block\n";
	my @target_block;
	$start = 'BEGIN TARGETS';
	$stop  = 'ENDBLOCK';
	$fileio->extract_text_block(\@ctl_file, \@target_block, $start, $stop);
	my @targets;
	foreach my $line (@target_block) {
		if ($line =~ /^\s*#/)   { next; } # discard comment line 
		chomp $line;
		push (@targets, $line);
	}

	# Get the target information
	#print "\n\t ### Getting target genome data\n";
	my %targets;
	my $genome_path = $self->{genome_use_path};
	foreach my $target_string (@targets) {
		my @split_string  = split("\/", $target_string);
		my $genome_group  = shift @split_string;
		my $organism      = shift @split_string;
		my $path		  = $genome_path . "/$genome_group/$organism/"; 
		my %data;
		$data{path}         = $path;
		$data{genome_group} = $genome_group;
		$data{organism}     = $organism;
		#print "\n\t Got genome in group $genome_group species $organism";
		$targets{$organism} = \%data;	
	}
	$self->{target_paths} = \%targets;	

}

#***************************************************************************
# Subroutine:  validate screen setup
# Description: read an input file to get parameters for screening
#***************************************************************************
sub validate_screen_setup {

	my ($self, $pipeline_ref, $queries_ref) = @_;

	# Get parameters inherited from Pipeline.pm
	my $process_id    = $self->{process_id};
	my $genome_path   = $self->{genome_use_path};
	my $output_path   = $self->{output_path};
	unless ($genome_path and $process_id and $output_path) { die; }

	# Get parameters parsed from ctl file
	my $db_name            = $self->{db_name};
	my $tblastn_min        = $self->{bit_score_min_tblastn};
	my $blastn_min         = $self->{bit_score_min_blastn};
	my $target_paths_ref   = $self->{target_paths}; # Get the target paths 
	my $query_aa_fasta     = $self->{query_aa_fasta};
	my $query_nt_fasta     = $self->{query_nt_fasta};
	my $query_glue         = $self->{query_glue};
	
	# Sanity checking - check that control file properly formed
	print "\n\n\t ### Validating screening parameters";
	unless ($db_name)          { die "\n\t No screening DB specified in ctl file"; }
	unless ($target_paths_ref) { die "\n\t No targets for screening defined\n\n";  }
	if ($query_aa_fasta) {
		unless ($blastn_min) { # Set to default minimum
			$self->{bit_score_min_tblastn} = $default_blastn_min;
		}
	}
	if ($query_nt_fasta) {
		unless ($tblastn_min) { # Set to default minimum
			$self->{bit_score_min_tblastn} = $default_tblastn_min;
		}
	}
	unless ($query_aa_fasta or $query_nt_fasta) {
		die "\n\t No probe set specified in control file\n\n";
	}
	#print "\n\t ORF lib, $blast_orf_lib_path, \n\t UTR lib, $blast_utr_lib_path";

	#else {
	#	die "\n\t No reference library for reciprocal BLAST defined\n\n"; 
	#}	
}

############################################################################
# EOF
############################################################################
