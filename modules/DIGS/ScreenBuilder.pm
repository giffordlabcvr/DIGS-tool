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
my $devtools  = DevTools->new();
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
		verbose              => $parameter_ref->{verbose},

		# Paths
		genome_use_path      => $parameter_ref->{genome_use_path},
		blast_bin_path       => $parameter_ref->{blast_bin_path},
		
		# Member classes 
		blast_obj            => $parameter_ref->{blast_obj},
		
		# Data
		previously_executed_searches => $parameter_ref->{previously_executed_searches},

	};
	bless ($self, $class);
	return $self;
}

############################################################################
# TOP-LEVEL HANDLER
############################################################################

#***************************************************************************
# Subroutine:  setup_screen
# Description: Set up all the queries to execute, indexed by genome target
#***************************************************************************
sub setup_screen  {
	
	my ($self, $pipeline_obj, $queries_ref) = @_;
	
	# Import the probes
	my @probes;
	$self->setup_blast_probes(\@probes);
	
	# Set up the reference library for BLAST
	$self->setup_reference_libraries($pipeline_obj);
	
	# Set target sequence files for screening
	my %targets;
	my %target_groups;
	$self->set_targets(\%targets, \%target_groups);
	$pipeline_obj->{target_groups} = \%target_groups; 

	# Create the BLAST queries for this screen
	my $num_queries = $self->set_queries($pipeline_obj, \@probes, \%targets, $queries_ref);
	my @keys        = keys %targets;
	my $unique      = scalar @keys;
	
	# Show output about the number of oustanding targets and number of queries loaded
	print "\n\t  Targets:           $unique target files";
	unless ($num_queries) { print "\n\n\t ### No outstanding searches were loaded\n"; }
	else { print "\n\t  Searches to run    $num_queries\n"; }
	return $num_queries;
}

############################################################################
# MAIN FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  parse control file
# Description: read an input file to get parameters for screening
#***************************************************************************
sub parse_control_file {

	my ($self, $ctl_file, $pipeline_obj, $option) = @_;
	
	# Read input file
	my @ctl_file;
	my $valid = $fileio->read_file($ctl_file, \@ctl_file);

	# Parse the 'SCREENDB' block (screening DB name and MySQL connection details)
	$self->parse_screendb_block(\@ctl_file);

	# Parse the 'SCREENSETS' block (input data and parameters for screening)
	$self->parse_screensets_block(\@ctl_file);

	# READ the 'TARGETS' block (files to be screened)	
	my @targets;
	my $start_token = 'BEGIN TARGETS';
	my $stop_token  = 'ENDBLOCK';
	$self->parse_target_block(\@ctl_file, $start_token, $stop_token, \@targets);
	$self->{target_paths} = \@targets;

	# READ the 'EXCLUDE' block	(files in target path to be excluded from screen)	
	my @exclude;
	$start_token = 'BEGIN EXCLUDE';
	$stop_token  = 'ENDBLOCK';
	$self->parse_target_block(\@ctl_file, $start_token, $stop_token, \@exclude);
	$self->{exclude_paths} = \@exclude;

	# READ the 'SKIPINDEX' block if present (Targets to skip during indexing)
	my @skipindex;
	$start_token = 'BEGIN SKIPINDEX';
	$stop_token  = 'ENDBLOCK';
	$self->parse_target_block(\@ctl_file, $start_token, $stop_token, \@skipindex);
	$self->{skipindexing_paths} = \@skipindex;

	# Set parameters in pipeline object
	
	# Screening DB name and MySQL connection details
	$pipeline_obj->{db_name}                = $self->{db_name};
	$pipeline_obj->{mysql_server}           = $self->{mysql_server};
	$pipeline_obj->{mysql_username}         = $self->{mysql_username};
	$pipeline_obj->{mysql_password}         = $self->{mysql_password};

	# Input and output file paths 	
	$pipeline_obj->{output_path}            = $self->{output_path};
	$pipeline_obj->{blast_orf_lib_path}     = $self->{blast_orf_lib_path};
	$pipeline_obj->{blast_utr_lib_path}     = $self->{blast_utr_lib_path};
	$pipeline_obj->{target_paths}           = $self->{target_paths};
	$pipeline_obj->{skipindexing_paths}     = $self->{skipindexing_paths};
	
	# Set parameters for screening
	$pipeline_obj->{defragment_range}       = $self->{defragment_range};
	$pipeline_obj->{consolidate_range}      = $self->{consolidate_range};
	$pipeline_obj->{extract_buffer}         = $self->{extract_buffer};

	# Set numthreads in the BLAST object 
	my $num_threads = $self->{num_threads};
	unless ($num_threads) { $num_threads = 1; }  # Default setting
	$pipeline_obj->{blast_obj}->{num_threads} = $num_threads;

	if ($option eq 6) {
	
		# READ the 'NOMENCLATURE' block
		$start_token  = 'BEGIN NOMENCLATURE';
		$stop_token   = 'ENDBLOCK';
		$self->parse_nomenclature_block(\@ctl_file, $start_token, $stop_token);
	
		# Set paths for applying nomenclature 	
		$pipeline_obj->{new_track_path}        = $self->{new_track_path};
		$pipeline_obj->{translation_path}      = $self->{translation_path};
		$pipeline_obj->{tax_level}             = $self->{tax_level};
		$pipeline_obj->{organism_code}         = $self->{organism_code};
		$pipeline_obj->{locus_class}           = $self->{locus_class};
		$pipeline_obj->{nomenclature_version}  = $self->{nomenclature_version};
		$pipeline_obj->{nomenclature_organism} = $self->{nomenclature_organism};
		$pipeline_obj->{genome_structure}      = $self->{genome_structure};
	}

	# Set the thresholds for filtering results
	$pipeline_obj->{seq_length_minimum}     = $self->{seq_length_minimum};

	# Set the bit score minimum for a tBLASTn screen
	if ($self->{bitscore_min_tblastn}) {
		$pipeline_obj->{bitscore_minimum} = $self->{bitscore_min_tblastn};
	}
	# Set the bit score minimum for a BLASTn screen
	if ($self->{bitscore_min_blastn}) {
		$pipeline_obj->{bitscore_minimum} = $self->{bitscore_min_blastn};	
	}

}

#***************************************************************************
# Subroutine:  setup_reference_libraries
# Description: handler for reference library set up
#***************************************************************************
sub setup_reference_libraries {

	my ($self, $pipeline_obj) = @_;

	# Get reference library params
	my $reference_type;
	my $ref_fasta;
	my @ref_fasta;
	my $num_fasta;
	
	# Load a peptide sequence library
	if ($self->{reference_aa_fasta}) {
 		$self->{reference_library_type} = 'aa';
		$ref_fasta = $self->{reference_aa_fasta};
		$reference_type = 'amino acid';
		$self->create_reference_library($pipeline_obj, $ref_fasta, $reference_type);
	}
	
	# Load a nucleotide sequence library	
	if ($self->{reference_na_fasta}) {
 		$self->{reference_library_type} = 'na';
		$ref_fasta = $self->{reference_na_fasta};
		$reference_type = 'nucleic acid';
		$self->create_reference_library($pipeline_obj, $ref_fasta, $reference_type);
	}
}

#***************************************************************************
# Subroutine:  create_reference_library
# Description: 
#***************************************************************************
sub create_reference_library {
	
	my ($self, $pipeline_obj, $ref_fasta, $reference_type) = @_;

	my @path = split(/\//, $ref_fasta);
	my $lib_file = pop @path;

	my $nonunique = 0;
	my %nonunique;
	my @fasta;
	$self->read_fasta($ref_fasta, \@fasta);
	my $num_fasta = scalar @fasta;
	unless ($num_fasta) { die "\n\t  Reference library: $reference_type FASTA not found'\n\n"; }

	print "\n\t  Reference library: $num_fasta $reference_type sequences";
	my $i = 0;
	#$devtools->print_array(\@fasta); die;
	my @references;
	my %refseq_ids; # Hash to check probe names are unique		
	foreach my $seq_ref (@fasta) {

		#$devtools->print_hash($seq_ref);
		$i++;
		my $header = $seq_ref->{header};

		$header    =~ s/\s+/_/g;
		my %header_data;
		my $valid = $self->parse_fasta_header_data($header, \%header_data);
		unless ($valid) { die; }
		my $name      = $header_data{name};
		my $gene_name = $header_data{gene_name};
		my $refseq_id = $name . "_$gene_name";
		my $seq    = $seq_ref->{sequence};

		if ($refseq_ids{$refseq_id}) {
			$nonunique++;
			if   ($nonunique{$refseq_id}) { $nonunique{$refseq_id}=1; }
			else                          { $nonunique{$refseq_id}++; } 
			#print "\n\t\t  Warning: non-unique reference name '$refseq_id'";
			#sleep 1;
		}
		else {
			$refseq_ids{$refseq_id} = 1;
		}
		my $fasta = ">$name" . "_$gene_name" . "\n$seq\n\n";
		push (@references, $fasta);
	}

	if ($nonunique) {
		print "\n\t   - Warning: $nonunique non-unique names identified in reference library";
		sleep 1;
		#print "\n\t\t  Warning: non-unique reference name '$refseq_id'";
	}

	# Format the library for BLAST
	if ($num_fasta) {
		if ($self->{reference_library_type} eq 'aa') {
			$self->create_blast_lib(\@references, 'aa');
			$pipeline_obj->{aa_reference_library} = $lib_file;
		}
		if ($self->{reference_library_type} eq 'na') {
			$self->create_blast_lib(\@references, 'na');
			$pipeline_obj->{na_reference_library} = $lib_file;
		}
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

	# Set up peptide probes
	my $probe_type;
	my $input_path;
	my $got_probes = undef;
	if ($self->{query_aa_fasta}) {
 		$self->{probe_library_type} = 'aa';
		$input_path = $self->{query_aa_fasta};
		$probe_type = 'amino acid FASTA';
		$self->get_fasta_probes($probes_ref, $input_path, $probe_type);
		$got_probes = scalar @$probes_ref;;
	}

	# Set up nucleotide probes
	if ($self->{query_na_fasta}) {
 		$self->{probe_library_type} = 'na';
		$input_path = $self->{query_na_fasta};
		$probe_type = 'nucleic acid FASTA';
		$self->get_fasta_probes($probes_ref, $input_path, $probe_type);
		$got_probes = scalar @$probes_ref;;
	}
	
	unless ($got_probes) { 
		$devtools->print_hash($self);
		die "\n\t No path to probes setting has been loaded, check control file\n\n\n";
	}
}

#***************************************************************************
# Subroutine:  get_fasta_probes
# Description: 
#***************************************************************************
sub get_fasta_probes {
	
	my ($self, $probes_ref, $query_fasta, $probe_type) = @_;

	# Read FASTA probe library
	my @fasta;
	$self->read_fasta($query_fasta, \@fasta);
	my $num_fasta = scalar @fasta;
	print "\n\t  Probe sequences:   $num_fasta $probe_type sequences";
	my $i = 0;
	my $type = $self->{probe_library_type};
	my %probe_ids; # Hash to check probe names are unique
	foreach my $seq_ref (@fasta) {
		$i++;
		my $header  = $seq_ref->{header};
		my %header_data;
		my $valid = $self->parse_fasta_header_data($header, \%header_data);
		if ($valid) {
			my $name      = $header_data{name};
			my $gene_name = $header_data{gene_name};
			my $probe_id  = $name . "_$gene_name";
			my $seq       = $seq_ref->{sequence};
			
			my %probe;
			$probe{probe_name}      = $name;
			$probe{probe_gene}      = $gene_name;
			$probe{probe_id}        = $probe_id;
			$probe{sequence}        = $seq;
			
			if ($probe_ids{$probe_id}) {
				print "\n\t   - Warning: non-unique probe name '$probe_id' in probe set";
			}
			else {
				$probe_ids{$probe_id} = 1;
			}
			
			if ($type eq 'aa') {
				$probe{probe_type}  = 'ORF';
				$probe{blast_alg}   = 'tblastn';
			}
			if ($type eq 'na') {
				$probe{probe_type}  = 'UTR';
				$probe{blast_alg}   = 'blastn';
			}
			push(@$probes_ref, \%probe);	

		}
		else {
			print "\n\t Couldn't extract data from header '$header'"; 
		}
	}
	unless ($i) {
		print "\n\t No Probes were loaded";
		print "\n\t Check path is correct\n"; 
		print "\n\t Check file is in unix text format (no mac linebreaks)\n\n\n"; 
	}
}

#***************************************************************************
# Subroutine:  set_targets
# Description: set up the target files for screening
#***************************************************************************
sub set_targets {
	
	my ($self, $targets_ref, $target_groups_ref) = @_;

	# Initialise target sequence library
	my $genome_obj = TargetDB->new($self);
	my $genome_use_path  = $self->{genome_use_path};
	my $target_paths_ref = $self->{target_paths};
	unless ($target_paths_ref) { die; } 
	
	# Iterate through the list of paths 
	my %paths;
	my %target_data;	
	my @targets;	
	foreach my $path (@$target_paths_ref) {
				
		my $full_path = $genome_use_path . "/$path";	
		my $exists = $fileio->check_directory_exists($full_path);
		my @leaves;
		# If the path is a directory, get paths to all the files in it and its subdirectories
		if ($exists) {
			$fileio->read_directory_tree_leaves_simple($full_path, \@leaves);
		}
		# If the path is not to a directory, process as a file 
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
		push (@targets, @leaves);
				
	}

	# Show error and exit if no targets found
	my $targets = scalar @targets;
	unless ($targets) {
		$self->show_no_targets_found_error($target_paths_ref);
	}
	
	# Record the 'group'
	$self->set_target_groups($targets_ref, $target_groups_ref);
	
}

#***************************************************************************
# Subroutine:  show_no_targets_found_error
# Description: show 'no_targets_found' error
#***************************************************************************
sub show_no_targets_found_error {

	my ($self, $target_paths_ref) = @_;

	print "\n\n\t No target databases found";
	print " - check target paths are correctly specified in control file\n";
	print "\n\t  \$DIGS_GENOMES path is set to '$ENV{DIGS_GENOMES}'\n";
	print "\n\t  TARGETS block from control file has these paths:";
	my $i = 0;
	foreach my $path (@$target_paths_ref) {
		$i++;
		print "\n\t\t PATH 1: '$path'";
	}
	print "\n\n";
	exit;

}

#***************************************************************************
# Subroutine:  set_target_groups
# Description: Record the top-level 'group' part of the path to a target file
#***************************************************************************
sub set_target_groups {

	my ($self, $targets_ref, $target_groups_ref) = @_;
	
	# Iterate through targets
	my @keys = keys %$targets_ref;
	foreach my $target_name (@keys) {
			
		# Get target data
		my $target_ref   = $targets_ref->{$target_name};
		my $group        = $target_ref->{group};
		my $organism     = $target_ref->{organism};
		my $datatype     = $target_ref->{datatype};
		my $version      = $target_ref->{version};
		my $target_path  = $target_ref->{path};
		my $target_name  = $target_ref->{file};		
		
		# Sanity checking	
		unless ($organism and $version and $datatype and $target_name) { die; }
		
		# Create a unique key for this genome
		my @genome = ( $organism , $datatype, $version );
		my $target_id = join ('|', @genome);
		
		# Set a key to get the top level group name in the target path
		$target_groups_ref->{$target_id} = $group;
	}
}

#***************************************************************************
# Subroutine:  read genome files
# Description: processes the top level (leaves) of the genome directory
#***************************************************************************
sub read_genome_files {
	
	my ($self, $leaves_ref, $targets_ref, $exclude_ref) = @_;

	my $genome_use_path  = $self->{genome_use_path};
	my $exclude_paths_ref = $self->{exclude_paths};
	my %excluded;
	foreach my $file_ref (@$leaves_ref) {

		my $file = $file_ref->{file};
		my $path = $file_ref->{path};
		$path =~ s/\/\//\//g; # Convert any double backslashes to single		
		my $file_type = $fileio->get_infile_type($file);
			
		# Use files that are of the correct type
		if ($file_type eq 'fa' or $file_type eq 'fas' 
		or  $file_type eq 'fasta' or $file_type eq 'fna' ) {
			
			my @path = split(/\//, $path);
			my $file     = pop @path;
			my $version  = pop @path;
			my $type     = pop @path;
			my $organism = pop @path;
			my $group    = pop @path;
			unless ($organism and $type and $version) { die; }
			my @target = ( $organism , $type, $version, $file );
			my $target_id = join ('|', @target);

			my @key_path = ( $group , $organism, $type, $version, $file);
			my $key_path = join('/', @key_path);
			
			# Exclude paths from the exclude block
			my $exclude = undef;
			foreach my $path (@$exclude_paths_ref) {
				$path =~ s/_\///; # Remove a trailing forwardslash
				if ($key_path =~ m/$path/) { 
					unless ($excluded{$path}) {
						print "\n\t    Excluding: $path";
						print "\n\t    Matches:   $key_path";
						$excluded{$path} = 1;
						last;
					}
				}
				$exclude = 'true';
			}
			#$devtools->print_hash($exclude_ref); die;
			
			unless ($exclude) { 

				# Store using key
				my %data;
				$data{file}      = $file;
				$data{path}      = $path;
				$data{organism}  = $organism;
				$data{version}   = $version;
				$data{datatype}  = $type;
				$data{group}     = $group;
				$targets_ref->{$target_id} = \%data;	
				#print "\n\t STORING TARGET $path";
			}
		}
	}
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
	my $done_ref = $self->{previously_executed_searches};

	# Get relevant member variables and objects
	my $path;
	my $num_probes = scalar @$probes_ref;
	unless ($num_probes) { die "\n\t no probes found\n\n\n"; }
	
	my %target_groups; # To record top level group name in the target path (e.g. Mammals)		
	my $outstanding;   # Number of searches still to be performed
	
	# Create the queries
	foreach my $probe_ref (@$probes_ref) {
	
		my $blast_alg       = $probe_ref->{blast_alg};
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
			my $group        = $target_ref->{group};
			my $organism     = $target_ref->{organism};
			my $datatype     = $target_ref->{datatype};
			my $version      = $target_ref->{version};
			my $target_path  = $target_ref->{path};
			my $target_name  = $target_ref->{file};		
		
			# Sanity checking	
			unless ( $organism and  $version and $datatype and
                     $target_name and $probe_name and $probe_gene ) {
			 		die;
			}
			
			# Create a unique key for this genome
			my @genome = ( $organism , $datatype, $version );
			my $target_id = join ('|', @genome);

			# Create a unique key for this query
			my @key = ( $target_id, $target_name, $probe_id );
			my $key = join ('|', @key);

			if ($done_ref->{$key}) { 
				next; # Skip queries that have been issued
			} 

			# Else store the query
			$outstanding++;
			my %query = %$probe_ref;
			$query{target_id}       = $target_id;		
			$query{organism}        = $organism;		
			$query{target_version}  = $version;
			$query{target_datatype} = $datatype;		
			$query{target_name}     = $target_name;		
			$query{target_path}     = $target_ref->{path};;		

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
		
	return $outstanding;
}

############################################################################
# INPUT FILE PARSING FUNCTIONS
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
	my $output_path            = $self->{output_path};

	# Nucleic acid FASTA input	
	my $query_na_fasta         = $self->{query_na_fasta};
	my $reference_na_fasta     = $self->{reference_na_fasta};
	my $blastn_min             = $self->{bitscore_min_blastn};

	# Amino acid FASTA input
	my $query_aa_fasta         = $self->{query_aa_fasta};
	my $reference_aa_fasta     = $self->{reference_aa_fasta};
	my $tblastn_min            = $self->{bitscore_min_tblastn};

	# Track based input
	my $query_na_track         = $self->{query_na_track};
	my $query_na_track_genome  = $self->{query_na_track_genome}; 

	# Screen parameters
	my $extract_buffer         = $self->{extract_buffer};
	my $defragment_range       = $self->{defragment_range};
	
	
	if ($self->{threadhit_gap_buffer}) {
		$defragment_range = $self->{threadhit_gap_buffer}; # Deprecated
	}

	unless ($output_path) {
		print "\n\t  Warning no output path defined, results folder will be created in current directory\n\n\n";
	}
	
	# Validation for a amino acid, FASTA-based screen
	if ($query_aa_fasta) { # If a set of protein probes has been specified
		# Check if BLAST bitscore or evalue minimum set
		unless ($tblastn_min) { # Set to default minimum
			print "\n\t  Warning: no bitscore minimum defined for tblastn\n\n";
			sleep 1;
		}
		unless ($reference_aa_fasta) { # Set to default minimum
		  die "\n\t Control file error: no AA reference library defined for AA query set\n\n\n";
		}
	}
	# Validation for a nucleic acid, FASTA-based screen
	if ($query_na_fasta) {
		unless ($blastn_min) { # Set to default minimum
			print "\n\t Warning: no bitscore minimum defined for blastn\n\n";
			sleep 1;
		}
		unless ($reference_na_fasta) { # Set to default minimum
		  die "\n\t Control file error: no NT reference library defined for NT query set\n\n\n";
		}
	}
	# Validation for a track-based screen
	if ($query_na_track) {
		unless ($blastn_min) { # Set to default minimum
			print "\n\t  Warning: no bitscore minimum defined for blastn\n\n";
			sleep 1;
		}
		unless ($reference_na_fasta) { # Set to default minimum
		  die "\n\t Control file error: no NT reference library defined for NT query set\n\n\n";
		}
	}

	unless ($query_aa_fasta or $query_na_fasta or $query_na_track) {
		die "\n\t Control file error: no path to query sequence input is defined\n\n\n";
	}
	unless ($defragment_range) {
		die "\n\t Control file error: 'Screensets' block parameter 'defragment_range' is undefined. \n\n\n";
	}

}

#***************************************************************************
# Subroutine:  parse_target_block
# Description: get paths to the target sequence databases for DIGS
#***************************************************************************
sub parse_target_block {

	my ($self, $file_ref, $start, $stop, $targets_ref) = @_;

	unless ($start and $stop and $targets_ref) { die; } # Sanity checking

	# READ the 'TARGETS' block
	my @target_block;
	$fileio->extract_text_block($file_ref, \@target_block, $start, $stop);
	my $screenset_lines = scalar @target_block;
	
	# Parse the target strings
	my $targets = 0;
	my @targets;
	foreach my $line (@target_block) {
		
		chomp $line;
		#print "\n\t LINE $line";
		$line =~ s/\s+//g; # remove whitespace
		if ($line =~ /^\s*$/)   { next; } # discard blank line
		if ($line =~ /^\s*#/)   { next; } # discard comment line 
		push (@$targets_ref, $line);
		$targets++;
	}
}

#***************************************************************************
# Subroutine:  parse_nomenclature_block
# Description: read nomenclature block
# Path to file with new track(s) sorted by scaffold/chromosome & start position
# Path to file with the master (current namespace) track
# Path to translation table
# Translation system to use
# Genome ID params: Organism and assembly version
#***************************************************************************
sub parse_nomenclature_block {

	my ($self, $file_ref, $start, $stop) = @_;

	# Extract the block	
	my %params;
	$fileio->read_standard_field_value_block($file_ref, $start, $stop, \%params);
	#$devtools->print_hash(\%params);

	# Check that required parameters are set	
	my $new_track_path   = $params{new_track_path};
	my $namespace_path   = $params{namespace_path};
	my $translation_path = $params{translation_path};
	my $tax_level        = $params{tax_level};
	my $organism         = $params{nomenclature_organism};
	my $version          = $params{nomenclature_version};
	my $organism_code    = $params{organism_code};
	my $locus_class      = $params{locus_class};

	# Check what params we got
	unless ($new_track_path)      { die "\n\t No input track path defined\n\n"; }
	unless ($translation_path)    { die "\n\t No translation path defined\n\n"; }
	unless ($tax_level)           { die "\n\t No tax_level defined\n\n";        }
	unless ($organism_code)       { die "\n\t No orgaism code defined \n\n";    }
	unless ($locus_class)         { die "\n\t No locus class defined\n\n";      }
	unless ($organism)            { die "\n\t No organism defined\n\n";         }
	unless ($version)             { die "\n\t No version defined\n\n";          }

	# Set nomenclature parameters
	$self->{new_track_path}        = $new_track_path;
	$self->{namespace_path}        = $namespace_path;
	$self->{translation_path}      = $translation_path;
	$self->{tax_level}             = $tax_level;
	$self->{nomenclature_organism} = $organism;  # Genome assembly ID
	$self->{nomenclature_version}  = $version;   # Genome assembly ID
	$self->{organism_code}         = $organism_code;
	$self->{locus_class}           = $locus_class;

}

#######################################################################
# READ FASTA
############################################################################

#***************************************************************************
# Subroutine:  read_fasta
# Description: read a fasta file into an array of hashes. 
# Arguments:   $file: the name of the file to read
#              $array_ref: reference to the hash array to copy to
#***************************************************************************
sub read_fasta {

	my ($self, $file, $array_ref, $identifier) = @_;
	
	unless ($identifier) { $identifier = 'SEQ'; }

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
				my $alias_id = $alias_stem . '_' . $i;
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

############################################################################
# UTILITY FUNCTIONS ASSOCIATED  WITH SETTING UP DIGS
############################################################################

#***************************************************************************
# Subroutine:  create output directories
# Description: create a unique 'report' directory for this process
#***************************************************************************
sub create_output_directories {
	
	my ($self, $pipeline_obj) = @_;

	# Create a unique ID and report directory for this run
	my $process_id  = $self->{process_id};
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
# Subroutine:  extract_track_sequences
# Description: extract FASTA nucs from a genome assembly using an input track 
#***************************************************************************
sub extract_track_sequences {
	
	my ($self, $extracted_ref, $track_path) = @_;

	# Get paths, objects, data structures and variables from self
	my $blast_obj   = $self->{blast_obj};

	# Try to read the tab-delimited infile
	print "\n\n\t #### WARNING: This function expects a tab-delimited data table with column headers in order!";
	my $question1 = "\n\n\t Please enter the path to the file with the table data and column headings\n\n\t";
	my $query_na_track_genome_path = $console->ask_question($question1);
	
	# Read FASTA probe library
	my @track;
	$fileio->read_file($track_path, \@track);

	# Iterate through the tracks extracting
	my $i = 0;	
	my @probe_fasta;
	foreach my $line (@track) {
		
		$i++;
		
		chomp $line; # remove newline
		my @line = split("\t", $line);

		my $j;
		#foreach my $element (@line) {
		#	$j++;
		#	print "\n\t ELEMENT $j : $element"
		#}
		#die;
		
		my $name          = $line[0];
		my $scaffold      = $line[2];
		
		# TODO - hacky - resolve
		if ($scaffold =~ /Unk/) { next; }
		
		my $subject_start = $line[3];
		my $subject_end   = $line[4];
		my $gene          = $line[5];
		my $id            = $line[0];
		my $orientation;

		if ($subject_start < $subject_end) {
			$orientation = '+';
		}
		elsif ($subject_start > $subject_end) {
			$orientation = '-';
			my $switch     = $subject_start;
			$subject_start = $subject_end;
			$subject_end   = $switch;
		}
		
		# Extract the sequence
		my %data;
		$data{subject_start} = $subject_start;
		$data{subject_end}   = $subject_end;
		$data{orientation}   = $orientation;
		$data{scaffold}      = $scaffold;

		my $target_path = "$query_na_track_genome_path" . "/$scaffold" . '.fa';;	
		my $sequence = $blast_obj->extract_sequence($target_path, \%data);
		if ($sequence) {	
			#print "\n\t got seq $sequence \n\n";
			my %probe;
			$probe{probe_name}      = $name;
			$probe{probe_gene}      = $gene;
			$probe{probe_id}        = $name . "_$gene";
			$probe{sequence}        = $sequence;
			$probe{probe_type}      = 'UTR';
			$probe{blast_alg}       = 'blastn';
			
			push(@$extracted_ref, \%probe);
			
			my $header = "$name" . "_$gene";
			$header =~ s/\(/\./g;
			$header =~ s/\)//g;

			print "\n\t\t Getting probe $i: $header";

			my $digs_fasta = ">$header" . "\n$sequence\n";
			push (@probe_fasta, $digs_fasta);;
		
		}
		else {
			die "\n\t Sequence extraction failed";
		}
	}	
	
	my $outfile = 'extracted.DIGS.fna';
	$fileio->write_file($outfile, \@probe_fasta);
}

############################################################################
# EOF
############################################################################
