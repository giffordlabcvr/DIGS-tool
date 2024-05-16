#!/usr/bin/perl -w
############################################################################
# Module:      ScreenBuilder.pm
# Description: Functions for setting up a screen using the DIGS tool
#              based on settings received in a DIGS file
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
		report_dir           => $parameter_ref->{report_dir},
		
		# Member classes 
		blast_obj            => $parameter_ref->{blast_obj},
		
		# Data
		previously_executed_searches => $parameter_ref->{previously_executed_searches},

	};
	bless ($self, $class);
	return $self;
}

############################################################################
# SETUP FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  setup_screen
# Description: Set up all the queries to execute, indexed by genome target
#***************************************************************************
sub setup_screen  {
	
	my ($self, $digs_obj, $queries_ref) = @_;
	
	# Import the probes
	my @probes;
	$self->setup_blast_probes(\@probes);
	
	# Set up the reference library for BLAST
	$self->setup_reference_libraries($digs_obj);
	
	# Set target sequence files for screening
	my %targets;
	my $num_targets = $self->set_targets(\%targets);

	# Show error and exit if no targets found
	unless ($num_targets) {
		$self->show_no_targets_found_error();
	}

	# Record the target group designations
	my %target_groups;
	$self->set_target_groups(\%targets, \%target_groups);
	$digs_obj->{target_groups} = \%target_groups; 

	# Show output about the number of oustanding targets and number of queries loaded
	my @keys        = keys %targets;
	my $unique      = scalar @keys;
	print "\n\t  Targets:           $unique target files";

	# Create the BLAST queries for this screen
	my $num_queries = $self->set_queries($digs_obj, \@probes, \%targets, $queries_ref);

	unless ($num_queries) {
		print "\n\n\t ### No outstanding queries were loaded\n";
	}
	else { print "\n\t  Searches to run    $num_queries\n"; }

	return $num_queries;
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
		#$devtools->print_hash($self);
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
	unless ($num_fasta) {
		print "\n\t No Probes were loaded";
		print "\n\t Check path is correct\n"; 
		print "\n\t Check file is in unix text format (no mac linebreaks)\n\n\n"; 
	}
	
	# Show probes if verbose flag set
	if ($self->{verbose}) {
		
		#$devtools->print_array(\@fasta);
		my $i;
		print "\n\n";
		foreach my $seq_ref (@fasta) {
			
			$i++;	
			my $header = $seq_ref->{header};
			my $seq_id = $seq_ref->{sequence_id};
			my $sequence = $seq_ref->{sequence};
			print "\n\t sequence $i: $header (alias $seq_id)";
		
		}
		print "\n\n";
	
	}	
	
	print "\n\t  Probe sequences:   $num_fasta $probe_type sequences\n\n";
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
				print "\t   - Warning: non-unique probe name '$probe_id' in probe set\n";
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
		print "\n\t No Probes were loaded - please check FASTA header format";
		print "\n\t Check path is correct\n"; 
		print "\n\t DIGS probes should be in FASTA format\n\n\n"; 
	}
}

#***************************************************************************
# Subroutine:  setup_reference_libraries
# Description: handler for reference library set up
#***************************************************************************
sub setup_reference_libraries {

	my ($self, $digs_obj, $consolidate) = @_;

	# DEV $devtools->print_hash($digs_obj); die;
	
	# Get reference library params
	my $reference_type;
	my $ref_fasta;
	my @ref_fasta;
	my $num_fasta;

	if ($consolidate) {
	
		unless ($digs_obj->{consolidated_reference_aa_fasta} or $digs_obj->{consolidated_reference_na_fasta} ) {
			die "\n\t No library defined for consolidate function"
		}
	}
		
	# Load a peptide sequence library
	if ($digs_obj->{reference_aa_fasta} or $digs_obj->{consolidated_reference_aa_fasta}) {
 		$self->{reference_library_type} = 'aa';
 		unless ($consolidate) {
 			$ref_fasta = $digs_obj->{reference_aa_fasta};
		}
		else {
 			$ref_fasta = $digs_obj->{consolidated_reference_aa_fasta};
		}
		$reference_type = 'amino acid';
		$self->create_reference_library($digs_obj, $ref_fasta, $reference_type);
	}
	
	# Load a nucleotide sequence library	
	if ($digs_obj->{reference_na_fasta} or $digs_obj->{consolidated_reference_na_fasta}) {
 		$self->{reference_library_type} = 'na';
 		unless ($consolidate) {
 			$ref_fasta = $digs_obj->{reference_na_fasta};
		}
		else {
 			$ref_fasta = $digs_obj->{consolidated_reference_na_fasta};
		}
		$reference_type = 'nucleic acid';
		$self->create_reference_library($digs_obj, $ref_fasta, $reference_type);
	}

	#$devtools->print_hash($digs_obj); die;

}

#***************************************************************************
# Subroutine:  create_reference_library
# Description: Set up the references for BLAST-based classification
#***************************************************************************
sub create_reference_library {
	
	my ($self, $digs_obj, $ref_fasta, $reference_type) = @_;

	unless ($ref_fasta) { 
		print "\n\t  Reference library: $reference_type FASTA not found'\n\n"; 
		print "\n\t  If you're running a nucleotide search, check that all amino acid (aa) related variables are commented out'\n\n"; 
		die;
	}
	
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
			if ($self->{verbose}) {
				print "\n\t\t  Warning: non-unique reference name '$refseq_id'";
			}
		}
		else {
			$refseq_ids{$refseq_id} = 1;
		}
		
		my $fasta = ">$refseq_id" . "\n$seq\n\n";
		push (@references, $fasta);
	}

	if ($nonunique) {
		print "\n\n\t  ### Warning: $nonunique non-unique names identified in reference library\n";
		sleep 1;
		#print "\n\t\t  Warning: non-unique reference name '$refseq_id'";
	}

	# Format the library for BLAST
	if ($num_fasta) {
		if ($self->{reference_library_type} eq 'aa') {
			$self->create_blast_lib(\@references, 'aa');
			$digs_obj->{aa_reference_library} = $lib_file;
		}
		if ($self->{reference_library_type} eq 'na') {
			$self->create_blast_lib(\@references, 'na');
			$digs_obj->{na_reference_library} = $lib_file;
		}
	}

	# Set the paths to the BLAST-formatted libraries
	$digs_obj->{blast_utr_lib_path} = $self->{blast_utr_lib_path};
	$digs_obj->{blast_orf_lib_path} = $self->{blast_orf_lib_path};

	#$devtools->print_hash($digs_obj); die;

}

#***************************************************************************
# Subroutine:  set_targets
# Description: set up the target files for screening
#***************************************************************************
sub set_targets {
	
	my ($self, $targets_ref) = @_;

	# Initialise target sequence library
	my $genome_obj = TargetDB->new($self);
	my $genome_use_path   = $self->{genome_use_path};
	my $target_paths_ref  = $self->{target_paths};
	unless ($target_paths_ref) { $devtools->print_hash($self); die; }  # Sanity checking
	
	# Iterate through the list of paths 
	my %paths;
	my %target_data;	
	my $num_targets;
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
		my $num_leaves = $self->read_genome_files(\@leaves, $targets_ref);

		# Keep count
		if ($num_leaves) {
			if ($num_targets) { 
				$num_targets = $num_targets + $num_leaves;
			}
			else {
				$num_targets = $num_leaves;
			}
		}
	}
	return $num_targets;
}

#***************************************************************************
# Subroutine:  set_target_groups
# Description: Record the top-level 'group' part of the path to a target file
#***************************************************************************
sub set_target_groups {

	my ($self, $targets_ref, $target_groups_ref) = @_;
	
	# Iterate through targets
	my @keys = keys %$targets_ref;
	my %shown;
	foreach my $target_id (@keys) {
			
		# Get target data
		my $target_ref   = $targets_ref->{$target_id};
		my $group        = $target_ref->{group};
		my $organism     = $target_ref->{organism};
		my $datatype     = $target_ref->{datatype};
		my $version      = $target_ref->{version};
		my $target_path  = $target_ref->{path};
		my $target_name  = $target_ref->{file};		
		
		# Sanity checking	
		unless ($organism and $version and $datatype and $target_name) { die; }
		
		# Create a unique key for this genome
		#my @genome = ( $organism , $datatype, $version );
		#my $new_target_id = join ('|', @genome);
		#print "\n\t ### old '$target_id' \n\t ### new '$new_target_id'";
		
		# Set a key to get the top level group name in the target path
		unless ($shown{$organism}) {
			$shown{$organism} = 1;
			if ($self->{verbose}) {
				print "\n\t  Based on your target directory structure, organism $organism is a member of group '$group'";
	
			}
		}
		$target_groups_ref->{$target_id} = $group;
	}
}

#***************************************************************************
# Subroutine:  read_genome_files
# Description: select target files from an array 
#***************************************************************************
sub read_genome_files {
	
	my ($self, $files_array_ref, $targets_ref) = @_;

	my $exclude_paths_ref = $self->{exclude_paths};
	my $verbose           = $self->{verbose};
	#$devtools->print_hash($exclude_paths_ref); die; # DEBUG
	
	my $count = 0;
	foreach my $file_ref (@$files_array_ref) {

		# Test whether this file has a FASTA file extension
		my $file = $file_ref->{file};
		my $is_fasta = $self->does_file_have_fasta_extension($file);		
		unless ($is_fasta) {
			next; # Skip everything that isn't explicitly labelled as FASTA
		}

		my $path = $file_ref->{path};
		$path =~ s/\/\//\//g; # Convert any double backslashes to single		
		$path =~ s/\/\//\//g; # Convert any double backslashes to single
			
		my @path = split(/\//, $path);
		pop @path;
		my $version  = pop @path;
		my $type     = pop @path;
		my $organism = pop @path;
		my $group    = pop @path;
		unless ($organism and $type and $version) { die; }
		my @target = ( $organism , $type, $version, $file );
		my $target_id = join ('|', @target);

		my @key_path = ( $group , $organism, $type, $version, $file);
		my $key_path = join('/', @key_path);
		#print "\n\t ######  PATH $key_path";
		unless ($exclude_paths_ref->{$key_path}) {

			if ($verbose) {
				#print "\n\t\t  Target '$key_path' added";
			}

			# Store using target_id as a key
			my %data;
			$data{file}      = $file;
			$data{path}      = $path;
			$data{organism}  = $organism;
			$data{version}   = $version;
			$data{datatype}  = $type;
			$data{group}     = $group;
			$targets_ref->{$target_id} = \%data;	
			$count++;

		}
		elsif ($verbose) {
			print "\n\t\t  Target '$key_path' EXCLUDED";
		}				
	}
	return $count;
}

#***************************************************************************
# Subroutine:  show_no_targets_found_error
# Description: show 'no_targets_found' error
#***************************************************************************
sub show_no_targets_found_error {

	my ($self) = @_;

	my $target_paths_ref = $self->{target_paths};
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
# Subroutine:  does_file_have_fasta_extension
# Description: use file extension to determine if a file is FASTA format
#***************************************************************************
sub does_file_have_fasta_extension {
	
	my ($self, $file_name) = @_;

	my $has_fasta_extension = undef;

	my $file_extension = $fileio->get_infile_type($file_name);	
	if ($file_extension) {

		# Make extension lowercase (in case it includes capitals)
		$file_extension = lc $file_extension;
		
		# Test if the extension is one of a few possible FASTA ones				
		if ($file_extension eq 'fa' 
		or  $file_extension eq 'fas' 
		or  $file_extension eq 'fasta'
		or  $file_extension eq 'faa' 
		or  $file_extension eq 'fna' ) {
			 $has_fasta_extension = 'true';
		} 
	}
	return $has_fasta_extension;
}

#***************************************************************************
# Subroutine:  set_queries
# Description: set up the individual BLAST searches 
#***************************************************************************
sub set_queries {
	
	my ($self, $digs_obj, $probes_ref, $targets_ref, $queries_ref) = @_;

	my $db = $digs_obj->{db};
	unless ($db) { die; }
	
	# Get data from self
	my $report_dir = $digs_obj->{report_dir};
	my $tmp_path   = $digs_obj->{tmp_path};
	unless ($report_dir) { die; } 
	unless ($tmp_path)   { die; } 

	# Get the target database information
	my @target_names = sort keys %$targets_ref;
	my $num_targets = scalar @target_names;
	unless ($num_targets) { die "\n\t No target databases found\n\n\n";	}

	# Work out the current state with respect to searches performed
	my $done_ref = $self->{previously_executed_searches};
	my @query_keys = keys %$done_ref;
	my $num_previous_queries = scalar @query_keys;
	print "\n\t  Previous queries:  $num_previous_queries previous queries";

	# Get relevant member variables and objects
	my $path;
	my $num_probes = scalar @$probes_ref;
	unless ($num_probes) { die "\n\t no probes found\n\n\n"; }
	
	my %target_groups; # To record top level group name in the target path (e.g. Mammals)		
	my $outstanding;   # Number of searches still to be performed
	
	# Create the queries
	my $skipped_because_done = '0';
	my $num_queries = '0';
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
			
			$num_queries++;			

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
				$skipped_because_done++;
				next; # Skip queries that have already been performed
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
	
	print "\n\t  Skipped in set:    $skipped_because_done (of $num_queries)";
	return $outstanding;
}

############################################################################
# INPUT FILE PARSING FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  parse control file
# Description: read an input file to get parameters for screening
#***************************************************************************
sub parse_control_file {

	my ($self, $ctl_file, $digs_obj) = @_;
	
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
	my %exclude;
	$self->set_exclude_targets(\@exclude, \%exclude);
	$self->{exclude_paths} = \%exclude;

	# READ the 'SKIPINDEX' block if present (Targets to skip during indexing)
	my @skipindex;
	$start_token = 'BEGIN SKIPINDEX';
	$stop_token  = 'ENDBLOCK';
	$self->parse_target_block(\@ctl_file, $start_token, $stop_token, \@skipindex);
	my %skipindex;
	$self->set_skipindex_targets(\@skipindex, \%skipindex);
	#$devtools->print_hash(\%skipindex); die;
	
	$self->{skipindexing_paths} = \%skipindex;
	#$devtools->print_hash($self); die;
	
	# Set parameters for DIGS based on parameters stored in $self after after parsing
	# Screening DB name and MySQL connection details
	$digs_obj->{db_name}                = $self->{db_name};
	$digs_obj->{mysql_server}           = $self->{mysql_server};

	# Input and output file paths 	
	$digs_obj->{output_path}            = $self->{output_path};
	$digs_obj->{target_paths}           = $self->{target_paths};
	$digs_obj->{skipindexing_paths}     = $self->{skipindexing_paths};

	# Set paths to query/probe files
	$digs_obj->{query_na_fasta}         = $self->{query_na_fasta};
	$digs_obj->{reference_na_fasta}     = $self->{reference_na_fasta};
	$digs_obj->{query_aa_fasta}         = $self->{query_aa_fasta};
	$digs_obj->{reference_aa_fasta}     = $self->{reference_aa_fasta};

	# Set screening extract and defragment parameters
	$digs_obj->{defragment_range}       = $self->{defragment_range};
	$digs_obj->{extract_buffer}         = $self->{extract_buffer};

	# Set minimum length thresholds for extracting hits
	$digs_obj->{seq_length_minimum}     = $self->{seq_length_minimum};

	# Set the bit score minimum for extracting hits in a tBLASTn-based forward screen
	if ($self->{bitscore_min_tblastn}) {
		$digs_obj->{bitscore_minimum} = $self->{bitscore_min_tblastn};
	}
	# Set the bit score minimum for extracting hits in a BLASTn-based forward screen
	if ($self->{bitscore_min_blastn}) {
		$digs_obj->{bitscore_minimum} = $self->{bitscore_min_blastn};	
	}

	# Set parameters for forward BLAST (probe versus target database) 
	my $fwd_num_threads = $self->{fwd_num_threads};
	unless ($fwd_num_threads)  { $fwd_num_threads = 1; }  # Set to one  setting
	$digs_obj->{fwd_num_threads} = $fwd_num_threads;
	$digs_obj->{fwd_word_size}   = $self->{fwd_word_size};
	$digs_obj->{fwd_evalue}      = $self->{fwd_evalue};
	$digs_obj->{fwd_penalty}     = $self->{fwd_penalty};
	$digs_obj->{fwd_reward}      = $self->{fwd_reward};
	$digs_obj->{fwd_gapopen}     = $self->{fwd_gapopen};
	$digs_obj->{fwd_gapextend}   = $self->{fwd_gapextend};
	$digs_obj->{fwd_dust}        = $self->{fwd_dust};
	$digs_obj->{fwd_softmasking} = $self->{fwd_softmasking};
	$digs_obj->{fwd_seg}         = $self->{fwd_seg};

	# Set parameters for forward BLAST (probe versus target database) 
	my $rev_num_threads = $self->{rev_num_threads};
	unless ($rev_num_threads)  { $rev_num_threads = 1; }  # Set to one  setting
	$digs_obj->{rev_num_threads} = $rev_num_threads;
	$digs_obj->{rev_word_size}   = $self->{rev_word_size};
	$digs_obj->{rev_evalue}      = $self->{rev_evalue};
	$digs_obj->{rev_penalty}     = $self->{rev_penalty};
	$digs_obj->{rev_reward}      = $self->{rev_reward};
	$digs_obj->{rev_gapopen}     = $self->{rev_gapopen};
	$digs_obj->{rev_gapextend}   = $self->{rev_gapextend};
	$digs_obj->{rev_dust}        = $self->{rev_dust};
	$digs_obj->{rev_softmasking} = $self->{rev_softmasking};
	$digs_obj->{rev_seg}         = $self->{rev_seg};

	# Set parameters for consolidation step
	$digs_obj->{consolidate_range}                = $self->{consolidate_range};
	$digs_obj->{consolidated_reference_aa_fasta}  = $self->{consolidated_reference_aa_fasta};
	$digs_obj->{consolidated_reference_na_fasta}  = $self->{consolidated_reference_na_fasta};
	
	# Capture list of any target databases to skip when formatting for BLAST
	$digs_obj->{skipindexing_paths} = \%skipindex;
	
	# DEV $devtools->print_hash($digs_obj); die;

}

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

	unless ($db_name)  {
		die "\n\t Control file error: 'db_name' undefined in 'SCREENDB' block\n\n\n";
	}
	unless ($server)  {
		die "\n\t Control file error: 'mysql_server' undefined in 'SCREENDB' block\n\n\n";
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
	my $consolidate_range      = $self->{consolidate_range};

	
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
	my @targets;
	foreach my $line (@target_block) {
		
		chomp $line;
		$line =~ s/\s+//g; # remove whitespace
		if ($line =~ /^\s*$/)   { next; } # discard blank line
		if ($line =~ /^\s*#/)   { next; } # discard comment line 
		push (@$targets_ref, $line);
	}	
}

#***************************************************************************
# Subroutine:  set_exclude_targets
# Description: set up hash to record paths of targets to exclude from screening
#***************************************************************************
sub set_exclude_targets {
	
	my ($self, $exclude_array_ref, $exclude_hash_ref) = @_;

	# Iterate through the list of paths 
	my %paths;
	my $genome_use_path = $self->{genome_use_path};
	unless ($genome_use_path) { die; }
	
	foreach my $path (@$exclude_array_ref) {
				
		my $full_path = $genome_use_path . "/$path";	
		my $exists = $fileio->check_directory_exists($full_path);

		my @leaves;
		# If the path is a directory, get paths to all the files in it and its subdirectories
		if ($exists) {
			$fileio->read_directory_tree_leaves_simple($full_path, \@leaves);
		}
		# If the path is not to a directory, process as a file 
		else {
		
			$path =~ s/\/\//\//g; # Convert any double backslashes to single		
			$path =~ s/\/\//\//g; # Convert any double backslashes to single		
			
			my @path = split(/\//, $path);
			my $file = pop @path;
			my %file;
			$file{file} = $file;
			$file{path} = $path;
			push (@leaves, \%file);
		}

		# Record in a hash if its a FASTA file
		foreach my $leaf (@leaves) {
			my $path = $leaf->{path};
			$path =~ s/\/\//\//g; # Convert any double backslashes to single		
			$path =~ s/\/\//\//g; # Convert any double backslashes to single		
			my @path = split(/\//, $path);
			my $file = pop @path;
			my $is_fasta = $self->does_file_have_fasta_extension($file);		
			if ($is_fasta) {
				
				my $version  = pop @path;
				my $type     = pop @path;
				my $organism = pop @path;
				my $group    = pop @path;
				my @key_path = ( $group , $organism, $type, $version, $file);
				my $key_path = join('/', @key_path);
				$exclude_hash_ref->{$key_path} = 1;
			}
		}		
	}
}

#***************************************************************************
# Subroutine:  set_skipindex_targets
# Description: set up hash to record paths of targets to exclude from screening
#***************************************************************************
sub set_skipindex_targets {
	
	my ($self, $skipindex_array_ref, $skipindex_hash_ref) = @_;

	# Iterate through the list of paths 
	my %paths;
	my $genome_use_path = $self->{genome_use_path};
	unless ($genome_use_path) { die; }
	
	foreach my $path (@$skipindex_array_ref) {
			
		my $full_path = $genome_use_path . "/$path";	
		my $exists = $fileio->check_directory_exists($full_path);

		my @leaves;
		# If the path is a directory, get paths to all the files in it and its subdirectories
		if ($exists) {
			$fileio->read_directory_tree_leaves_simple($full_path, \@leaves);
		}
		# If the path is not to a directory, process as a file 
		else {
		
			$path =~ s/\/\//\//g; # Convert any double backslashes to single		
			$path =~ s/\/\//\//g; # Convert any double backslashes to single		
			
			my @path = split(/\//, $path);
			my $file = pop @path;
			my %file;
			$file{file} = $file;
			$file{path} = $path;
			push (@leaves, \%file);
		}

		# Record in a hash if its a FASTA file
		foreach my $leaf (@leaves) {
			my $path = $leaf->{path};
			$path =~ s/\/\//\//g; # Convert any double backslashes to single		
			$path =~ s/\/\//\//g; # Convert any double backslashes to single		
			my @path = split(/\//, $path);
			my $file = pop @path;
			my $is_fasta = $self->does_file_have_fasta_extension($file);		
			if ($is_fasta) {
				
				my $version  = pop @path;
				my $type     = pop @path;
				my $organism = pop @path;
				my $group    = pop @path;
				my @key_path = ( $group , $organism, $type, $version);
				my $key_path = join('/', @key_path);
				$skipindex_hash_ref->{$key_path} = 100;
			}
		}		
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

	my ($self, $digs_obj, $file_ref, $start, $stop) = @_;

	return;
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

	# Set paths for applying nomenclature 	
	#$digs_obj->{new_track_path}        = $self->{new_track_path};
	#$digs_obj->{translation_path}      = $self->{translation_path};
	#$digs_obj->{tax_level}             = $self->{tax_level};
	#$digs_obj->{organism_code}         = $self->{organism_code};
	#$digs_obj->{locus_class}           = $self->{locus_class};
	#$digs_obj->{nomenclature_version}  = $self->{nomenclature_version};
	#$digs_obj->{nomenclature_organism} = $self->{nomenclature_organism};
	#$digs_obj->{genome_structure}      = $self->{genome_structure};

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
				$header = $self->clean_fasta_header($header);
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
		$header = $self->clean_fasta_header($header);
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
# Subroutine:  clean_fasta_header
# Description: 
#***************************************************************************
sub clean_fasta_header {
	
	my ($self, $header) = @_;

	# Remove illegal characters from the header line: these include:
	# / : * ? " < > |   because we need to write files using header elements
	# '                 because quotes interfere with SQL statements
	$header =~ s/\|//g;
	$header =~ s/\(//g;
	$header =~ s/\)//g;
	$header =~ s/\///g;
	$header =~ s/\*//g;
	$header =~ s/\?//g;
	$header =~ s/://g;
	$header =~ s/"//g;
	$header =~ s/<//g;
	$header =~ s/>//g;
	$header =~ s/'//g;
	$header =~ s/\s+//g;

	return $header;
}

############################################################################
# UTILITY FUNCTIONS ASSOCIATED  WITH SETTING UP DIGS
############################################################################

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

	# Execute command
	my $makedb_cmd;
	if ($type eq 'aa') {
		$makedb_cmd = "$blast_program -in $lib_path -dbtype prot > /dev/null";
		$self->{blast_orf_lib_path} = $lib_path; 
	}
	elsif ($type eq 'na') {
		$makedb_cmd = "$blast_program -in $lib_path -dbtype nucl> /dev/null";
		$self->{blast_utr_lib_path} = $lib_path;
	}
	my $result = system $makedb_cmd;
	if ($result) {
		print "\n\t Failed to format reference library for BLAST! \n\n";
		print "\n\t $makedb_cmd \n\n"; exit;	
	}
}

############################################################################
# EOF
############################################################################

