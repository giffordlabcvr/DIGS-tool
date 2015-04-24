#!/usr/bin/perl -w
############################################################################
# Module:      ScreenBuild.pm
# Description: Set up a paired BLAST screen in the DIGS framework
# History:     December 2012: Created by Robert Gifford 
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
use DIGS::ScreeningDB;

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
		
		# Paths
		genome_use_path      => $parameter_ref->{genome_use_path},
		blast_bin_path       => $parameter_ref->{blast_bin_path},
		
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
	
	my ($self, $pipeline_obj, $queries_ref) = @_;

	# Create the output directories
	$self->create_output_directories($pipeline_obj);
	
	# Set up the reference library
	$self->setup_reference_library($pipeline_obj);

	# Set up the probes
	my @probes;
	$self->setup_blast_probes(\@probes);

	# Set up the target sequences to screen
	my %targets;
	$self->set_targets(\%targets);
	
	# Create the list of BLAST queries for screening
	my $num = $self->set_queries($pipeline_obj, \@probes, \%targets, $queries_ref);
	return $num; # Return the number of paired BLAST queries created
}

############################################################################
# SET REFERENCE LIBRARY FOR RECIPROCAL BLAST
############################################################################

#***************************************************************************
# Subroutine:  setup_reference_library
# Description: handler for reference library set up
#***************************************************************************
sub setup_reference_library {
	
	my ($self, $pipeline_obj) = @_;

	if ($self->{reference_aa_fasta}) {
		$self->load_aa_fasta_reference_library();
	}
	elsif ($self->{reference_na_fasta}) {
		$self->load_nt_fasta_reference_library();
	}
	elsif ($self->{reference_glue}) {
		$self->load_glue_reference_library();
	}
	else { 
		die "\n\t No reference sequences found for screen\n\n\n";
	}
	# Set the paths to the BLAST-formatted libraries
	$pipeline_obj->{blast_utr_lib_path} = $self->{blast_utr_lib_path};
	$pipeline_obj->{blast_orf_lib_path} = $self->{blast_orf_lib_path};
}

#***************************************************************************
# Subroutine:  load_aa_fasta_reference_library
# Description: load a library of protein sequences in FASTA format as references 
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
		print "\n\n\t   '$num_fasta' FASTA formatted protein sequences in reference library";
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
				my $aa_seq    = $seq_ref->{sequence};
				my $fasta = ">$name" . "_$gene_name" . "\n$aa_seq\n\n";
				push (@ref_aa_fasta, $fasta);
			}
		}
	}

	# Create the libraries
	if ($num_fasta) { $self->create_blast_aa_lib(\@ref_aa_fasta); }
}

#***************************************************************************
# Subroutine:  load_nt_fasta_reference_library
# Description: load a library of nucleotide sequences in FASTA format as references 
#***************************************************************************
sub load_nt_fasta_reference_library {
	
	my ($self) = @_;

	# Get reference library params
	my $ref_nt_fasta = $self->{reference_na_fasta};
	
	# Format a reference NT library	
	print "\n\t ### Loading FASTA reference sequences";
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
		print "\n\n\t   '$num_fasta' FASTA formatted nucleotide sequences in reference library";
		my $i = 0;
		foreach my $seq_ref (@fasta) {
			$i++;
			my $header  = $seq_ref->{header};
			my %header_data;
			my $mode = $self->{ref_fasta_header_mode};
			my $valid = $self->parse_fasta_header_data($header, \%header_data);
			if ($valid) {
				my $name     = $header_data{name};
				my $gene_name = $header_data{gene_name};
				my $nt_seq   = $seq_ref->{sequence};
				my $fasta = ">$name" . "_$gene_name" . "\n$nt_seq\n\n";
				push (@ref_nt_fasta, $fasta);
			}
		}
	}
	
	# Create the libraries
	if ($num_fasta) { $self->create_blast_nt_lib(\@ref_nt_fasta); }
}

#***************************************************************************
# Subroutine:  load_glue_reference_library 
# Description: 
#***************************************************************************
sub load_glue_reference_library {
	
	my ($self, $queries_ref) = @_;

	my $report_dir  = $self->{output_path};
	my $refseq_lib = RefSeqLibrary->new($self);
	my %refseq_library;
	my $path = $self->{reference_glue};
	$refseq_lib->load_glue_reference_library(\%refseq_library, $path);

	# Set up library
	my $parser_obj  = RefSeqParser->new();
	my $i = 0;
	my @refseq_names = keys %refseq_library;
	my $num_refseqs = scalar @refseq_names;
	print "\n\n\t GLUE reference library: '$num_refseqs' sequences";
	my @ref_fasta_nt;
	my @ref_fasta_aa;
	foreach my $refseq_name (@refseq_names) {
		
		$i++;
		my $refseq = $refseq_library{$refseq_name};	
	
		# Get the translated ORFs
		my %orf_sequences;
		$refseq->get_translated_orfs(\%orf_sequences); # Get orfs nucleic acid seqs 
		
		# Get the untranslated ORFs
		my %utr_sequences;
		$refseq->get_utrs(\%utr_sequences); # Get orfs nucleic acid seqs 
		#$devtools->print_hash(\%utr_sequences); die;

		my @orf_names = keys %orf_sequences;
		foreach my $orf_name (@orf_names) {
			my $orf_seq = $orf_sequences{$orf_name};
			my $name = $refseq_name . "_$orf_name";
			my $fasta = ">$name\n$orf_seq\n\n";
			push (@ref_fasta_aa, $fasta);
		}
		
		# Iterate through UTRs and add those
		my @utr_names = keys %utr_sequences;
		foreach my $utr_name (@utr_names) {
			my $utr_seq = $utr_sequences{$utr_name};
			my $name = $refseq_name . "_$utr_name";
			my $fasta = ">$name\n$utr_seq\n\n";
			my $utr_len = length $utr_seq;
			unless ($utr_len) { next; }
			push (@ref_fasta_nt, $fasta);
		}
	}

	# Create the libraries
	my $num_ref_fasta_nt = 0;
	my $num_ref_fasta_aa = 0;
	$num_ref_fasta_nt = scalar @ref_fasta_nt;
	$num_ref_fasta_aa = scalar @ref_fasta_aa;
	print "\n\n\t '$num_ref_fasta_nt' NA references loaded";
	print "\n\t '$num_ref_fasta_aa' AA references loaded\n\n";

	if ($num_ref_fasta_nt) { $self->create_blast_nt_lib(\@ref_fasta_nt); }
	if ($num_ref_fasta_aa) { $self->create_blast_aa_lib(\@ref_fasta_aa); }
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
	
	# Copy file to the report directory
	my $aa_lib_path = $report_dir . "/reference_lib_aa.fas";
	$fileio->write_file($aa_lib_path, $aa_lib_ref); 
	
	# Set path to blast binary
	my $blast_program = 'makeblastdb';
	my $blast_bin_dir = $self->{blast_bin_path};
	my $bin_path;
	if ($blast_bin_dir) {
		 $bin_path = $self->{blast_bin_path} . $blast_program;
	}
	else { $bin_path = $blast_program; }

	# Execute command
	my $makedb_cmd = "$bin_path -in $aa_lib_path -dbtype prot > /dev/null";
	my $result = system $makedb_cmd;
	if ($result) {
		#print "\n\t $makedb_cmd \n\n"; die;	
		die "\n\t Failed to format reference library for BLAST! \n\n";
	}
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

	# Copy file to the report directory
	my $nt_lib_path = $report_dir . "/reference_lib_nt.fas";
	$fileio->write_file($nt_lib_path, $nt_lib_ref);

	# Set path to blast binary
	my $blast_program = 'makeblastdb';	
	my $blast_bin_dir = $self->{blast_bin_path};
	my $bin_path;
	if ($blast_bin_dir) {
		$bin_path = $self->{blast_bin_path} . $blast_program;
	}
	else { $bin_path = $blast_program; }

	# Execute command
	my $makedb_cmd = "$bin_path -in $nt_lib_path -dbtype nucl> /dev/null";
	#print $makedb_cmd;
	my $result = system $makedb_cmd;
	if ($result) {
		#print "\n\t $makedb_cmd \n\n"; die;	
		die "\n\t Failed to format reference library for BLAST! \n\n";
	}
	$self->{blast_utr_lib_path} = $nt_lib_path;
}

############################################################################
# SET UP PROBES
############################################################################

#***************************************************************************
# Subroutine:  setup_blast_probes
# Description: handler for setting up probes for a screen
#***************************************************************************
sub setup_blast_probes {
	
	my ($self, $probes_ref) = @_;

	print "\n\n\t ### Setting up the sequence 'probes' for screening";
	if ($self->{query_aa_fasta}) {
		$self->load_aa_fasta_probes($probes_ref);
	}
	elsif ($self->{query_na_fasta}) {
		$self->load_nt_fasta_probes($probes_ref);
	}
	elsif ($self->{query_glue}) {
		$self->load_glue_query($probes_ref);
	}
	else { 
		die "\n\t No probes found for screen\n\n\n";
	}
	#$devtools->print_array(\@probes); die;	# DEBUG
}

#***************************************************************************
# Subroutine:  load_glue_query
# Description: 
#***************************************************************************
sub load_glue_query {
	
	my ($self, $probes_ref) = @_;

	my $refseq_lib = RefSeqLibrary->new($self);
	my %refseq_library;
	my $path = $self->{query_glue};
	$refseq_lib->load_glue_reference_library(\%refseq_library, $path);

	# Get cutoffs
	my $bit_score_min_tblastn = $self->{bit_score_min_tblastn};
	my $bit_score_min_blastn  = $self->{bit_score_min_blastn};
	
	my $i = 0;
	my @query_fasta;
	my @na_query_fasta;
	my @refseqs = keys %refseq_library;
	foreach my $refseq_name (@refseqs) {
		
		my $refseq = $refseq_library{$refseq_name};

		# Get the translated ORFs
		my %orf_sequences;
		$refseq->get_translated_orfs(\%orf_sequences); # Get orfs nucleic acid seqs 
		
		# Get the untranslated regions
		my %utr_sequences;
		$refseq->get_utrs(\%utr_sequences); # Get orfs nucleic acid seqs 
		#$devtools->print_hash(\%utr_sequences); die;

		# ADD THE ORF SEARCHES
		my @orf_names = keys %orf_sequences;
		unless ($bit_score_min_tblastn) { die; }
		foreach my $orf_name (@orf_names) {
			my $aa_seq   = $orf_sequences{$orf_name};
			$self->add_aa_probe($probes_ref, $refseq_name, $orf_name, $aa_seq);
		}
	
		# Iterate through UTRs and add those
		my @utr_names = keys %utr_sequences;
		unless ($bit_score_min_blastn) { die; }
		foreach my $utr_name (@utr_names) {
			my $utr_seq = $utr_sequences{$utr_name};
			$self->add_na_probe($probes_ref, $refseq_name, $utr_name, $utr_seq);
		}
	}
	#$devtools->print_array($probes_ref); die;	# DEBUG
}

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
		print "\n\n\t   '$num_fasta' FASTA formatted protein sequences will be used as probes";
		my $i = 0;
		foreach my $seq_ref (@fasta) {
			$i++;
			my $header  = $seq_ref->{header};
			my %header_data;
			my $valid = $self->parse_fasta_header_data($header, \%header_data);
			if ($valid) {
				my $name     = $header_data{name};
				my $gene_name = $header_data{gene_name};
				my $aa_seq   = $seq_ref->{sequence};
				#$devtools->print_hash(\%header_data);
				$self->add_aa_probe($probes_ref, $name, $gene_name, $aa_seq);
			}
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
	my $query_na_fasta = $self->{query_na_fasta};

	# Read FASTA nt probe library
	if ($query_na_fasta) {
		my @fasta;
		$seqio->read_fasta($query_na_fasta, \@fasta);
		#unless ($status) { die "\n\t Input error: couldn't open FASTA probe library\n\n"; }
		my $num_fasta = scalar @fasta;
		print "\n\n\t   '$num_fasta' FASTA formatted nucleotide sequences will be used as probes";
		my $i = 0;
		foreach my $seq_ref (@fasta) {
			$i++;
			my $header  = $seq_ref->{header};
			my %header_data;
			my $valid = $self->parse_fasta_header_data($header, \%header_data);
			if ($valid) {
				my $name     = $header_data{name};
				my $gene_name = $header_data{gene_name};
				my $utr_seq   = $seq_ref->{sequence};
				#$devtools->print_hash(\%header_data);
				$self->add_na_probe($probes_ref, $name, $gene_name, $utr_seq);
			}
		}
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
		$valid = undef;
	}
	
	$data_ref->{name}     = $name;
	$data_ref->{gene_name} = $gene_name;

	# DEBUG	
	#print "\n\t # HEADER $header";	
	#print "\n\t # NAME   $name";	
	#print "\n\t # ORF    $gene_name";	

	return $valid;
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
# CREATE THE LIST OF PAIRED BLAST QUERIES FOR SCREENING
############################################################################

#***************************************************************************
# Subroutine:  set_queries
# Description: set up the individual BLAST searches 
#***************************************************************************
sub set_queries {
	
	my ($self, $pipeline_obj, $probes_ref, $targets_ref, $queries_ref) = @_;

	print "\n\n\t ### Creating the BLAST queries\n";
	my $db = $pipeline_obj->{db};
	unless ($db) { die; }
	
	# Get data from self
	my $report_dir = $self->{report_dir};
	my $tmp_path   = $self->{tmp_path};
	unless ($report_dir and $tmp_path) { die; } 
	#$devtools->print_hash($targets_ref); #die; # DEBUG

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

	my $i;
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

			# Else store the query
			#print "\n\t\t #~#~# Setting query: '$probe_id' vs '$target_name'";
			$i++;
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

	unless ($i) {
		print "\n\n\t ### No screening queries were loaded\n";
	}
	else {
		print "\n\t\t $i paired BLAST queries were loaded\n";
	}

	#$devtools->print_hash($queries_ref); die; # DEBUG
	return $i;
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

	# Initialise target sequence library
	my $genome_obj = TargetDB->new($self); 
	$genome_obj->refresh_genomes($targets_ref);
	
	# Iterate through targets set target file paths
	my %paths;
	my %target_data;
	my $genome_use_path  = $self->{genome_use_path};
	my $target_paths_ref = $self->{target_paths};
	unless ($target_paths_ref) { die; } 
	print "\n\n\t ### Getting the target sequences to screen\n";
	foreach my $path (@$target_paths_ref) {
		
		my $full_path = $genome_use_path . "/$path";	
		#print "\n\t Opening $full_path";
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
	my @keys = keys %$targets_ref;
	my $unique_targets = scalar @keys;
	unless ($unique_targets) { die "\n\t No target databases found\n\n\n";	}
	
	print "\n\n\t '$unique_targets' FASTA formatted target files identified in target directory";
	#$devtools->print_hash($targets_ref); die; # DEBUG
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
		#if ($file_type eq 'fa' or $file_type eq 'fas' or $file_type eq 'fasta') {
		if ($file_type eq 'fa') {
			
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
			print "\n\t Target: $target_id";

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

############################################################################
# PARSING ETC
############################################################################

#***************************************************************************
# Subroutine:  create output directories
# Description: create a unique 'report' directory for this process
#***************************************************************************
sub create_output_directories {
	
	my ($self, $pipeline_obj) = @_;

	# Get process ID to create unique output directory
	my $process_id   = $self->{process_id};
	print "\n\n\t ### Screening process ID is '$process_id'\n";
	
	# Create a unique ID and report directory for this run
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
# Subroutine:  parse control file
# Description: read an input file to get parameters for screening
#***************************************************************************
sub parse_control_file {

	my ($self, $ctl_file, $pipeline_obj) = @_;
	
	# Read input file
	my @ctl_file;
	my $valid = $fileio->read_file($ctl_file, \@ctl_file);

	# Parse the 'SCREENDB' block
	$self->parse_screendb_block(\@ctl_file);

	# Parse the 'SCREENSETS' block
	$self->parse_screensets_block(\@ctl_file);

	# READ the 'TARGETS' block
	$self->parse_targets_block(\@ctl_file);
	
	# READ the 'SCREENSQL' block
	$self->parse_screensql_block(\@ctl_file);

	# READ the 'CONSOLIDATION' block
	$self->parse_consolidation_block(\@ctl_file);
	
    # READ the 'UCSC' block
    $self->parse_UCSC_block(\@ctl_file);
	
	# Transfer parameters from this object to the pipeline object
	$pipeline_obj->{mysql_server}           = $self->{mysql_server};
	$pipeline_obj->{mysql_username}         = $self->{mysql_username};
	$pipeline_obj->{mysql_password}         = $self->{mysql_password};
	$pipeline_obj->{db_name}                = $self->{db_name};
	$pipeline_obj->{server}                 = $self->{mysql_server};
	$pipeline_obj->{password}               = $self->{mysql_password};
	$pipeline_obj->{username}               = $self->{mysql_username};
	$pipeline_obj->{output_path}            = $self->{output_path};
	$pipeline_obj->{blast_orf_lib_path}     = $self->{blast_orf_lib_path};
	$pipeline_obj->{blast_utr_lib_path}     = $self->{blast_utr_lib_path};
	$pipeline_obj->{seq_length_minimum}     = $self->{seq_length_minimum};
	$pipeline_obj->{seq_length_minimum}     = $self->{seq_length_minimum};
	$pipeline_obj->{bit_score_min_tblastn}  = $self->{bit_score_min_tblastn};
	$pipeline_obj->{bit_score_min_blastn}   = $self->{bit_score_min_blastn};
	$pipeline_obj->{blast_orf_lib_path}     = $self->{blast_orf_lib_path};
	$pipeline_obj->{blast_utr_lib_path}     = $self->{blast_utr_lib_path};
	$pipeline_obj->{redundancy_mode}        = $self->{redundancy_mode};
	$pipeline_obj->{threadhit_probe_buffer} = $self->{threadhit_probe_buffer};
	$pipeline_obj->{threadhit_gap_buffer}   = $self->{threadhit_gap_buffer};
	$pipeline_obj->{threadhit_max_gap}      = $self->{threadhit_max_gap};
	$pipeline_obj->{reference_glue}         = $self->{reference_glue};
	$pipeline_obj->{query_glue}             = $self->{query_glue};
	$pipeline_obj->{target_paths}           = $self->{target_paths};
	$pipeline_obj->{length_threshold}       = $self->{length_threshold_between_ORFs};
	$pipeline_obj->{genome_structure}       = $self->{genome_structure};
	$pipeline_obj->{ucsc_extracted}			= $self->{UCSC_extracted};
	$pipeline_obj->{ucsc_loci}				= $self->{UCSC_loci};
	# Screen SQL
	my $select_list     = $self->{select_list};
	my $where_statement = $self->{where_statement};
	if ($select_list) {
		#print "\n\t Loaded Select list '$select_list'";
		$self->{select_list}  = lc $select_list;
	}
	if ($where_statement) {
		#print "\n\t WHERE statement '$where_statement'";
		$self->{where_statement} = lc $where_statement;
	}

}

#***************************************************************************
# Subroutine:  parse_screendb_block
# Description: read an input file to get parameters for screening
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
	my $tblastn_min        = $self->{bit_score_min_tblastn};
	my $blastn_min         = $self->{bit_score_min_blastn};
	my $query_aa_fasta     = $self->{query_aa_fasta};
	my $reference_aa_fasta = $self->{reference_aa_fasta};
	my $query_na_fasta     = $self->{query_na_fasta};
	my $reference_na_fasta = $self->{reference_na_fasta};
	my $query_glue         = $self->{query_glue};
	my $reference_glue     = $self->{reference_glue};
	my $redundancy_mode    = $self->{redundancy_mode};
	my $threadhit_probe_buffer = $self->{threadhit_probe_buffer};
	my $threadhit_gap_buffer   = $self->{threadhit_gap_buffer};
	my $threadhit_max_gap      = $self->{threadhit_max_gap};
	my $output_path      = $self->{output_path};

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
		# Validate query FASTA sequence data
		$self->validate_reflib_fasta('aa', $reference_aa_fasta);
	}
	# Validation for a nucleic acid, FASTA-based screen
	if ($query_na_fasta) {
		unless ($tblastn_min) { # Set to default minimum
			$self->{bit_score_min_tblastn} = $default_tblastn_min;
		}
		unless ($reference_na_fasta) { # Set to default minimum
		  die "\n\t Control file error: no NT reference library defined for NT query set\n\n\n";
		}
		# Validate query FASTA sequence data
		$self->validate_reflib_fasta('na', $reference_aa_fasta);
	}
	# Validation for a GLUE-based screen
	if ($query_glue) {
		# TODO
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
		if ($line =~ /^\s*#/)   { next; } # discard comment line 
		chomp $line;
		push (@targets, $line);
		$targets++;
	}
	$self->{target_paths} = \@targets;
}

#***************************************************************************
# Subroutine:  parse_consolidation_block
# Description: read an input file to get parameters for screening
#***************************************************************************
sub parse_consolidation_block {

	my ($self, $file_ref) = @_;

	# Parse the 'CONSOLIDATION' block
	my $start = 'BEGIN CONSOLIDATION';
	my $stop  = 'ENDBLOCK';
	my $db_block = $fileio->read_standard_field_value_block($file_ref, $start, $stop, $self);
	unless ($db_block)  {
		print "\n\t Warning no 'CONSOLIDATION' block found\n\n\n";
		return;
	}

	# Get the 'CONSOLIDATION' block values and validate
	my $length_threshold    = $self->{length_threshold_between_ORFs};
	my $genome_structure    = $self->{genome_structure};

	unless ($length_threshold)  {
		die "\n\t Control file error: 'length_threshold_between_ORFs' undefined in 'CONSOLIDATION' block\n\n\n";
	}
	unless ($genome_structure)  {
		die "\n\t Control file error: 'genome_structure' undefined in 'CONSOLIDATION' block\n\n\n";
	}
}

#***************************************************************************
# Subroutine:  parse_UCSC_block
# Description: read an input file to get parameters for screening
#***************************************************************************
sub parse_UCSC_block {

    my ($self, $file_ref) = @_;
    # Parse the 'CONSOLIDATION' block
    my $start = 'BEGIN UCSC';
    my $stop  = 'ENDBLOCK';
    my $db_block = $fileio->read_standard_field_value_block($file_ref, $start, $stop, $self);
    unless ($db_block)  {
        print "\n\t Warning no 'UCSC' block found\n\n\n";
        return;
    }    
    # Get the 'CONSOLIDATION' block values and validate
    my $ucsc_extracted    	= $self->{UCSC_extracted};
    my $ucsc_loci		    = $self->{UCSC_loci};
    unless ($ucsc_extracted or $ucsc_loci)  {
        die "\n\t Control file error: 'UCSC table' undefined in 'UCSC' block\n\n\n";
    }    
}


#***************************************************************************
# Subroutine:  parse_screensql_block
# Description: read an input file to get parameters for screening
#***************************************************************************
sub parse_screensql_block {

	my ($self, $file_ref) = @_;

	# READ the 'SCREENSQL' block
	my $start = 'BEGIN SCREENSQL';
	my $stop  = 'ENDBLOCK';
	my $sql_block = $fileio->read_sql_block($file_ref, $start, $stop, $self);
	#$devtools->print_hash($self); die;
}


#***************************************************************************
# Subroutine:  validate_reflib_fasta
# Description: TODO 
#***************************************************************************
sub validate_reflib_fasta {

	my ($self, $seq_type, $file_path) = @_;



}

############################################################################
# EOF
############################################################################
