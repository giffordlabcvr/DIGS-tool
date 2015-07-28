#!/usr/bin/perl -w
############################################################################
# Module:      Retrieve.pm
# Description: Genome screening pipeline using reciprocal BLAST
# History:     December 2009: Created by Robert Gifford 
############################################################################
package Retrieve;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::FileIO;
use Base::DevTools;
use Base::Console;
use Base::Sequence;    # For performing basic sequence manipulations

# Program components
use DIGS::ScreenBuild; # Functions to set up screen

############################################################################
# Globals
############################################################################

# Default minimum length of BLAST hit to extract 
my $default_min_seqlen = 100; # minimum sequence length of BLAST hit to extract
my $process_dir_warn   = 100; # folders in process directory to trigger warning 

# Base objects
my $fileio    = FileIO->new();
my $devtools  = DevTools->new();
my $console   = Console->new();
my $seq_obj   = Sequence->new();
my $verbose   = undef;
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create new Retrieve 'object'
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {
		
		# Flags
		process_id             => $parameter_ref->{process_id},
		
		# Paths and member variables
		blast_bin_path         => $parameter_ref->{blast_bin_path},
		genome_use_path        => $parameter_ref->{genome_use_path},
		output_path            => $parameter_ref->{output_path},

		# Database variables
		db                     => $parameter_ref->{db},
	
		# Member classes 
		blast_obj              => $parameter_ref->{blast_obj},
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# SECTION: Handler functions 
############################################################################

#***************************************************************************
# Subroutine:  run_data_retrieval_functions
# Description: handler for data retrieval in the DIGS framework
#***************************************************************************
sub run_data_retrieval_functions {

	my ($self, $ctl_file) = @_;

	# Get information from screening db
	$self->initialise_screening_db_ranges();
	$devtools->print_hash($self); die; 
	
	die;

	# Do dialogue for selecting sequences
	$self->do_retrieve_data_dialogue();
	
	# Do dialogue for selecting sequences
	$self->do_align_data_dialogue();


	# Analyse GLUE alignment
	#$glue_obj->{glue_msa_path};
	#$glue_obj->analyse_glue_msa();
	
	# Initialise usage statement to print if usage is incorrect
	#my($USAGE) = "\n\tUsage:\n";
	#  $USAGE  .= "\n\tSingle extract: ";
	#  $USAGE  .= "\n\t\t$0 [start] [stop] [scaffold] [path-to-file]";
	#  $USAGE  .= "\n\tExtract multiple from file:";
	#  $USAGE  .= "\n\t\t$0 [data file] -f1";
	#  $USAGE  .= "\n\tExtract multiple from file and extend:";
	#  $USAGE  .= "\n\t\t$0 [data file] [upstream] [downstream] -f2";
	#  $USAGE  .= "\n\tExtract using SQL query, and extend:";
	#  $USAGE  .= "\n\t\t$0 [data file] [upstream] [downstream] -f3";
	#  $USAGE  .= "\n\n";

	# Viable program options for command line executions 
	my %allowed_options;
	$allowed_options{f} = 1;
	$allowed_options{1} = 1;
	$allowed_options{2} = 1;
	$allowed_options{3} = 1;

	############################################################################
	# Begin code
	############################################################################

	# Collect command options if received and set run parameters accordingly
	my $options = pop @ARGV;
	my @options;
	if ($options) {
		@options = split('', $options);
	}
	my %options;
	foreach my $option (@options) {
		if ($option eq '-') { next; }
		#unless ($allowed_options{$option}) { die $USAGE; };
		$options{$option} = 1;
	}

	# Show the about box
	$console->refresh();
	show_title();

	if ($options{f}) {
		
		my $file_path = $ARGV[0];
		#unless ($file_path) { die $USAGE; }

		if ($options{1}) {
			extract_multiple_sequences($file_path);
		}
		
		elsif ($options{2}) {
			my $file_path  = $ARGV[0];
			my $upstream   = $ARGV[1];
			my $downstream = $ARGV[2];
			unless ($file_path and $upstream and $downstream) {
		#		die $USAGE; 
			}
			extend_and_extract_multiple_sequences($file_path, $upstream, $downstream);
		}
		elsif ($options{3}) {
			my $file_path  = $ARGV[0];
			my $upstream   = $ARGV[1];
			my $downstream = $ARGV[2];
			unless ($file_path and $upstream and $downstream) {
		#		die $USAGE; 
			}
			extract_using_sql_query($file_path, $upstream, $downstream);
		}
	}

	else {
		my $start    = $ARGV[0];
		my $stop     = $ARGV[1];
		my $scaffold = $ARGV[2];
		my $chunk    = $ARGV[3];
		unless ($start and $stop and $scaffold and $chunk) {
		#	die $USAGE; 
		}
		extract($start, $stop, $scaffold, $chunk);
	}

	# Exit program
	print "\n\n\t...process completed\n\n";
	exit;
}

#***************************************************************************
# Subroutine:  do_retrieve_data_dialogue
# Description: 
#***************************************************************************
sub do_retrieve_data_dialogue {

	my ($self, $ctl_file) = @_;

	my $db = $self->{db};

	# Get data from the Extracted table or the Loci table
	my @choices = qw [ e l ];
	my $question = "\n\t Use the Extracted table or the Loci table?";
	my $choice = $console->ask_simple_choice_question($question, \@choices);	
	my $table; 
	if ($choice eq 'e') {     # Extracted
		$table = $db->{extracted_table};
	}
	elsif ($choice eq 'l') {  # Loci
		$table = $db->{loci_table};
		die;
	}

	# Select an approach for retrieving a set of homologous sequences
	# By one or more of:
	my $genes_ref = $self->{assigned_genes};
	my %gene_names;
	my $i = 0;
	my @where; 
	foreach my $gene_ref (@$genes_ref) {
		$i++;
		my $gene_name = $gene_ref->{assigned_gene};
		print "\n\t $i: $gene_name";
		$gene_names{$i} = $gene_name;
	}
	$question = "\n\n\t Which gene to use? ";
	$choice = $console->ask_list_question($question, $i);
	my $assigned_gene = $gene_names{$choice};
	push (@where, "WHERE  assigned_gene = '$assigned_gene' ");

	# Get the reference sequences
	my $names_ref     = $self->{assigned_names};
	$i = 0;
	my %refseq_names;
	foreach my $name_ref (@$names_ref) {
		$i++;
		my $assigned_name = $name_ref->{assigned_name};
		print "\n\t $i: $assigned_name";
		$refseq_names{$i} = $assigned_name;
	}
	$question = "\n\n\t Which refseq to use? ";
	$choice = $console->ask_list_question($question, $i);
	my $assigned_name = $refseq_names{$choice};
	push (@where, "AND assigned_name = '$assigned_name' ");

	# Get the organisms
	$question = "\n\n\t Constrain by host? ";
	my %organisms;
	my $constrain = $console->ask_yes_no_question($question);
	if ($constrain eq 'y') {
		my $organisms_ref = $self->{organisms};
		$i = 0;
		foreach my $organism_ref (@$organisms_ref) {
			$i++;
			my $organism = $organism_ref->{organism};
			print "\n\t $i: $organism";
		}
		$question = "\n\n\t Which organism to use? ";
		$choice = $console->ask_list_question($question, $i);
		my $organism = $organisms{$choice};	
		push (@where, "AND organism = '$organism' ");
	}

	# Get the organisms
	$question = "\n\n\t Set the bitscore cutoff ";
	my $low  = $self->{bitscore_low};
	my $high = $self->{bitscore_high};
	unless ($low and $high) { die; }
	my $cutoff = $console->ask_int_with_bounds_question($question, $low, $high);	
	push (@where, " AND bit_score > $cutoff ");
	my $where = join ("\n", @where);
	
	# Create SQL to get sequences
	my @fields = qw [ record_id assigned_name assigned_gene 
                      organism target_name scaffold orientation
                      extract_start extract_end
                      sequence ];
	my @data;
	$table->select_rows(\@fields, \@data, $where);

	my @fasta_sequences;
	my @sequence_data;
	$self->{fasta_sequences} = \@fasta_sequences;
	$self->{sequence_data}   = \@sequence_data;

}

#***************************************************************************
# Subroutine:  do_align_data_dialogue
# Description: 
#***************************************************************************
sub do_align_data_dialogue {

	my ($self, $ctl_file) = @_;

	# Do you want to align the sequences?
	my $fasta_sequences = $self->{fasta_sequences};
	my $sequence_data   = $self->{sequence_data};

	# Instantiate main program classes using global settings
	my $glue_obj = GLUE->new($self);
	#$glue_obj->{refseq_path};

	$glue_obj->{query_fasta} = $fasta_sequences;
	$glue_obj->{data_path} = $sequence_data;
	die;

	# Set up the GLUE process	
	my $builder  = GLUE_Build->new();
	$builder->setup_glue_process($glue_obj);
		
	# Create GLUE alignment
	my $success = $glue_obj->create_glue_alignment();
	unless ($success) {
		print "\n\t No sequences were aligned\n\n"; 
		exit;
	}

}

############################################################################
# SECTION: Retrieve data from screening database
############################################################################

#***************************************************************************
# Subroutine:  retrieve 
# Description: retrieve data from screening database
#***************************************************************************
sub retrieve {

	my ($self, $ctl_file) = @_;

	# Get params from self	
	my $select = $self->{select_list};
	my $where  = $self->{where_statement};
	unless ($select) { die; }

	my $db = $self->{db};
	my @select;
	my @tmp  = split(/,/, $select);
	foreach my $field (@tmp) {
		$field =~ s/\s+//;
		push (@select, $field);
	}
	#$devtools->print_array(\@select); die;

	# Get params from self	
	my @sequences;
	$db->retrieve_sequences(\@sequences, \@select, $where);
	my $seqfile = 'sequences.fa';
	$fileio->write_file($seqfile, \@sequences);

}


#***************************************************************************
# Subroutine:  initialise_screening_db_ranges
# Description: 
#***************************************************************************
sub initialise_screening_db_ranges {

	my ($self, $ctl_file) = @_;

	# get database objects
	my $db = $self->{db};
	unless ($db) { die; }
	
	#  get tables
	my $extracted_table = $db->{extracted_table};
	my $loci_table      = $db->{loci_table};
	
	#  reference sequence
	my @assigned_names;
	my @fields1 = qw [ assigned_name ];
	$extracted_table->select_distinct(\@fields1, \@assigned_names);
	$self->{assigned_names} = \@assigned_names;

	#  gene
	my @assigned_genes;
	my @fields2 = qw [ assigned_gene ];
	$extracted_table->select_distinct(\@fields2, \@assigned_genes);
	$self->{assigned_genes} = \@assigned_genes;
	
	#  organism
	my @organisms;
	my @fields3 = qw [ organism ];
	$extracted_table->select_distinct(\@fields3, \@organisms);
	$self->{organisms} = \@organisms;

	#  bitscore
	my @fields4 = qw [ bit_score ];
	my @bitscores;
	my $where = " ORDER BY bit_score DESC ";
	$extracted_table->select_rows(\@fields4, \@bitscores, $where);
	my $data_ref = pop @bitscores;
	my $bitscore_low  = $data_ref->{bit_score};
	$data_ref = shift @bitscores;
	my $bitscore_high  = $data_ref->{bit_score};
	#print "\n\t low $bitscore_low high $bitscore_high";
	$self->{bitscore_low}  = $bitscore_low;
	$self->{bitscore_high} = $bitscore_high;

}


#***************************************************************************
# Subroutine:  extract
# Description: extract a single sequence  
#***************************************************************************
sub extract {
		
	my ($start, $stop, $scaffold, $chunk_path) = @_;
	
	my $extract_obj = Extract->new();
	my $sequence = $extract_obj->extract_single_match($start, $stop, $scaffold, $chunk_path);
	unless ($sequence) {
		print "\n\t scaffold $scaffold : $start : $stop\n\n";
		print "\n\t path $chunk_path";
		print "\n\t NOT FOUND!!";
		die;
	}
	my $fasta = ">$scaffold $start-$stop\n$sequence\n\n";
	$fileio->write_text_to_file('extracted.txt', $fasta);
}

#***************************************************************************
# Subroutine:  extract_multiple_sequences
# Description: get a set of sequences using coordinates in a file
#***************************************************************************
sub extract_multiple_sequences {

	my ($file_path) = @_;
	
	my %genome_data;
	index_genome_data(\%genome_data);
	#$devtools->print_hash(\%genome_data); # DEBUG

	my @file_data;
	$fileio->read_input_file($file_path, \@file_data);
	my $i;
	my $output_path = $file_path . '.extracted.txt';
	foreach my $line (@file_data) {
		
		print $line;
		chomp $line;
		
		$i++; # Increment line count
		if ($line =~ /^\s*$/) { next; } # discard blank line

		# Get data from the line
		my @line = split("\t", $line);
		my $chunk       = $line[0];
		my $scaffold    = $line[1];
		my $start       = $line[2];
		my $stop        = $line[3];
		my $orientation = $line[4];
		unless ($start and $stop and $scaffold and $chunk) {
			die "\n\t Input data error: missing parameters at line $i\n\n\n"; 
		}
		
		my $chunk_path = $genome_data{$chunk};
		unless ($chunk_path) {
			die "\n\t No path found for chunk '$chunk'  $i\n\n\n"; 
		}

		my $extract_obj = Extract->new();
		my $sequence = $extract_obj->extract_single_match($start, $stop, $scaffold, $chunk_path);
		unless ($sequence) {
			print "\n\t scaffold $scaffold : $start : $stop\n\n";
			print "\n\t path $chunk_path";
			print "\n\t NOT FOUND!!";
			die;
		}
		if ($orientation eq '-ve') {
			$sequence = $seq_obj->reverse_and_complement($sequence);
		}
		my $fasta = ">$chunk $scaffold $start-$stop\n$sequence\n\n";
		$fileio->append_text_to_file($output_path, $fasta);
	}
	
}

#***************************************************************************
# Subroutine:  extend_and_extract_multiple_sequences
# Description: get a set of sequences using coordinates in a file,
#  extending by the amounts specified in $upstream and $downstream
#***************************************************************************
sub extend_and_extract_multiple_sequences {

	my ($file_path, $upstream, $downstream) = @_;
	
	my %genome_data;
	index_genome_data(\%genome_data);
	#$devtools->print_hash(\%genome_data); # DEBUG

	my @file_data;
	$fileio->read_input_file($file_path, \@file_data);
	my $i;
	my $output_path = $file_path . '.extracted.txt';
	foreach my $line (@file_data) {
		
		print $line;
		chomp $line;
		
		$i++; # Increment line count
		if ($line =~ /^\s*$/) { next; } # discard blank line

		# Get data from the line
		my @line = split("\t", $line);
		my $chunk       = $line[0];
		my $scaffold    = $line[1];
		my $start       = $line[2];
		my $stop        = $line[3];
		my $orientation = $line[4];
		unless ($start and $stop and $scaffold and $chunk) {
			die "\n\t Input data error: missing parameters at line $i\n\n\n"; 
		}
		
		my $chunk_path = $genome_data{$chunk};
		unless ($chunk_path) {
			die "\n\t No path found for chunk '$chunk'  $i\n\n\n"; 
		}

		# Adjust coordinates
		my $new_start = $start - $upstream;
		if ($new_start < 0) { $new_start = 1; }
		my $new_stop  = $stop + $downstream;

		my $extract_obj = Extract->new();
		my $sequence = $extract_obj->extract_single_match($new_start, $new_stop, 
		                                                $scaffold, $chunk_path);
		unless ($sequence) {
			print "\n\t scaffold $scaffold : $start : $stop\n\n";
			print "\n\t path $chunk_path";
			print "\n\t NOT FOUND!!";
			die;
		}
		
		if ($orientation eq '-ve') {
			$sequence = $seq_obj->reverse_and_complement($sequence);
		}
		my $fasta = ">$chunk $scaffold $start-$stop\n$sequence\n\n";
		#print "\n$fasta\n\n\n";
		$fileio->append_text_to_file($output_path, $fasta);
	}
	
}

#***************************************************************************
# Subroutine:  extract_using_sql_query
# Description: 
#***************************************************************************
sub extract_using_sql_query {
	
	my ($query_file, $upstream, $downstream) = @_;

	my @file;
	$fileio->read_input_file($query_file, \@file);
	my $query = join(" ", @file);
	my @data;
	sql_retrieve($query, \@data);
	my  $i = 0;
	my $sql_outfile = $query_file . '.out.txt';
	$fileio->write_output_file($sql_outfile, \@data);

	my $sequence_obj = Sequence->new();
	my %genome_data;
	index_genome_data(\%genome_data);
	#$devtools->print_hash(\%genome_data); # DEBUG
	
	# Iterate through sql query results, doing extractions
	my $outfile = $query_file . '.extracted.txt';
    foreach my $line (@data) {
		
		$i++;
		my @line = split("\t", $line);
		my $chunk    = $line[0];
		my $scaffold = $line[1];
		my $start    = $line[2];
		my $stop     = $line[3];
		my $orientation  = $line[4];
		print "DOING $i $chunk\n";
		
		my $chunk_path = $genome_data{$chunk};
		unless ($chunk_path) {
			die "\n\t No path found for chunk '$chunk'  $i\n\n\n"; 
		}
		
		
		#if ($orientation eq '-ve') {
		#	$start = $start + 30;
		#	$stop = $stop - 9;
		#	if ($stop < 1) { $stop = 1; }
		#}
		#else {
		$start = $start - $upstream;
		if ($start < 1) { $start = 1; }
		$stop = $stop + $downstream;
		#}
		
		my $extract_obj = Extract->new();
		my $sequence = $extract_obj->extract_single_match($start, $stop, $scaffold, $chunk_path);
		unless ($sequence) {
			print "\n\t scaffold $scaffold : $start : $stop\n\n";
			print "\n\t path $chunk_path";
			print "\n\t NOT FOUND!!";
			die;
		}
		if ($orientation eq '-ve') {
			$sequence = $sequence_obj->reverse_and_complement($sequence);
		}
		my $fasta = ">$chunk $scaffold $start-$stop\n$sequence\n\n";
		$fileio->append_text_to_file($outfile, $fasta);
							
	}
	#$devtools->print_hash(\%data);
}

#***************************************************************************
# Subroutine:  sql_retrieve
# Description: generic fxn to execute SQL query and retrieve tab delimited 
#              results
#***************************************************************************
sub sql_retrieve {

	my ($query, $data_ref) = @_;

	# Execute the query
    print "\n##################################";
    print "\n\n$query\n";	
    print "\n##################################\n\n";
    #my $sth = $dbh->prepare($query);
    #unless ($sth->execute()) { print $query; exit;}	
    #my @data;
	#while (my $row = $sth->fetchrow_arrayref) {

     #   my @row = @$row;
      #  my $data = join("\t", @row);
	#	#print "\n\t$data\n";
	#	push (@$data_ref, "$data\n");
    #}
}

#***************************************************************************
# Subroutine:  index_genome_data 
# Description:  
#***************************************************************************
sub index_genome_data {

	my ($genomes_ref) = @_;
	
	die;
	my $genome_path;
	my @file_data;
	$fileio->read_input_file($genome_path, \@file_data);
	my $i;
	foreach my $line (@file_data) {
		
		#print $line;
		chomp $line;
		$i++; # Increment line count
		if ($line =~ /^\s*$/) { next; } # discard blank line

		# Get data from the line
		my @line = split("\t", $line);
		my $species    = $line[0];
		my $chunk_name = $line[1];
		my $source     = $line[2];
		my $version    = $line[3];
		$chunk_name =~ s/\s+//g;
		$source  =~ s/\s+//g;
		$version =~ s/\s+//g;
		$species =~ s/\s+//g;
		my $path = "/genomes/macrolineage/Mammalia/$species/$source/$version/$chunk_name";
		$genomes_ref->{$chunk_name} = $path;
	}
}


############################################################################
# EOF
############################################################################
