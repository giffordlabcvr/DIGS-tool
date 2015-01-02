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
	#$devtools->print_hash($self); die; 

	# Do dialogue for selecting sequences
	$self->do_retrieve_data_dialogue();
	
	# Do dialogue for selecting sequences
	$self->do_align_data_dialogue();

	# Analyse GLUE alignment
	#$glue_obj->{glue_msa_path};
	#$glue_obj->analyse_glue_msa();

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
	$glue_obj->{refseq_path};

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

############################################################################
# EOF
############################################################################
