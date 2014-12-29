#!/usr/bin/perl -w
############################################################################
# Script:       GLUE.pm 
# Description:  Create GLUE MSAs
# History:      Rob Gifford, March 2013: Creation
############################################################################
package GLUE;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::IO;
use Base::BioIO;
use Base::SeqIO;
use Base::FileIO;
use Base::DevTools;
use Base::Sequence;

# Component Classes
use GLUE::RefSeqParser;
use GLUE::RefSeq;
use GLUE::RefSeqAlignment;
use GLUE::RefSeqAlignment;
use GLUE::DataTool;
#use Phylogeny::Tree;

############################################################################
# Globals
############################################################################

# Create Base objects
my $io       = IO->new();
my $fileio   = FileIO->new();
my $bioio    = BioIO->new();
my $seqio    = SeqIO->new();
my $devtools = DevTools->new();
my $console  = Console->new();
my $writer   = HTML_Utilities->new();

# For displaying alignment progress
my $increment = 50;

1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create new GLUE program class 
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Store sequence IDs in the order they were submitted
	my @seq_ids;
	my $sequences = $parameter_ref->{sequences};
	foreach my $seq_ref (@$sequences) {
		my $seq_id = $seq_ref->{header};
		push (@seq_ids, $seq_id);
	}

	# Set defaults
	unless ($parameter_ref->{output_type}) {
		$parameter_ref->{output_type} = 'text';
	}

	# Set parameters for the VGlue object
	my $self = {
		
		# Query data set
		sequence_ids         => \@seq_ids,
		sequences            => $parameter_ref->{sequences},
		raw_tab_data         => $parameter_ref->{raw_tab_data},
		glue_msa_obj         => $parameter_ref->{glue_msa_obj},
		set_name             => $parameter_ref->{set_name},
		
		# Paths and constants
		process_id            => $parameter_ref->{process_id},
		refseq_lib_path       => $parameter_ref->{refseq_lib_path},
		refseq_use_path       => $parameter_ref->{refseq_use_path},
		output_path           => $parameter_ref->{output_path},
		tmp_path              => $parameter_ref->{tmp_path},
		report_dir            => $parameter_ref->{report_dir},
		html_header           => $parameter_ref->{html_header},

		# Component interface classes
		muscle_obj            => $parameter_ref->{muscle_obj},
		blast_obj             => $parameter_ref->{blast_obj},
		lap_obj               => $parameter_ref->{lap_obj},
		phylogeny_obj         => $parameter_ref->{phylogeny_obj},
		
		# Flags
		output_type               => $parameter_ref->{output_type},
		override_blast            => $parameter_ref->{override_blast},
		lap_align                 => $parameter_ref->{lap_align},
		muscle_align              => $parameter_ref->{muscle_align},
		phylogeny                 => $parameter_ref->{phylogeny},
		mut_profile               => $parameter_ref->{mut_profile},
		msa_qa_stats              => $parameter_ref->{msa_qa_stats},
		mut_profile               => $parameter_ref->{mut_profile},
		mut_frequencies           => $parameter_ref->{mut_frequencies},
		stratify_field            => $parameter_ref->{stratify_field},
		derive_polymorphism_list  => $parameter_ref->{derive_polymorphism_list},
		polymorphism_threshold    => $parameter_ref->{polymorphism_threshold},
		show_synonymous_changes   => $parameter_ref->{show_synonymous_changes},
		all_mutations             => $parameter_ref->{all_mutations},
		show_mutlist_seqs         => $parameter_ref->{show_mutlist_seqs},
			
		# DB connection variables
		server                => $parameter_ref->{mysql_server},
		username              => $parameter_ref->{mysql_username},
		password              => $parameter_ref->{mysql_password},
	};	
	
	bless ($self, $class);
	return $self;
}

############################################################################
# Create GLUE MSA
############################################################################

#***************************************************************************
# Subroutine:  create_glue_alignment
# Description: pairwise align query sequences to references 
#***************************************************************************
sub create_glue_alignment {
	
	my ($self) = @_;

	# Get settings/paths/data from self	
	my $output_type   = $self->{output_type};
	my $flat_dir      = $self->{refseq_use_path};
	my $refseq        = $self->{refseq};
	my $report_dir    = $self->{report_dir};
	my $html_header   = $self->{html_header};
	my $refseq_name   = $refseq->{name};
	unless ($refseq)      { die "\n\t NO refseq found \n\n"; } # Sanity checking
	unless ($refseq_name) { die "\n\t NO refseq name found\n\n"; } # Sanity checking
	unless ($flat_dir)    { die; } # Sanity checking
	#print "NAME $refseq_name";

	my @html;
	my $indent = ' ' x 10;
	if ($output_type eq 'html') { 
		# Define the page in the correct format and write it
		unless ($html_header) { die "\n\t No HTML header found\n\n"; } # Sanity checking
		my @page_header;
		$fileio->read_file($html_header, \@page_header);
		my $header = join ("\n", @page_header);
		print $header;
	}

	# Set up rolling output
	my $sequences_ref  = $self->{sequences};
	unless ($sequences_ref) {die; }
	my $num_seqs       = scalar @$sequences_ref;
	my $next_i = $self->set_up_rolling_output($num_seqs, $output_type);
	my $i = 0;
	my $start = $i + 1;
	#$devtools->print_array($sequences_ref); die;
	
	# Align each of the input sequences in turn
	my $num_aligned_sequences = 0;
	my $num_failed = 0;
	my @valid;
	my @failed;
	foreach my $seq_ref (@$sequences_ref) {
		
		# Show rolling output
		$i++;
		$next_i = $self->show_rolling_output($start, $next_i, $i, $num_seqs, $output_type);
		my $header = $seq_ref->{header};
		my $id     = $seq_ref->{sequence_id};
		$seq_ref->{raw_seq_length} = length $seq_ref->{sequence};
		
		# Do BLAST align
		$self->glue_blast_align($seq_ref); 
		unless ($seq_ref->{blast_valid}) {
			$num_failed++;
			my $message = " BLAST-align fail $num_failed: sequence '$header'";
			$io->show_output_message($message, $self->{output_type});
			next;
		}
		
		# Do MUSCLE align
		#$self->glue_muscle_align($seq_ref);
		
		# Do LAP ORF align
		#$self->glue_lap_align($seq_ref);
		
		# Reconcile alignments
		$self->reconcile_alignments($seq_ref, \@valid, \@failed); 
	
		if ($self->{mut_profile}) {
			$refseq->{start} = $seq_ref->{aln_start};
			$refseq->{stop}  = $seq_ref->{aln_stop};
			$refseq->compare_to_aligned_sequence($seq_ref);
			#$devtools->print_hash($seq_ref); die;
		}

		# Upload sequence & GLUE MSA data to databae
		if ($self->{upload_to_db}) { # Load data to sequence database
			$self->load_data_to_sequence_db($seq_ref);
		}
	}

	# Store the aligned and failed, unaligned sequences
	$self->{aln_fail_sequences} = \@failed;
	my $num_valid = scalar @valid;
	unless ($num_valid) {
		return 0;	
	}
	
	$self->create_msa_object(\@valid);
	my $glue_msa_obj = $self->{glue_msa_obj};

	# Write the MSA
	my $glue_msa_file = 'aligned_sequences.glu';
	my $glue_msa_path = $report_dir  . $glue_msa_file;
	$glue_msa_obj->write_glue_msa($glue_msa_path);
	my $message ="\n\t # GLUE MSA written to file '$glue_msa_path'";
	$io->show_output_message($message, $self->{output_type});
	$self->{glue_msa_file} = $glue_msa_file; 
	$self->{glue_msa_path} = $glue_msa_path; 
	if ($self->{mut_profile}) {
		$glue_msa_obj->{mutation_profiled} = 'true';
	}

	return 1;	
}

#***************************************************************************
# Subroutine:  analyse_msa
# Description: calculate statistics on GLUE MSA
#***************************************************************************
sub analyse_glue_msa {

	my ($self) = @_;

	# MSA quality analysis processes
	if ($self->{msa_qa_stats}) { 
		$self->calculate_coverage();
		#$self->calculate_shannon_e();
	}
	# MSA mutation frequencies
	if ($self->{mut_frequencies}) { 
		$self->derive_mutation_frequencies();
		if ($self->{stratify_field}) {
			$self->derive_stratified_mutation_frequencies();
		}
		# Create 'typical' mutation list
		if ($self->{derive_polymorphism_list}) { 
			$self->derive_polymorphism_list();
		}	
	}
	# Create phylogeny
	if ($self->{phylogeny}) { 
		$self->create_phylogeny();
	}

	# Write the report
	my $report_obj = GLUE_Report->new($self);
	$report_obj->write_msa_report();
}

#***************************************************************************
# Subroutine:  load_data_to_sequence_db
# Description: add sequence data to DB
#***************************************************************************
sub load_data_to_sequence_db {

	my ($self, $seq_ref) = @_;

	# Get the table references
	my $db_obj = $self->{sequence_db};
	my $mutation_table = $db_obj->{mutation_table};
	my $genotype_table = $db_obj->{genotype_table};
	my $sequence_table = $db_obj->{sequence_table};
	my $location_table = $db_obj->{location_table};
	unless ($sequence_table and $location_table) { 
		die "\n\t Sequence DB has not been correctly loaded\n\n\n";
	}	

	my %data;
	$data{'sequence_id'}     = $seq_ref->{header};
	my $data_ref = $seq_ref->{data};
	my @fields = keys %$data_ref;
	foreach my $field (@fields) {

		my $value = $data_ref->{$field};
		if ($field eq 'isolation_country') {
			$data{'location_id'} = $value;
			$data{'location_type'} = 'country';
		}
		elsif ($field eq 'sequence_id') {
			$data{header} = $value;
		}
		else {
			$data{$field} = $value;
		}
	}
	#$devtools->print_hash(\%data); 

	# Get the sequence length
	my $sequence = $seq_ref->{sequence};
	my $seq_len  = length $sequence;
	unless ($seq_len) { die; }

	# Data
	unless ($data{isolation_date}) {
		$data{isolation_date} = 'NULL';
	}
	$data{'sample_id'} = 'NULL';
	if ($data{'isolate'}) {
		my $sample_id = $data{'isolate'};
		$sample_id =~ s/\s+//g;
		$data{'sample_id'} = $sample_id;
	}
	$data{'sequence'}        = $sequence;
	#$data{'sequence_date'}   = $seq_ref->{sequence_date};
	$data{'sequence_length'} = $seq_len;
	$data{'isolation_date_class'} = 'Genbank';
	$data{'genotype_method'} = 'Genbank';
	if ($data{'genotype'} eq 'NULL') {
		if ($data{'serotype'}) {
			unless ($data{'serotype'} eq 'NULL') { 
				$data{'genotype'} = $data{'serotype'};
				$data{'genotype_method'} = 'serology';
			}
		}
	}
	unless ($data{'sample_id'}) {
		my $isolate = $data{'isolate'};
		$data{'sample_id'} = $isolate;
	}
	unless ($data{'location_id'}) {
		$data{'location_id'} = 'NK';
	}
	$data{'score'}           = 'NA';
	$data{'start'}           = $seq_ref->{aln_start};
	$data{'stop'}            = $seq_ref->{aln_stop};
	#$devtools->print_hash(\%data); # die;
	$sequence_table->insert_row(\%data);
	$location_table->insert_row(\%data);
	$genotype_table->insert_row(\%data);

	#$devtools->print_hash($seq_ref); die;
	$self->load_mutation_data($seq_ref);
	#$self->load_insertion_data($seq_ref);
	#$devtools->print_hash(\%fields); die;
	#$devtools->print_array(\@parsed); die;
}


#***************************************************************************
# Subroutine:  load_mutation_data
# Description: add sequence data to DB
#***************************************************************************
sub load_mutation_data {

	my ($self, $seq_ref) = @_;

	# Get the table references
	my $db_obj = $self->{sequence_db};
	my $mutation_table = $db_obj->{mutation_table};
	unless ($mutation_table) {
		die "\n\t Sequence DB has not been correctly loaded\n\n\n";
	}	

	my $mutations_ref = $seq_ref->{nonsyn_mutations};
	my %data;
	foreach my $mutation_ref (@$mutations_ref) {
		#$devtools->print_hash($mutation_ref); die;
		$data{sequence_id} = $seq_ref->{header};
		$data{gene_id}     = $mutation_ref->{gene};
		$data{level}       = '1.0';
		$data{reference}   = $mutation_ref->{refseq_aa};
		$data{position}    = $mutation_ref->{position};
		$data{mutation}    = $mutation_ref->{aa};
		#$devtools->print_hash(\%data); die; 
		$mutation_table->insert_row(\%data);
	}
}

############################################################################
# MSA QA stats 
############################################################################

#***************************************************************************
# Subroutine:  calculate_coverage 
# Description: calculate MSA coverage 
#***************************************************************************
sub calculate_coverage {

	my ($self) = @_;
	
	# Get the GLUE alignment object
	my $alignment = $self->{glue_msa_obj};
	unless ($alignment) { die; }  # Sanity checking
	
	# Calculating coverage in the alignment 
	my $message = "\n\t # Calculating MSA coverage";
	$io->show_output_message($message, $self->{output_type}); 
	$alignment->calculate_coverage();
}

#***************************************************************************
# Subroutine:  calculate_shannon_entropy 
# Description: calculate shannon entropy 
#***************************************************************************
sub calculate_shannon_e {

	my ($self) = @_;

	# Get the GLUE alignment object
	my $alignment = $self->{glue_msa_obj};
	unless ($alignment) { die; }  # Sanity checking

	# Calculating Shannon entropy in the alignment 
	my $message = "\n # Calculating Shannon entropy";
	$io->show_output_message($message, $self->{output_type}); 
	my $align_file = 'MSA.fas';
	my $output_path = $self->{report_dir};
	my $align_path = $output_path . $align_file;	
	$alignment->write_glue_alignment($align_path, 'no insertions');
	print "\n\t # Calulating Shannon entropy for GLUE MSA '$align_file'";
	my @fasta;
	$fileio->read_file($align_path, \@fasta);
	#$alignment->derive_shannon_entropy(\@fasta);
	my $command = "rm $align_path";
	system $command;
}

############################################################################
# Creating MSA mutation profiles and deriving mutation frequencies
############################################################################

#***************************************************************************
# Subroutine:  derive mutation profile
# Description: 
#***************************************************************************
sub derive_mutation_profile {

	my ($self) = @_;

	# Get the GLUE alignment object
	my $alignment = $self->{glue_msa_obj};
	unless ($alignment) { die "\n\t No alignment found\n\n"; }  # Sanity checking

	# Do mutation profiling
	my $message = "\n\t # Profiling sequence mutations";
	$io->show_output_message($message, $self->{output_type}); 
	$alignment->profile_sequences();  # Get mutations in each sequence
	
}
	
#***************************************************************************
# Subroutine:  derive mutation frequencies
# Description: 
#***************************************************************************
sub derive_mutation_frequencies {

	my ($self) = @_;

	# Get the GLUE alignment object
	my $alignment = $self->{glue_msa_obj};
	unless ($alignment) { die; }  # Sanity checking

	# Do mutation profiling
	unless ($alignment->{mutation_profiled}) {
		my $message = "\n\t # Profiling sequence mutations";
		$io->show_output_message($message, $self->{output_type}); 
		$alignment->profile_sequences();  # Get mutations in each sequence
	}
	#$alignment->describe(); die;
	
	# Do mutation frequency calculation 
	my $message = "\n\t # Calculating mutation frequencies for entire data set";
	$io->show_output_message($message, $self->{output_type}); 
	$alignment->get_gene_counts();    # Get gene counts
	$alignment->derive_codon_counts();
	$alignment->derive_mutation_counts();
	$alignment->derive_mutation_frequencies();
	$alignment->create_mutation_frequency_table();
	
}

#***************************************************************************
# Subroutine:  derive stratified mutation frequencies
# Description: Derive stratified mutation frequencies by a given field
#***************************************************************************
sub derive_stratified_mutation_frequencies {

	my ($self) = @_;

	my $alignment = $self->{glue_msa_obj};
	unless ($alignment) { die; }  # Sanity check

	# Do mutation frequency calculation 
	my $field = $self->{stratify_field};
	unless ($field) { die; }  # Sanity check
	
	# Derive stratified prevalences for the target field
	# Create the mutation frequencies table with total + stratified frequences
	# Work out mutation frequencies, stratified according to options 
	my $message = "\n\t # Calculating mutation frequencies stratified by field '$field'";
	$io->show_output_message($message, $self->{output_type}); 

	$message = "\n\t # Deriving codon counts for field '$field'";
	$io->show_output_message($message, $self->{output_type}); 
	$alignment->derive_codon_counts($field);
	
	$message = "\n\t # Deriving mutation counts for field '$field'";
	$io->show_output_message($message, $self->{output_type});
	$alignment->derive_mutation_counts($field);
	
	$message = "\n\t # Deriving mutation frequencies for field '$field'";
	$io->show_output_message($message, $self->{output_type});
	$alignment->derive_mutation_frequencies($field);
	$alignment->create_mutation_frequency_table($field);
}

#***************************************************************************
# Subroutine:  derive polymorphism list
# Description: list of polymorphisms in MSA, based on a theshold frequency
#***************************************************************************
sub derive_polymorphism_list {

	my ($self) = @_;

	my $alignment = $self->{glue_msa_obj};
	unless ($alignment) { die; }  # Sanity check
	my $threshold = $self->{polymorphism_threshold};
	unless ($threshold) { die; }  # Sanity check
	my $message = " # Deriving polymorphism list with threshold frequency '$threshold'";
	$io->show_output_message($message, $self->{output_type}); 
	my %typical;	
	$alignment->derive_typical_list(\%typical, $threshold);
	$alignment->{typical_list} = \%typical;

}

############################################################################
# Phylogeny
############################################################################

#***************************************************************************
# Subroutine:  create_phylogeny 
# Description: create a distance tree using submitted sequences and a 
#              reference sequence set
#***************************************************************************
sub create_phylogeny {
	
	my ($self) = @_;

	# Override (if option for phylogenetic analysis not set)
	unless ($self->{phylogeny}) { return; }
	
	# Get relevant objects from self
	my $phylogeny_obj = $self->{phylogeny_obj};
	my $output_path   = $self->{report_dir};
	my $alignment     = $self->{glue_msa_obj};
	$io->show_output_message("\n\t # Creating phylogeny", $self->{output_type});	
	
	# Load the reference MSA and stack on the query sequence MSA
	my $fasta_path;
	if ($self->{include_glue_refset}) {
		#my $rss_path      = $self->{rss_path}; # Get path to the RSS
		#my $refseq        = $self->{refseq};
		#my $refseq_name   = $refseq->{name};
		#my $path = $rss_path . $refseq_name . '.fas';
		#$alignment->load_reference_set($path);
		# Stack with reference GLUE aligns
		# Write the alignment in minimal FASTA format (no leading or trailing)
		#my $fasta_file = 'GLUE_MSA_stacked.txt';
		#$fasta_path = $output_path . $fasta_file;	
		#my $flag = 'true';
		#$alignment->write_glue_alignment($fasta_path, $flag);
		#my $stacked_aln_ref = $alignment->{sequences};
		#$seqio->read_fasta($fasta_path, \@sequences);
		#print "\n\t ALIGN FILE $align_file";
		#$devtools->print_array(\@stacked_aln); die;
		#$devtools->print_array(\@rss_seqs);
		die;
	}

	# Write sequences as PHYLIP	
	my @translation;
	my @phylip;
	my $datatool = DataTool->new();
	my $sequences = $self->{sequences};
	my %translation;
	$datatool->fasta_to_phylip($sequences, \@phylip, \@translation, \%translation);
	my $align_file = $output_path . 'FASTME_input_MSA.phy'; 
	$fileio->write_file($align_file, \@phylip);
	
	# Build the tree
	my %tree_output;
	$phylogeny_obj->build_nj_tree($align_file, \%tree_output); # Create the tree
	#$devtools->print_hash(\%tree_output); die;

	# Get the newick tree	
	my $newick = $tree_output{newick};

	# Get the newick and back translate to the orginal seq IDs
	my $back_newick = $bioio->back_translate_newick($newick, \%translation);
	$alignment->{newick} = $back_newick;

	# Write tree to file
	my $refseq_name = $self->{refseq}->{name};
	my $newick_file = 'submitted_sequences.tre';
	my $newick_path = $output_path . $newick_file;
	my $tree_string = "tree $refseq_name = [&R] " . $newick;
	$fileio->write_text_to_file($newick_path, $tree_string);
	my $back_newick_file = 'submitted_sequences_backn.tre';
	my $back_newick_path = $output_path . $back_newick_file;
	my $back_tree_string = "tree $refseq_name = [&R] " . $back_newick;
	$fileio->write_text_to_file($back_newick_path, $back_tree_string);
	$alignment->{tree_file} = $newick_file;
	$alignment->{tree_path} = $newick_path;
}

############################################################################
# SECTION: FXNS FOR INSTANTIATING & WRITING GLUE MSA OBJECTS 
############################################################################

#***************************************************************************
# Subroutine:  create_msa_object
# Description: 
#***************************************************************************
sub create_msa_object {

	my ($self, $alignment_seqs) = @_;
	
	my %msa_data;
	
	# Process sequences
	my $lowest_start = undef;
	my $highest_stop = 1;
	my %hash_sequences;
	my @sequence_ids;
	my %header_id_link;
	my %field_tracker;
	my %fields;
	$fields{sequence_id} = 1;
	my %states_hash;
	$states_hash{total} = 1;
	$field_tracker{sequence_id} = \%states_hash;

	# Get alignment coordinates
	foreach my $seq_ref (@$alignment_seqs) {
		my $seq_id    = $seq_ref->{sequence_id};
		my $header    = $seq_ref->{header};
		my $sequence  = $seq_ref->{sequence};
		my $aln_start = $seq_ref->{aln_start};
		my $aln_stop  = $seq_ref->{aln_stop};
		unless ($seq_id and $header and $sequence and $aln_start and $aln_stop) {
			$devtools->print_hash($seq_ref); die;
		}
		unless ($lowest_start) { 
			$lowest_start = $aln_start; 
		}
		if ($aln_start < $lowest_start) {
			$lowest_start = $aln_start;
		}
		if ($aln_stop > $highest_stop) {
			$highest_stop = $aln_stop;
		}
		push (@sequence_ids, $seq_id);
		$hash_sequences{$seq_id} = $seq_ref;
		$header_id_link{$seq_id} = $header;
		my $data_ref = $seq_ref->{data};	
		#$devtools->print_hash($seq_ref); die;
		#$devtools->print_hash($data_ref); die;

		if ($data_ref) {
			my @fields = keys %$data_ref;
			foreach my $field (@fields) {
				$fields{$field} = 1;
				my $state = $data_ref->{$field};
				#print "\n\t # Field $field = $state";
				if ($field_tracker{$field}) {
					my $states_ref = $field_tracker{$field};
					$states_ref->{$state} = 1;
				}
				else {
					my %states;
					$states{$state} = 1;
					$field_tracker{$field} = \%states;
				}
			}		
		}
	}

	# Set start and stop coordinates and pad sequence
	foreach my $seq_ref (@$alignment_seqs) {
		
		my $seq_id    = $seq_ref->{sequence_id};
		my $header    = $seq_ref->{header};
		my $sequence  = $seq_ref->{sequence};
		my $aln_start = $seq_ref->{aln_start};
		my $aln_stop  = $seq_ref->{aln_stop};
		unless ($seq_id and $header and $sequence) {
			die;
		}
		#my $aln_seq = $seq_ref->{aln_seq};
	
		# Add the leading gap
		my $padded = $seq_ref->{padded};
		my $aln_seq;
		unless ($padded) {
			my $leader_len = $aln_start - $lowest_start;
			my $leader = '-' x $leader_len;
			$aln_seq = $leader . $sequence;
			#print "<br> leader length = $aln_start - $start = $leader_len";
			# Add the trailing gap
			my $trailer_len = $highest_stop - $aln_stop;
			my $trailer = '-' x $trailer_len;
			$aln_seq = $aln_seq . $trailer;
			#print "<br> trailer length = $stop - $aln_stop = $trailer_len";
		}
		else {
			$aln_seq = $sequence;
		}
		# Update the sequence settings
		#$seq_ref->{aln_seq}   = $aln_seq;
		$seq_ref->{sequence}  = $aln_seq;
		$seq_ref->{aln_start} = $lowest_start;
		$seq_ref->{aln_stop}  = $highest_stop;
	}

	my @fields_array = keys %fields;
	#print "\n\t LOWEST  START = $lowest_start, HIGHEST STOP = $highest_stop";
	$msa_data{lowest_start} = $lowest_start;
	$msa_data{highest_stop} = $highest_stop;
	$msa_data{msa_length}   = $highest_stop - $lowest_start + 1;

	# Refseq
	my $refseq = $self->{refseq};
	unless ($refseq) { die; }
	$refseq->{start}        = $lowest_start;
	$refseq->{stop}         = $highest_stop;
	
	# Store data
	$msa_data{refseq}         = $refseq;
	$msa_data{sequences}      = $alignment_seqs;
	$msa_data{hash_sequences} = \%hash_sequences;
		
	# Linking sequences to data
	$msa_data{sequence_ids}   = \@sequence_ids;
	$msa_data{header_id_link} = \%header_id_link;
	$msa_data{fields_array}   = \@fields_array;
	$msa_data{field_tracker}  = \%field_tracker;
	
	# Create GLUE MSA object
	my $alignment = RefSeqAlignment->new(\%msa_data);
	
	# Store GLUE MSA object
	$self->{glue_msa_obj}     = $alignment;

	#$devtools->print_hash(\%field_tracker); die;	
	#$devtools->print_array(\@sequence_ids); exit;
	#$devtools->print_hash(\%hash_sequences); exit;
	#$alignment->describe(); die; #DEBUG
}

############################################################################
# SECTION: FXNS FOR CREATING MSAs
############################################################################

#***************************************************************************
# Subroutine:  glue_blast_align 
# Description: create GLUE pairwise BLAST alignment
#***************************************************************************
sub glue_blast_align {

	my ($self, $seq_ref) = @_;

	if ($self->{override_blast}) { return; }
		
	# Get data from self
	my $blast_aligner  = $self->{blast_obj};
	my $refseq         = $self->{refseq};
	unless ($refseq) {die;} 

	# Set up for the alignments
	$seq_ref->{name}   = $refseq->{name};
	my $query_seq      = $seq_ref->{sequence};
	my $raw_seq_len    = length $query_seq;
	my $refseq_seq     = $refseq->{sequence};
	unless ($refseq_seq and $query_seq) {die;} 

	my $set = undef;
	if ($set) {
		
		# Make a file for BLAST
		#my $fasta      = ">$header\n$sequence";
		#my $query_file = $result_path . '/' . $id . '.fas';
		#$fileio->write_text_to_file($query_file, $fasta);

		# Do the BLAST
		#my $result_file = $id . '.blast_result';
		#my $path = $result_path . $result_file;
		#$blast_obj->{blast_bin_path} = './bin/blast/';
		#$blast_obj->blast($blast_alg, $lib_path, $query_file, $path, $settings_ref);
		#$blast_obj->parse_tab_format_step_one($path, $results_ref);
		#my $result_path = $blast_tool->assign($seq_ref, \@hits, \%settings);
		#$devtools->print_array(\@hits); #die;
		#$seq_ref->{result_path} = $result_path;
		#$seq_ref->{result} = \@hits;
	}
	else {
		# Do BLASTn pairwise align - uses XML BLAST output at present
		my $result_path = $self->{report_dir} . 'blast_result_';
		$blast_aligner->do_pairwise($query_seq, $refseq_seq, $result_path, $seq_ref);
		my $blast_ref = $seq_ref->{blast_align};
		my $valid = 1;
		unless ($blast_ref->{aln_seq}) {
			my $header   = $seq_ref->{header};
			my $sequence = $seq_ref->{sequence};
			$valid= undef;
		}
		$seq_ref->{blast_valid} = $valid;	
	}

}

#***************************************************************************
# Subroutine:  glue_muscle_align
# Description: Do MUSCLE pairwise align
#***************************************************************************
sub glue_muscle_align {

	my ($self, $seq_ref) = @_;
	unless ($self->{muscle_align}) { return; }
	#$devtools->print_hash($seq_ref); die;
	
	# Get data from self
	my $refseq         = $self->{refseq};
	my $query_seq      = $seq_ref->{sequence};
	my $refseq_seq     = $refseq->{sequence};
	$seq_ref->{name}   = $refseq->{name};
	
	# Write the raw sequences file with the refseq first
	my $muscle_aligner = $self->{muscle_obj};
	my @pair;
	my $refseq_fasta = ">REF\n$refseq_seq\n";
	my $query_fasta  = ">QUERY\n$query_seq\n";
	push(@pair, $refseq_fasta);
	push(@pair, $query_fasta);
	
	my $result_path = $self->{report_dir};
	$muscle_aligner->do_pairwise_aln(\@pair, $result_path, $seq_ref);
	my $muscle_ref = $seq_ref->{muscle_align};
	#$devtools->print_hash($muscle_ref);
	#print "\n\n";
	my $aln_seq = $muscle_ref->{aln_seq};
	
	my $valid = 1;
	unless ($aln_seq) {
		my $header = $seq_ref->{header};
	 	my $message = "Failed to align sequence '$header'";
		$io->show_output_message($message, $self->{output_type});
		$valid = undef;
	}
	$seq_ref->{muscle_valid} = $valid;	
}


#***************************************************************************
# Subroutine:  glue_lap_align 
# Description: Do Lap or tBLASTn align on ORFs in the reference  sequence 
#***************************************************************************
sub glue_lap_align {

	my ($self, $seq_ref) = @_;
	
	unless ($self->{lap_align}) { return; }
	if ($self->{override_blast}) { 
		die "\n\t Can't use 'override BLAST' setting with LAP"; 
	}
	unless ($seq_ref->{blast_valid}) { return; } 

	# Get/create flags, data, objects for doing the align
	my $output_type   = $self->{output_type},
	my $lap_aligner   = $self->{lap_obj};
	my $query_seq     = $seq_ref->{sequence};
	my $refseq        = $self->{refseq};
	my $refseq_seq    = $refseq->{sequence};
	my $seq_obj       = Sequence->new();
	unless ($refseq_seq and $query_seq) { die; } # Sanity checking
	
	# Get alignment details from BLASTn align
	my $blast_ref       = $seq_ref->{blast_align};
	my $blast_aln_start = $blast_ref->{aln_start};
	my $blast_aln_stop  = $blast_ref->{aln_stop};

	# CREATE THE TEMPLATE FOR LAP
	my $genes_ref = $refseq->{genes};	
	my %align_genes;
	foreach my $gene_ref (@$genes_ref) {
		my $exons_ref = $gene_ref->{exons};
		my @exon_starts = keys %$exons_ref;
		my $num_exons = scalar @exon_starts;
		if ($num_exons > 1) { next; }
		my $gene_start = $gene_ref->{start};
		my $gene_stop  = $gene_ref->{stop};
		my $gene_name  = $gene_ref->{name};
		#print "\n\t $gene_name: $aln_start-$aln_stop and $gene_start-$gene_stop";
		# To do (add conditions)
		$align_genes{$gene_name} = $gene_ref;
		#$devtools->print_hash($gene_ref);
	}
	
	# ALIGN EACH ORF
	my %translated_orfs;
	$refseq->get_translated_orfs(\%translated_orfs);
	my @orf_names = keys %translated_orfs;
	my %orfs_by_start;
	my %orf_aligns;
	foreach my $orf_name (@orf_names) {
		unless ($align_genes{$orf_name}) { next; }
		my $aa_seq = $translated_orfs{$orf_name};
		my %alignment;
		my $result_path = $self->{report_dir} . 'lap_result';
		$lap_aligner->do_align($query_seq, $aa_seq, $result_path, \%alignment);
		my $orf_start     = $align_genes{$orf_name}->{start}; 
		my $orf_end       = $align_genes{$orf_name}->{stop}; 
		my $subject_start = $alignment{subject_start}; 
		my $subject_end   = $alignment{subject_end}; 
		print "\n\t $orf_name ";
		#$subject_start    = $subject_start * 3;
		#$subject_end      = $subject_end * 3;
		#my $actual_start  = $subject_start + $orf_start;
		#my $actual_end    = $subject_end + $orf_end;
		#$alignment{actual_start} = $actual_start; 
		#$alignment{actual_end}   = $actual_end; 
		$orf_aligns{$orf_name} = \%alignment; 
		#print $query_seq; 
	}
	#die;

	# LAP 
	my %consolidated;
	my $seq_counter       = 0;
	my $last_aln_stop     = undef;
	my $last_aln_start    = undef;
	my $first_start       = undef;
	foreach my $gene_ref (@$genes_ref) {
		
		my $orf_name = $gene_ref->{name};
		unless ($align_genes{$orf_name}) { next; }
		#print "\n\t ORF GENE $orf_name"; 
		
		my $orf_start = $gene_ref->{start}; 
		my $orf_stop  = $gene_ref->{stop}; 
		my $orf_align = $orf_aligns{$orf_name};
		unless ($seq_counter) {
			$seq_counter = $orf_start;
			unless ($seq_counter) { die; } # Sanity checking
			$first_start = $orf_start;
		}
		else {
			$seq_counter = $last_aln_stop;
		}
		#$devtools->print_hash($orf_align); die;
		#my @keys = sort by_number keys %$align_sequence;
		#$align_sequence = $orf_align->{subject_hash};
		my $align_sequence = $orf_align->{aln_seq};
		my @align = split('', $align_sequence);
		my $i = 0;
		do {
			#print "\n\t $orf_name: seq_counter ($seq_counter), ($i), ($orf_stop)";
			if ($seq_counter >= $orf_start) {
				my $nt = $align[$i];
				#print " adding ($nt)";
				$consolidated{$seq_counter} = $nt; 
			}
			else {
				#print " adding (gap)";
				$consolidated{$seq_counter} = '-'; 
			}
			$i++;
			$seq_counter++;
			if ($i > 10000) { die; }
		} until ($seq_counter eq $orf_stop);
		
		# Store exit varaibles
		$seq_counter++;
		$last_aln_start = $orf_start;
		$last_aln_stop  = $orf_stop;
	}
	#die;
	#$devtools->print_hash(\%consolidated); die;
		

	my $seq_id = $seq_ref->{header};
	my $sequence;
	my @keys = sort by_number keys %consolidated;
	my $i = undef;
	foreach my $key (@keys) {
		unless ($i) { $i = $key }
		else        { $i++;     }
		#print "\n position $i:";
		unless ($i eq $key) { die "\n\t mistep $seq_id: ($i-$key)\n\n\n"; }
		my $nt = $consolidated{$key};
		#print "\t '$nt' ($key);";
		unless ($nt) { die "\n\t no $nt at $i\n\n\n"; }
		$sequence .= $nt;
	}
	#print "\n\t >GLUE\n$sequence"; exit;
	#my $aln_start = $orf_align->{subject_start};
	#my $aln_stop  = $orf_align->{subject_stop};
	#my $orf_end   = $gene_ref->{stop}; 
	#print "\n\t $orf_name";
	#$devtools->print_hash(\%alignment); die;
	#$lap_aligner->consolidate_align($refseq, \%orf_aligns, \%consolidated_align);
	my %consolidated_align;
	$consolidated_align{aln_seq}   = $sequence;
	$consolidated_align{aln_start} = $first_start;
	$consolidated_align{aln_stop}  = $last_aln_stop;
	$seq_ref->{lap_align} = \%consolidated_align;
	$seq_ref->{lap_valid} = 1;
}

############################################################################
# BASIC FXNS
############################################################################

#***************************************************************************
# Subroutine:  set_up_rolling_output 
# Description: 
#***************************************************************************
sub set_up_rolling_output {

	my ($self, $num_seqs, $output_type) = @_;
	
	# Set value for next increment
	my $next_increment;
	if   ($num_seqs < $increment) { 
		$next_increment = $num_seqs 
	}
	else { 
		$next_increment = $increment; 
	}
	if ($output_type eq 'html') { 
		print "Aligning $num_seqs sequences ";
		print " (please allow time for completion)";
	}
	else { print "\n\n\t ### Aligning $num_seqs sequences \n"; }
	return $next_increment;
}

#***************************************************************************
# Subroutine:  show_rolling_output 
# Description: 
#***************************************************************************
sub show_rolling_output {

	my ($self, $start, $next_increment, $counter, $num_seqs, $output_type) = @_;
	
	if ($counter eq 1) {
		if ($output_type eq 'html') { 
			print "</p><p>Aligning $counter-$next_increment ."; 
		}
		else { 
			print "\n\t Aligning $counter-$next_increment .";   
		}
	}
	elsif ($counter eq $next_increment and $counter < $num_seqs) {
		# Need a new line for rolling output
		$next_increment = $counter + $increment;
		if ($next_increment < $num_seqs) {
			if ($output_type eq 'html') { 
				print "</p><p>Aligning $counter-$next_increment ."; 
			}
			else { 
				print "\n\t Aligning $counter-$next_increment .";   
			}
		}
		else {
			if ($output_type eq 'html') { 
				print "</p><p>Aligning $counter-$num_seqs ."; 
			}
			else { 
				print "\n\t Aligning $counter-$num_seqs .";   
			}
		}
	}
	else {  print ".";  }
	return $next_increment;
}

############################################################################
# DEV
############################################################################

#***************************************************************************
# Subroutine:  reconcile_alignments 
# Description: 
#***************************************************************************
sub reconcile_alignments {

	my ($self, $seq_ref, $valid_ref, $failed_ref) = @_;
	
	# BLAST ALIGN
	my $align_ref;
	if ($seq_ref->{blast_valid}) { 
		
		my $blast_ref = $seq_ref->{blast_align};
		my $aln_start = $blast_ref->{aln_start};
		my $aln_stop  = $blast_ref->{aln_stop};
		#$devtools->print_hash($blast_ref);
		#die;
		$seq_ref->{aln_seq}     = $blast_ref->{aln_seq};
		$seq_ref->{aln_start}   = $aln_start;
		$seq_ref->{aln_stop}    = $aln_stop;
		$seq_ref->{insertions}  = $blast_ref->{insertions};
		push(@$valid_ref, $seq_ref);	
	}
	else { 
		push(@$failed_ref, $seq_ref);	
	}

	# MUSCLE ALIGN
	if ($self->{muscle_align} and $seq_ref->{muscle_valid}) { 
		my $muscle_ref = $seq_ref->{muscle_align};
		if ($seq_ref->{muscle_valid}) { 
			my $aln_seq   = $muscle_ref->{aln_seq};
			my $aln_start = $muscle_ref->{aln_start};
			my $aln_stop  = $muscle_ref->{aln_stop};
			#$devtools->print_hash($muscle_ref); die;
			$seq_ref->{aln_start}   = $aln_start;
			$seq_ref->{aln_stop}    = $aln_stop;
			$seq_ref->{insertions}  = $muscle_ref->{insertions};
			$seq_ref->{aln_seq}     = $aln_seq;
			delete $seq_ref->{muscle_align};
			push(@$valid_ref, $seq_ref);	
		}
	}
	# LAP ALIGN
	#elsif ($self->{lap_align} and $seq_ref->{lap_valid}) { 
	#	my $lap_ref = $seq_ref->{lap_align};
	#	unless ($lap_ref->{aln_seq}) { die; }
	#	$seq_ref->{aln_seq}     = $lap_ref->{aln_seq};
	#	$seq_ref->{aln_start}   = $lap_ref->{aln_start};
	#	$seq_ref->{aln_stop}    = $lap_ref->{aln_stop};
	#	$seq_ref->{insertions}  = $lap_ref->{insertions};
	#	push(@$valid_ref, $seq_ref);	
	#}
	#else { die; }
	
	$seq_ref->{sequence} = $seq_ref->{aln_seq};
	#delete $seq_ref->{aln_seq};
	return;
}

#***************************************************************************
# Subroutine:  by number
# Description: by number - for use with perl 'sort'  (cryptic but works) 
#***************************************************************************
sub by_number { $a <=> $b }	

############################################################################
# EOF
############################################################################
