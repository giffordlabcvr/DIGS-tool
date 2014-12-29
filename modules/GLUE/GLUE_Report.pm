#!/usr/bin/perl -w
############################################################################
# Script:       GLUE_Report.pm 
# Description: 	GLUE_Report object for writing results of GLUE anlayses 
# History:      Rob Gifford, June 2010: Creation
############################################################################
package GLUE_Report;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::FileIO;
use Base::HTML_Utilities;
use Base::DevTools;

############################################################################
# Globals
############################################################################



my $io         = IO->new();
my $devtools   = DevTools->new();
my $fileio     = FileIO->new();
my $html_utils = HTML_Utilities->new();
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create a new GLUE_Report.pm object
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref, $site_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Member variables
	my $self = {
	
		# Query data set
		sequence_ids            => $parameter_ref->{sequence_ids},
		sequences               => $parameter_ref->{sequences},
		raw_tab_data            => $parameter_ref->{raw_tab_data},
		num_raw_seqs            => $parameter_ref->{num_raw_seqs},
		set_name                => $parameter_ref->{set_name},
		
		# Paths and constants
		process_id              => $parameter_ref->{process_id},
		output_path             => $parameter_ref->{output_path},
		report_dir              => $parameter_ref->{report_dir},
		tmp_path                => $parameter_ref->{tmp_path},
		refseq_lib_path         => $parameter_ref->{refseq_lib_path},
		refseq_use_path         => $parameter_ref->{refseq_use_path},
		tree_file               => $parameter_ref->{tree_file},
		tree_path               => $parameter_ref->{tree_path},
		
		# Settings
		output_type             => $parameter_ref->{output_type}, 
		report_style            => $parameter_ref->{report_style}, 
		msa_qa_stats            => $parameter_ref->{msa_qa_stats},
		mut_profile             => $parameter_ref->{mut_profile},
		stratify_field          => $parameter_ref->{stratify_field},
		mut_frequencies         => $parameter_ref->{mut_frequencies},
		show_mutlist_seqs       => $parameter_ref->{show_mutlist_seqs},
		show_synonymous_changes => $parameter_ref->{show_synonymous_changes},
		all_mutations           => $parameter_ref->{all_mutations},
		derive_polymorphism_list=> $parameter_ref->{derive_polymorphism_list},
		phylogeny               => $parameter_ref->{phylogeny},
	
		# Settings
		glue_msa_obj            => $parameter_ref->{glue_msa_obj},
		glue_msa_file           => $parameter_ref->{glue_msa_file},
		glue_msa_path           => $parameter_ref->{glue_msa_path},

		# Website settings
		site_ref                => $site_ref,
	
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# Functions
############################################################################

#***************************************************************************
# Subroutine:  write_msa_report
# Description: write an HTML report providing a summary of a GLUE MSA 
#***************************************************************************
sub write_msa_report {

	my ($self) = @_;


	# Create the file for phylotype program
	#$self->create_phylotype_annotation();

	# Set up the html
	if ($self->{output_type} eq 'html') {
		$self->write_top_link_section();
	}

	# Write alignment statistics
	if ($self->{msa_qa_stats}) {
		#$self->write_alignment_statistics();
	}

	# Write summary of the dataset sequence by sequence
	if ($self->{mut_profile}) {
		$self->create_nonsynonymous_diversity_table();
		#print "BREAK"; die;
	}

	# Write synonymous change summary
	if ($self->{show_synonymous_changes}) {
		$self->create_synonymous_diversity_table();
	}

	# Special case
	if ($self->{all_mutations}) {
		$self->create_full_diversity_table();
	}

	# Write mutation frequencies
	if ($self->{mut_frequencies}) {
		$self->create_mutation_frequency_table();
	}

	# Write summary of the dataset sequence by sequence
	if ($self->{show_mutlist_seqs}) {
		$self->write_mutations_by_list_and_seq();
	}

	# Write the typical polymorphism list	
	if ($self->{derive_polymorphism_list}) { 
		$self->write_polymorphism_list();
	}

	if ($self->{output_type} eq 'html') {
		$self->write_bottom_link_section();
	}
}

#***************************************************************************
# Subroutine:  create_nonsynonymous_diversity_table
# Description: write summary of each sequence in dataset, including subtype
#              and all mutations
#***************************************************************************
sub create_nonsynonymous_diversity_table {

	my ($self) = @_;

	# Get the GLUE alignment object
	my $alignment_obj = $self->{glue_msa_obj};
	unless ($alignment_obj) { die; }  # Sanity checking
	
	# Get all the required data structures
	my $refseq            = $alignment_obj->{refseq};
	my $sequences_ref     = $alignment_obj->{sequences};
	
	# Write the table columns based on which genes were present
	my @gene_names;
	$refseq->get_gene_names(\@gene_names);
	my @head_row_data = ( 'ID' ); 
	push (@head_row_data, @gene_names); 
	my $head_row = join("\t", @head_row_data);
	
	# Iterate through all sequences, writing data for those with mutations on the list 
	my @table;
	push (@table, "$head_row\n");
	foreach my $sequence_ref (@$sequences_ref) {
		
		# Get mutations and sequence ID for this sequence
		#$devtools->print_hash($sequence_ref); die; exit;
		#print "<br> Doing $sequence_id";
		my $sequence_id       = $sequence_ref->{header};
		my $mutations_ref     = $sequence_ref->{nonsyn_mutations};
		my $syn_mutations_ref = $sequence_ref->{syn_mutations};
		$sequence_id = $self->format_sequence_id($sequence_id);
		
		# Iterate through mutations
		my %table_muts;
		$self->format_table_muts(\%table_muts, $mutations_ref, $refseq);

		# Write the mutations for each gene
		my @seq_row_data = ( $sequence_id ); 
		foreach my $gene_name (@gene_names) {
			# Get the non-synonymous mutations
			my $gene_mutations_ref = $table_muts{$gene_name};	
			my $mutstring ;
			unless ( $gene_mutations_ref ) { $mutstring = ''; }
			else {  $mutstring = join(', ', @$gene_mutations_ref); }
			push ( @seq_row_data, $mutstring );
		}
		my $seq_row = join("\t", @seq_row_data);
		push (@table, "$seq_row\n"); 
	}

	# Write text table
	my $report_dir = $self->{report_dir};
	my $set_name   = $self->{set_name};
	my $path = $report_dir . "/mutations.txt";
	if ($set_name) {
		$path = $report_dir . "/$set_name" . '.mutations.txt';
	}
	$fileio->write_file($path, \@table);

	# Write html table
	if ($self->{output_type} eq 'html') {
		my @html_table;
		$html_utils->convert_text_table_to_html(\@table, \@html_table);
		my $html_ref = $self->{html_output};
		my $table_title = "\t\tGENETIC VARIATION BY SEQUENCE\n"; 
		push (@$html_ref, $table_title);
		push (@$html_ref, @html_table);
	}
}

#***************************************************************************
# Subroutine:  create_synonymous_diversity_table
# Description: write summary of each sequence in dataset, including subtype
#              and all mutations
#***************************************************************************
sub create_synonymous_diversity_table {

	my ($self) = @_;

	# Get the GLUE alignment object
	my $alignment_obj = $self->{glue_msa_obj};
	unless ($alignment_obj) { die; }  # Sanity checking
	
	# Get all the required data structures
	my $refseq            = $alignment_obj->{refseq};
	my $sequences_ref     = $alignment_obj->{sequences};
	
	# Write the table header
	my @gene_names;
	$refseq->get_gene_names(\@gene_names);
	my @head_row_data = ( 'ID' ); 
	
	# Set up columns based on which genes were present
	push (@head_row_data, @gene_names); 
	push ( @head_row_data, '#s' );
	push ( @head_row_data, '#ns' );
	push ( @head_row_data, '#c' );
	push ( @head_row_data, '#ts' );
	push ( @head_row_data, '#tv' );
	my $head_row = join("\t", @head_row_data);
	
	# Iterate through all sequences, writing data for those with mutations on the list 
	my @table;
	push (@table, "$head_row\n");
	foreach my $sequence_ref (@$sequences_ref) {
		
		# Get mutations and sequence ID for this sequence
		#my $sequence_id   = $sequence_ref->{sequence_id};
		my $sequence_id       = $sequence_ref->{header};
		my $syn_mutations_ref = $sequence_ref->{syn_mutations};
		$sequence_id = $self->format_sequence_id($sequence_id);

		# Iterate through mutations
		my %table_muts;
		foreach my $mutation_ref (@$syn_mutations_ref) { 
			
			my $refseq_codon = $mutation_ref->{ref_codon};
			my $position     = $mutation_ref->{position};
			my $codon        = $mutation_ref->{codon};
			my $gene_name    = $mutation_ref->{gene};
			my $mutation     = $refseq_codon . $position . $codon;
			my $f_mutation   = $mutation;
			
			# Group by gene
			if ($table_muts{$gene_name}) {
				my $gene_mutations_ref = $table_muts{$gene_name};
				push (@$gene_mutations_ref, $f_mutation);	
			}
			else {
				my @gene_mutations;
				push (@gene_mutations, $f_mutation);
				$table_muts{$gene_name} = \@gene_mutations;
			}
		}
			
		# Write the mutations for each gene
		my @seq_row_data = ( $sequence_id ); 
		foreach my $gene_name (@gene_names) {
			# Get the non-synonymous mutations
			my $gene_mutations_ref = $table_muts{$gene_name};	
			my $mutstring ;
			unless ( $gene_mutations_ref ) { $mutstring = ''; }
			else {  $mutstring = join(', ', @$gene_mutations_ref); }
			push ( @seq_row_data, $mutstring );
		}
		
		#$devtools->print_hash($sequence_ref); die;
		my $num_syn           = $sequence_ref->{num_syn};
		my $num_nonsyn        = $sequence_ref->{num_nonsyn};
		my $num_changes       = $sequence_ref->{num_changes};
		my $num_transitions   = $sequence_ref->{num_transitions};
		my $num_transversions = $sequence_ref->{num_transversions};
		unless ($num_syn)            { $num_syn           = '-'; }
		unless ($num_nonsyn)         { $num_nonsyn        = '-'; }
		unless ($num_changes)        { $num_changes       = '-'; }
		unless ($num_transitions)    { $num_transitions   = '-'; }
		unless ($num_transversions)  { $num_transversions = '-'; }
		push ( @seq_row_data, $num_syn );
		push ( @seq_row_data, $num_nonsyn );
		push ( @seq_row_data, $num_changes );
		push ( @seq_row_data, $num_transitions );
		push ( @seq_row_data, $num_transversions );
		
		# Create the row
		my $seq_row = join("\t", @seq_row_data);
		push (@table, "$seq_row\n"); 
	}
	
	# Write text table
	my $report_dir = $self->{report_dir};
	my $set_name   = $self->{set_name};
	my $path = $report_dir . "/syn_mutations.txt";
	if ($set_name) {
		$path = $report_dir . "/$set_name" . '.syn_mutations.txt';
	}
	$fileio->write_file($path, \@table);

	# Write html table
	if ($self->{output_type} eq 'html') {
		my @html_table;
		$html_utils->convert_text_table_to_html(\@table, \@html_table);
		my $html_ref = $self->{html_output};
		my $table_title = "\t\tSYNONYMOUS VARIATION BY SEQUENCE\n"; 
		push (@$html_ref, $table_title);
		push (@$html_ref, @html_table);
	}

}

############################################################################
#  GENETIC DIVERSITY BY SEQUENCE (MUTATION BY MUTATION)
############################################################################

#***************************************************************************
# Subroutine:  create_full_diversity_table
# Description:  
#***************************************************************************
sub create_full_diversity_table {

	my ($self) = @_;

	# Get the GLUE alignment object
	my $alignment_obj = $self->{glue_msa_obj};
	unless ($alignment_obj) { die; }  # Sanity checking
	
	# Get all the required data structures
	my $refseq            = $alignment_obj->{refseq};
	my $sequences_ref     = $alignment_obj->{sequences};

	# Write the table header
	my @head_row_data = qw [ ID Gene Con Pos AA Con Codon ];
	my @matrix;
	push (@matrix, 'A->T');
	push (@matrix, 'A->C');
	push (@matrix, 'A->G');
	push (@matrix, 'T->A');
	push (@matrix, 'T->C');
	push (@matrix, 'T->G');
	push (@matrix, 'C->A');
	push (@matrix, 'C->T');
	push (@matrix, 'C->G');
	push (@matrix, 'G->A');
	push (@matrix, 'G->T');
	push (@matrix, 'G->C');
	push (@head_row_data, @matrix);
	my $head_row = join("\t", @head_row_data); 

	# Initialize totals hash
	my %totals;
	foreach my $change (@matrix) {
		$totals{$change} = 0;
	}
	
	# Iterate through all sequences, writing data for those with mutations on the list 
	my @table;
	push (@table, "$head_row\n"); 
	foreach my $sequence_ref (@$sequences_ref) {
		
		# Get mutations and sequence ID for this sequence
		my $sequence_id       = $sequence_ref->{header};
		my $mutations_ref     = $sequence_ref->{nonsyn_mutations};
		my $syn_mutations_ref = $sequence_ref->{syn_mutations};
		$sequence_id = $self->format_sequence_id($sequence_id);

		# Iterate through mutations
		my %table_muts;
		my @all_mutations;
		push (@all_mutations, @$mutations_ref);
		foreach my $mutation_ref (@all_mutations) { 
			
			my $refseq_aa = $mutation_ref->{refseq_aa};
			my $position  = $mutation_ref->{position};
			my $aa        = $mutation_ref->{aa};
			my $gene_name = $mutation_ref->{gene};
			my $ref_codon = $mutation_ref->{ref_codon};
			my $codon     = $mutation_ref->{codon};
			$codon = $self->format_codon_change($codon, $ref_codon);
			my $mutation  = $refseq_aa . $position . $aa;
			my $mutation_key = $mutation_ref->create_mutation_key();
			
			# Write the mutation
			my @seq_row_data;
			push (@seq_row_data, $sequence_id); 
			push (@seq_row_data, $gene_name); 
			push (@seq_row_data, $refseq_aa); 
			push (@seq_row_data, $position); 
			push (@seq_row_data, $aa); 
			push (@seq_row_data, $ref_codon); 
			push (@seq_row_data, $codon); 
			foreach my $change (@matrix) {
				my $num_changes = $mutation_ref->{$change};
				unless ($num_changes) { 
					push (@seq_row_data, '-'); 
				}
				else {
					push (@seq_row_data, $num_changes); 
					$totals{$change} = $totals{$change} + $num_changes;
				}
			}
			my $seq_row = join("\t", @seq_row_data);
			push (@table, "$seq_row\n"); 
		}
	}
	
	# Write text table
	my $report_dir = $self->{report_dir};
	my $set_name   = $self->{set_name};
	my $path = $report_dir . "/full_diversity.txt";
	if ($set_name) {
		$path = $report_dir . "/$set_name" . '.full_diversity.txt';
	}
	$fileio->write_file($path, \@table);

	# Write html table
	if ($self->{output_type} eq 'html') {
		my @html_table;
		$html_utils->convert_text_table_to_html(\@table, \@html_table);
		my $html_ref = $self->{html_output};
		my $table_title = "\t\tCOMPLETE MUTATION DATA\n"; 
		push (@$html_ref, $table_title);
		push (@$html_ref, @html_table);
	}

}

############################################################################
# WRITE MUTATION FREQUENCIES
############################################################################

#***************************************************************************
# Subroutine:  create_mutation_frequency_table
# Description:  
#***************************************************************************
sub create_mutation_frequency_table {

	my ($self) = @_;

	# Get the GLUE alignment object
	my $alignment_obj = $self->{glue_msa_obj};
	my $stratify_field = $self->{stratify_field};
	unless ($alignment_obj) { die; }  # Sanity checking
	
	# Get alignment data
	my $freq_table_ref      = $alignment_obj->{frequency_table};
	my $refseq              = $alignment_obj->{refseq};
	my $field_tracker       = $alignment_obj->{field_tracker};
	my $hash_lists_ref      = $refseq->{hash_lists};
	my $listnames_ref       = $refseq->{listnames};
	my $refseq_genes        = $refseq->{genes};

	# Create table data rows
	#my $head_row = shift @$freq_table_ref;
	my @table;
	my $remove = shift @$freq_table_ref;
	my @header;
	$header[0] = 'Gene';
	$header[1] = 'Mutation type';
	$header[2] = 'Mutation';
	$header[3] = 'Reference AA';
	$header[4] = 'Position';
	$header[5] = 'AA';
	$header[6] = 'Frequency';
	my $header = join("\t", @header);

	my @fields = qw [ sequence_id ];
	if ($stratify_field) {
		push (@fields, $stratify_field);
	}
	foreach my $field (@fields) {
		
		if ($field eq 'sequence_id') {
			next;
		}
		else {
			#print "\n\t FIELD $field";
			my $states_ref = $field_tracker->{$field};
			my @states = keys %$states_ref;
			foreach my $state (@states) {
				my $column_name = "$field ($state)";
				#print "\n\t NAME ($column_name)";
				$header .= "\t$column_name";
			}
		}
	}
	push (@table, "$header\n");
	foreach my $row_string (@$freq_table_ref) {
		
		chomp $row_string;
		my @row_data = split("\t", $row_string);
		my $gene     = shift @row_data;
		my $ref_aa   = shift @row_data;
		my $position = shift @row_data;
		my $aa       = shift @row_data;
		my $mutation = $ref_aa . $position . $aa;
		my $mutation_key  = $gene . ':' . $ref_aa . ':' . $position . ':' . $aa;
		
		# Check what the deal is with this mutation and lists
		my $listname;
		my $list;
		my @lists;
		my $typical_flag = undef;
		my $on_list_flag = undef;
		while (( $listname, $list ) = each %$hash_lists_ref ) {
			
			my $list_hash = $list->{hash_list}; # IMPORTANT
			if ($list_hash->{$mutation_key} and $listname eq 'typical') {
				$typical_flag = 1;
			}
			elsif ($list_hash->{$mutation_key}) {
				#print "\n\t mutation '$mutation_key' IS on $listname list";
				$on_list_flag = 1;
				push (@lists, $listname);
			}
			else {
				#print "\n\t mutation '$mutation_key' not on $listname list";
			}
		}
		unless ($typical_flag or $on_list_flag) {
			#next;
			push (@lists, 'atypical');
		}
		my $list_string = join(',', @lists);
		
		# Create row
		my @new_row;
		push (@new_row, $gene);
		push (@new_row, $list_string);
		push (@new_row, $mutation);
		push (@new_row, $ref_aa);
		push (@new_row, $position);
		push (@new_row, $aa);
		foreach my $freq (@row_data) {
			push(@new_row, $freq); 
		}
		my $row = join("\t", @new_row);
		push (@table, "$row\n"); 
	}
	#unshift(@table, "$head_row\n");

	# Write text table
	my $report_dir = $self->{report_dir};
	my $set_name   = $self->{set_name};
	my $path = $report_dir . "/frequencies.txt";
	if ($set_name) {
		$path = $report_dir . "/$set_name" . '.frequencies.txt';
	}
	$fileio->write_file($path, \@table);

	# Write html table
	if ($self->{output_type} eq 'html') {
		my @html_table;
		my $width = 900;
		$html_utils->convert_text_table_to_html(\@table, \@html_table, $width);
		my $html_ref = $self->{html_output};
		my $table_title = "\t\tMUTATION FREQUENCIES\n"; 
		push (@$html_ref, $table_title);
		push (@$html_ref, @html_table);
	}
}

############################################################################
# SHOW SEQUENCES THAT HAVE LIST MUTATIONS
############################################################################

#***************************************************************************
# Subroutine:  write_mutations_by_list_and_seq 
# Description:  
#***************************************************************************
sub write_mutations_by_list_and_seq {

	my ($self) = @_;
	
	# Get the GLUE alignment object
	my $alignment_obj = $self->{glue_msa_obj};
	unless ($alignment_obj) { die; }  # Sanity checking
	
	# Get all the required data structures
	my $refseq            = $alignment_obj->{refseq};
	my $hash_lists_ref    = $refseq->{hash_lists};
	my $listnames_ref     = $refseq->{listnames};
	unless ($refseq) { die; }
	
	# Which lists to show
	my %display_lists;
	my $sequences_ref = $alignment_obj->{sequences};
	foreach my $sequence_ref (@$sequences_ref) {

		# Get mutations and sequence ID for this sequence
		my $sequence_id   = $sequence_ref->{header};
		my $mutations_ref = $sequence_ref->{nonsyn_mutations};
		
		# Iterate through lists
		my $listname;
		my $hash_list_ref;
		my $got_list_mutations = undef;
		while (( $listname, $hash_list_ref) = each %$hash_lists_ref) {
			
			# Don't include common polymorphisms
			my $hash_mutlist_ref = $hash_list_ref->{hash_list};
			if ($listname eq 'typical') { next; }
			my $got_list_mutations = undef;
			foreach my $mutation_ref (@$mutations_ref) { 
				my $mutation_key = $mutation_ref->create_mutation_key();
				if ($hash_mutlist_ref->{$mutation_key}) { 
					$display_lists{$listname} = 1; 
				}
			}	
		}
	}
	
	# Which lists to show
	my @display_lists = keys %display_lists;

	# Iterate through all sequences, writing data for those with mutations on the list 
	my $got_any_mutations = undef;
	my $got_seq_mutations = undef;
	my @table;
	foreach my $sequence_ref (@$sequences_ref) {

		my @seq_row; # Array to store row data
		
		# Get mutations and sequence ID for this sequence
		my $sequence_id   = $sequence_ref->{header};
		my $mutations_ref = $sequence_ref->{nonsyn_mutations};
		push (@seq_row, $sequence_id);

		# Iterate through lists
		foreach my $listname (@display_lists) {
			
			# Don't include common polymorphisms
			if ($listname eq 'typical') { next; }
			
			# Get the list data reference
			my $hash_list_ref    = $hash_lists_ref->{$listname};
			my $hash_mutlist_ref = $hash_list_ref->{hash_list};
			my @list_muts;
			my %table_muts;
			my $got_list_mutations = undef;
			foreach my $mutation_ref (@$mutations_ref) { 

				#print "<BR> $sequence_id:  $listname: $mutation_key";	
				my $mutation_key = $mutation_ref->create_mutation_key();
				unless ($hash_mutlist_ref->{$mutation_key}) { next; }
				my $refseq_aa = $mutation_ref->{refseq_aa};
				my $position  = $mutation_ref->{position};
				my $aa        = $mutation_ref->{aa};
				my $gene_name = $mutation_ref->{gene};
				my $mutation  = $refseq_aa . $position . $aa;
				push (@list_muts, $mutation);
				$got_any_mutations  = 'true';
				$got_list_mutations = 'true';
				$got_seq_mutations = 'true';
			}	
			
			my $list_muts_string = '';
			if ($got_list_mutations) { 
				$list_muts_string = join(',', @list_muts); 
			}	
			push (@seq_row, $list_muts_string);
		}
		if ($got_seq_mutations) {
			my $row = join("\t", @seq_row);
			push (@table, "$row\n"); 
		}
		$got_seq_mutations = undef;
	}
	
	unless ($got_any_mutations) {
		#push (@$html_ref, "<ul><p> There were no listed mutations </ul></p>");	
		#push (@$html_ref,'<div id="menubottom"></div>');
	}
	else {
		
		my @head_row;
		push (@head_row, 'Sequence ID');
		push (@head_row, @display_lists);
		my $head_row = join("\t", @head_row);
		unshift(@table, "$head_row\n");		

		# Write text table
		my $report_dir = $self->{report_dir};
		my $set_name   = $self->{set_name};
		my $path = $report_dir . "/list_muts_by_seq.txt";
		if ($set_name) {
			$path = $report_dir . "/$set_name" . '.list_muts_by_seq.txt';
		}
		$fileio->write_file($path, \@table);
	}

	# Write html table
	if ($self->{output_type} eq 'html') {
		my @html_table;
		$html_utils->convert_text_table_to_html(\@table, \@html_table);
		my $html_ref = $self->{html_output};
		my $table_title = "\t\tSEQUENCES WITH LIST MUTATIONS\n"; 
		push (@$html_ref, $table_title);
		push (@$html_ref, @html_table);
	}

}

#***************************************************************************
# Subroutine:  write_polymorphism_list 
# Description: 
#***************************************************************************
sub write_polymorphism_list {
	
	my ($self) = @_;

	# Get the GLUE alignment object
	my $alignment_obj = $self->{glue_msa_obj};
	unless ($alignment_obj) { die; }  # Sanity checking
	my $list_ref = $alignment_obj->{typical_list};
	#$devtools->print_hash($list_ref); die;
	
	# Write polymorphisms
	my @list;
	my $key;
	my $freq;
	while ( ( $key, $freq ) = each %$list_ref ) {
		
		my @key = split(':', $key);
		my $gene = $key[0];
		my $ref_aa = $key[1];
		my $position = $key[2];
		my $aa = $key[3];
		my $date = 'NULL';
		my $line = "typical\t$date\t$gene\t$ref_aa\t$position\t$aa\t$freq\n";
		push (@list, $line);	
	}

	# Write text table
	my $report_dir = $self->{report_dir};
	my $set_name   = $self->{set_name};
	my $path = $report_dir . "/typical_polymorphisms.txt";
	if ($set_name) {
		$path = $report_dir . "/$set_name" . '.typical_polymorphisms.txt';
	}
	$fileio->write_file($path, \@list);
	
}

############################################################################
# MSA QA 
############################################################################

#***************************************************************************
# Subroutine:  write_alignment_statistics
# Description: write summary stats for MSA
#***************************************************************************
sub write_alignment_statistics {

	my ($self) = @_;

	# Get the GLUE alignment object
	my $alignment_obj = $self->{glue_msa_obj};
	unless ($alignment_obj) { die; }  # Sanity checking
	
	# Create numseqs row
	my $submitted     = $self->{num_raw_seqs};
	my $aligned       = $alignment_obj->{num_seqs};
	my $total_nt      = $alignment_obj->{total_nt};
	my $msa_nt        = $alignment_obj->{msa_nt};
	my $insertions_nt = $alignment_obj->{insertion_nt};
	unless ($total_nt) { die; }

	# Fraction	
	my $msa_fraction   = $msa_nt / $total_nt * 100;
	my $f_msa_fraction = sprintf("%.2f", $msa_fraction);
	my $insertions_fraction   = $insertions_nt / $total_nt * 100;
	my $f_insertions_fraction = sprintf("%.2f", $insertions_fraction);

	# Write the table header
	my @table;
	my @head_row_data = ( 'Submitted', 'Aligned', 'Total NT'); 
	push (@head_row_data, 'NT in MSA' ); 
	push (@head_row_data, '%' ); 
	push (@head_row_data, 'NT in insertions' ); 
	push (@head_row_data, '%' ); 
	my $head_row = join("\t", @head_row_data);		
	push (@table, "$head_row\n");

	my @seqs_row  = ($submitted, $aligned, $total_nt);
	push (@seqs_row, $msa_nt);
	push (@seqs_row, $f_msa_fraction);
	push (@seqs_row, $insertions_nt);
	push (@seqs_row, $f_insertions_fraction);
	my $seq_row = join("\t", @seqs_row);		
	push (@table, "$seq_row\n");
	
	# Write text table
	my $report_dir = $self->{report_dir};
	my $set_name   = $self->{set_name};
	my $path = $report_dir . "/alignment_stats.txt";
	if ($set_name) {
		$path = $report_dir . "/$set_name" . '.alignment_stats.txt';
	}
	$fileio->write_file($path, \@table);
	
	if ($alignment_obj->{fields_array}) { 
		$self->show_data_linking_failures();
	}
}

#***************************************************************************
# Subroutine:  show_data_linking_failures
# Description: write a table showing the specific failures 
#***************************************************************************
sub show_data_linking_failures {

	my ($self) = @_;

	# Get the GLUE alignment object
	my $alignment_obj = $self->{glue_msa_obj};
	unless ($alignment_obj) { die; }  # Sanity checking
	
	# Record the results of this attempt to link
	my $linking_ref = $alignment_obj->{linking};
	my $num_raw_seqs      = $self->{num_raw_seqs};
	my $num_without_data  = $linking_ref->{num_without_data};
	my $num_without_seqs  = $linking_ref->{num_without_seqs};
	my $num_linked        = $linking_ref->{num_linked};
	my $num_unlinked      = $linking_ref->{num_unlinked};
	my $fields_ref        = $linking_ref->{fields};            
	#$devtools->print_hash($linking_ref); die;

	# Write the table header
	my @table;
	my @head_row_data;
	push (@head_row_data, 'Submitted' ); 
	push (@head_row_data, 'Linked' ); 
	push (@head_row_data, 'Unlinked' ); 
	push (@head_row_data, 'Seqs without data' ); 
	push (@head_row_data, 'Data without seqs' ); 
	my $head_row = join("\t", @head_row_data);		
	push (@table, "$head_row\n");

	my @row;
	push (@row, $num_raw_seqs);
	push (@row, $num_linked);
	push (@row, $num_unlinked);
	push (@row, $num_without_data);
	push (@row, $num_without_seqs);
	my $row = join("\t", @row);		
	push (@table, "$row\n");
	
	# Write text table
	my $report_dir = $self->{report_dir};
	my $set_name   = $self->{set_name};
	my $path = $report_dir . "/linking.txt";
	if ($set_name) {
		$path = $report_dir . "/$set_name" . '.linking.txt';
	}
	$fileio->write_file($path, \@table);

}

############################################################################
# HTML formatting
############################################################################

#***************************************************************************
# Subroutine:  write_top_link_section
# Description:  
#***************************************************************************
sub write_top_link_section {

	my ($self) = @_;

	# Get the GLUE alignment object
	my $glue_msa_obj = $self->{glue_msa_obj};
	unless ($glue_msa_obj) { die; }  # Sanity checking
	
	# Create array for the HTML output
	my @output;

	# Get website paths etc
	my $site_ref    = $self->{site_ref};
	#my $header_path = $site_ref->{header_path};
	my $header_path = './site/ssi/glue_report_head.html';
	# Write html table
	if ($self->{output_type} eq 'html') {
		my @page_header;
		$fileio->read_file($header_path, \@page_header);
		my $header = join ("\n", @page_header);
		push (@output, $header);
	}

	# create the linked divider for this section
	my $class   = 'h3';
	my $title   = "GLUE MULTIPLE SEQUENCE ALIGNMENT";
	#my $divider = $writer->create_divider('href',  $class, $title); 
	my $divider;
	push (@output, $divider);

	# Create link to alignment
	my $glue_msa_file = $self->{glue_msa_file};
	unless ($glue_msa_file) { die; }
	my $alignlink  = "<p>Click <a href='$glue_msa_file'>here</a>";
	$alignlink    .= " to retrieve the GLUE multiple sequence alignment (MSA)</p>";
	push (@output, $alignlink);

	# Create link to phylogeny (if set)
	if ($self->{phylogeny}) {
		
		# Link to fig tree applet
		# <applet ARCHIVE="figtreeapplet.jar" CODE="figtree/applet/FigTreeApplet.class" CODETYPE="java/application" WIDTH="600" HEIGHT="800"></applet>
		my $newick       = $glue_msa_obj->{newick};
		my $tree_path    = $glue_msa_obj->{tree_path};
		my $tree_file    = $glue_msa_obj->{tree_file};
		
		unless ($tree_path) { die; }
		$tree_path    =~ s/.//;
		my $tree_url  = "http://saturn.adarc.org/GLUE-dev/$tree_path";
		my $figtree  = '<applet ARCHIVE="../../../bin/figtree/figtreeapplet.jar"';
		$figtree    .= 'CODE="figtree/applet/FigTreeApplet.class" CODETYPE="java/applet"';
		$figtree    .= 'WIDTH="800" HEIGHT="700">';
		$figtree    .= '<param name="tree" value="' . $tree_url .'"/>';
		$figtree    .= '<param name="style" value="default" />';
		$figtree    .= 'Browser does not support Java</applet>';

		my $fileio = FileIO->new();
		my $report_dir  = $self->{report_dir};
		my $figtree_path = $report_dir . 'figtree.html';
		print "\n\t <BR><BR> write $figtree_path";
		$fileio->write_text_to_file($figtree_path, $figtree);
		#die;

		# Close summary table
		#my $table_close = $writer->close_table(); 
		#push (@$html_ref, $table_close);
		my $treelink  = "<p>Click <a href='$tree_file'>here</a> to retrieve the tree, ";
		$treelink .= "Click <a href='figtree.html' target='figtree'>here</a>";
		$treelink .= "  to view the tree</p>";
	
		#my $phylopath = $self->{phylotype_path};
		#my $phylolink  = "<p>Click <a href='$phylopath'>here</a>phylotyping file</p>";
		#push (@$html_ref, $phylolink);
		push (@output, $treelink);
	}

	$self->{html_output} = \@output;
}


#***************************************************************************
# Subroutine: 
# Description: 
#***************************************************************************
sub create_phylotype_annotation {

	my ($self) = @_;

	# Get all the required data structures
	my $hash_lists_ref    = $self->{hash_lists};
	my $sequences_ref     = $self->{aligned_sequences};
	#$devtools->print_hash($hash_lists_ref);

	my %list_columns;
	my %lists_present;
	my $listname;
	my $hash_list_ref;
	while (( $listname, $hash_list_ref) = each %$hash_lists_ref) {
		
		# Don't include the list of common polymorphisms
		if ($listname eq 'typical') { next; }
		
		my $got_list_mutations = undef;
		my %list_mutation_positions; #store positions of list mutations
		foreach my $sequence_ref (@$sequences_ref) {
		
			my $mutations_ref = $sequence_ref->{mutations};
			foreach my $mutation_ref (@$mutations_ref) { 

				my $mutation_key = $mutation_ref->create_mutation_key();
				unless ($hash_list_ref->{$mutation_key}) { next; }
				my $position_key = $mutation_ref->create_position_key();
				$got_list_mutations = 'true';
				$list_mutation_positions{$position_key} = 1;
			}	
		}
		$list_columns{$listname}  = \%list_mutation_positions;
		$lists_present{$listname} = $got_list_mutations;
	}

	# Create header row for phylotype file
	my @table;
	my @columns = qw [ sequence_id ];
	while (( $listname, $hash_list_ref) = each %$hash_lists_ref) {

		# Don't include list unless we've got some data for it
		unless ($lists_present{$listname}) { next; }
		push (@columns, $listname);
		my $list_columns_ref = $list_columns{$listname};
		my @positions = keys %$list_columns_ref;
		foreach my $position (@positions) {
			push (@columns, $position);
		}
		#push (@columns, "$listname all");
	}
	my $header_row = join(" , ", @columns);
	push (@table, "$header_row\n");
	#$devtools->print_array(\@table);
	
	# Write a row for each sequence
	foreach my $sequence_ref (@$sequences_ref) {
	
		my @row;
		my $sequence_id   = $sequence_ref->{header}; 
		#print "\n\t Doing '$sequence_id'";
		push (@row, "'$sequence_id'");

		my $mutations_ref = $sequence_ref->{hash_mutations};
		my $list_mutations_ref;
		#$devtools->print_hash($hash_lists_ref); exit;
		while (( $listname, $list_mutations_ref) = each %$hash_lists_ref) {
			
			if ($listname eq 'typical')  { next; }
			#print "\n\t\t list '$listname'";

			my @list_row;

			# Don't include list unless we've got some data for it
			unless ($lists_present{$listname}) { next; }
			my $list_columns_ref = $list_columns{$listname};
			my @positions = keys %$list_columns_ref;
		
			my $got_any_list_mutations = undef;
			my $mutation_string = '';
			my @by_position_values;
			foreach my $position (@positions) {
				
				#print "\n\t\t\t position '$position'";
				my $mutation_key;
				my $mutation_ref;
				my @position_aminos;
				my $got_any_position_mutations = undef;
				#$devtools->print_hash($mutations_ref); exit;
				my $ref_aa;
				my $pos;
				while ( ( $mutation_key, $mutation_ref ) = each %$list_mutations_ref) {
				
					unless ($mutations_ref->{$mutation_key}) { next; }
					#print "\n\t\t\t\t mutation key $mutation_key";
					#$devtools->print_hash($mutation_ref); exit;
					my $gene  = $mutation_ref->{gene_name};
					my $aa    = $mutation_ref->{aa};
					$ref_aa   = $mutation_ref->{reference_aa};
					$pos      = $mutation_ref->{position};
					my $position_key = $gene . ':' . $pos;
					if ($position_key eq $position) {
						#print "\n\t\t\t\t position mutation '$mutation_key'";
						push (@position_aminos, $aa);
						$got_any_position_mutations = 1;
						$got_any_list_mutations = 1;
					}
				}
				if ($got_any_position_mutations) {
					my $position_aminos = join('', @position_aminos); 
					my $value = $ref_aa . $pos . $position_aminos;
					push (@by_position_values, "'$value'");
				}
				else {
					push (@by_position_values, "'0'");
				}
			}

			
			if ($got_any_list_mutations) {
				push (@list_row, "'1'");
			}
			else {
				push (@list_row, "'0'");
			}
			push (@list_row, @by_position_values);
			#push (@list_row, $list_mut_string);
			push (@row, @list_row);
				
		}
		my $row = join(" , ", @row);
		push (@table, "$row\n");
	}

	# add query set values;
	my @rss_values;
	my $rss_ref    = $self->{rss};
	my $refset_ref = $rss_ref->{rss_array};	
	foreach my $sequence_ref (@$refset_ref) {
		
		my @row;
		my $sequence_id   = $sequence_ref->{header}; 
		push (@row, "'$sequence_id'");
		my $list_mutations_ref;
		while (( $listname, $list_mutations_ref) = each %$hash_lists_ref) {
			
			if ($listname eq 'typical')  { next; }
			my @list_row;
			# Don't include list unless we've got some data for it
			unless ($lists_present{$listname}) { next; }
			my $list_columns_ref = $list_columns{$listname};
			my @positions = keys %$list_columns_ref;
			my @by_position_values;
			foreach my $position (@positions) {
				push (@by_position_values, "'0'");
			}
			
			push (@list_row, "'0'");
			push (@list_row, @by_position_values);
			#push (@list_row, $list_mut_string);
			push (@row, @list_row);
				
		}
		my $row = join(" , ", @row);
		push (@table, "$row\n");
	}
	
	my $outfile = 'listtree.phylo';
	my $outfile_path = $self->{output_path} . $outfile;
	$fileio->write_output_file($outfile_path, \@table);
	$self->{phylotype_path} = $outfile;

	# Write clustering data out
	my $newick    = $self->{newick};
	#my $annotated_newick = $self->annotate_newick($newick, @table);
	#my $nexus_tree = "#NEXUS\nbegin trees;\ntree test = $annotated_newick\nend;";
	#my $tree_path = $self->{output_path} . 'phylo.tre';
	#$fileio->write_text_to_file($tree_path, $nexus_tree);

}

#***************************************************************************
# Subroutine:  write_bottom_link_section
# Description: 
#***************************************************************************
sub write_bottom_link_section {

	my ($self) = @_;

	# Write summary of the dataset sequence by sequence
	#$self->create_qa_section($alignment_obj, \@output);
	my $output_path = $self->{output_path};
	my $report_dir  = $self->{report_dir};
	my $footer_path  = $self->{footer_path};
	my $html_ref    = $self->{html_output};
	my @output;

	# Write the file
	my $report_path = $report_dir . 'GLUE_MSA_report.hmtl';
	$fileio->write_file($report_path, $html_ref);

	#print "\n\n";
	print "\n\n\t<br><br> Click <a href='$report_path'>here</a> to access the report";
	
	# Add the footer
	if ($footer_path) {
		$html_utils->write_ssi(\@output, $footer_path);
	}
}

#***************************************************************************
# Subroutine:  format_table_muts
# Description:  
#***************************************************************************
sub format_table_muts {

	my ($self, $table_muts, $mutations_ref, $refseq) = @_;

	# Get all the required data structures
	my $hash_lists_ref    = $refseq->{hash_lists};
	my $listnames_ref     = $refseq->{listnames};
	#$devtools->print_hash($hash_lists_ref); die;

	foreach my $mutation_ref (@$mutations_ref) { 
		
		my $refseq_aa = $mutation_ref->{refseq_aa};
		my $position  = $mutation_ref->{position};
		my $aa        = $mutation_ref->{aa};
		my $gene_name = $mutation_ref->{gene};
		my $mutation  = $refseq_aa . $position . $aa;
		my $mutation_key = $mutation_ref->create_mutation_key();
		#$devtools->print_hash($mutation_ref); die;
		#print "\n Got $mutation_key";

		# Format output using list data	
		my $f_mutation = $mutation;;
		my $on_list = undef;
		my $color;
		foreach my $listname (@$listnames_ref) {
		
			my $list_ref = $hash_lists_ref->{$listname};
			my $list     = $list_ref->{hash_list};
			#$devtools->print_hash($list_ref); die;
			unless ($list) { die; }
			$color = $list_ref->{color};
			unless ($color) { die; }
			if ($listname eq 'typical') {
				if ($list->{$mutation_key}) {
					# skip typical list
					$color = '';
					next; 
				}
				else {
					$color = 'darkred'; 
					$f_mutation = "<font color=$color>$mutation</font>";    
					$on_list = 'true';
				}
			}
			elsif ($list->{$mutation_key}) {
				$f_mutation = "<font color=$color>$mutation</font>";    
				$on_list = 'true';
				#die "\n\t LIST $listname";
				#print "<BR> $listname: $color, mutation '$f_mutation'";
			}
		}	
		unless ($on_list) {
			#$color = 'black'; 
			#$f_mutation = "<font color=$color>$mutation</font>";    
			$f_mutation = $mutation; # Default    
		}
		#print "<BR> color was '$color' pushing mutation '$f_mutation'";

		# Group by gene
		if ($table_muts->{$gene_name}) {
			my $gene_mutations_ref = $table_muts->{$gene_name};
			push (@$gene_mutations_ref, $f_mutation);	
		}
		else {
			my @gene_mutations;
			push (@gene_mutations, $f_mutation);
			$table_muts->{$gene_name} = \@gene_mutations;
		}
	}
}

#***************************************************************************
# Subroutine:  format_sequence_id 
# Description:  
#***************************************************************************
sub format_sequence_id {

	my ($self, $sequence_id) = @_;

	# Fixing seq id length for table display
	my $id_len = length $sequence_id;
	if ($id_len > 25) {
		my @id = split ('', $sequence_id);
		my $new_id;
		my $i = 0;
		do {
			$new_id .= $id[$i];
			$i++;
		} until ($i eq 29);
		$sequence_id = $new_id . '...';	
	}

	return $sequence_id;
}

#***************************************************************************
# Subroutine:  format_codon_change 
# Description:  
#***************************************************************************
sub format_codon_change {

	my ($self, $codon, $ref_codon) = @_;

	my $f_codon;
	my @codon     = split ('', $codon);
	my @ref_codon = split ('', $ref_codon);
	my $i = 0;
	foreach my $nt (@codon) {
		my $f_nt   = $nt;
		my $ref_nt = $ref_codon[$i];
		if ($nt ne $ref_nt) {
			$f_nt = "<font color=red>$nt</font>";
		}
		$f_codon .= $f_nt;
		$i++;
	
	} 

	return $f_codon;
}

############################################################################
# EOF
############################################################################
