#!/usr/bin/perl -w
############################################################################
# Module:      RV_Views.pm
# Description: fxns for profiling reference sequence libraries and 
#              screening databases
# History:     November 2011: Created by Robert Gifford 
############################################################################
package RV_Views;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;
use DBI;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base Modules
use Base::Console;
use Base::DevTools;
use Base::FileIO;
use Base::BioIO;
use Base::SeqIO;

# Component Modules
use Component::GLUE::Sequence;
use Component::GLUE::RefSeq;
use Component::GLUE::RefSeqParser;
use Component::GLUE::RefSeqLibrary;

# Program Modules
use Program::Pipeline::DB;

############################################################################
# Globals
############################################################################

# TODO: replace this
my $indent  = ' ' x 10;

# Global site paths
my $site_name     = 'vglue_sandbox';                       
my $site_path     = './site/';                       
my $main_path     = './site/main/';                       
my $db_pages_path = './site/main/db_pages/';                       
my $url           = "http://saturn.adarc.org/vglue_sandbox/";

# Base objects
my $fileio    = FileIO->new();
my $devtools  = DevTools->new();
my $seqio     = SeqIO->new();
my $seq_obj   = Sequence->new();
my $console   = Console->new();
my $writer    = HTML_Utilities->new();
1;


############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Parameters
#***************************************************************************
sub new {

	my ($invocant, $parameters) = @_;
	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {
		
		# Flags
		mode                  => $parameters->{mode}, 
		
		# Site paths
		site_path             => $site_path,
		main_path             => $main_path,
		output_path           => $parameters->{output_path},
		
		# RefSeq paths
		refseq_lib_path       => $parameters->{refseq_lib_path},
		refseq_use_path       => $parameters->{refseq_use_path},
		
		# BLAST paths
		blast_db_path         => $parameters->{blast_db_path},
		blast_orf_lib_path    => $parameters->{blast_orf_lib_path},
		blast_env_lib_path    => $parameters->{blast_env_lib_path},
		blast_nt_orf_lib_path => $parameters->{blast_nt_orf_lib_path},
		blast_utr_lib_path    => $parameters->{blast_utr_lib_path},
		blast_genome_lib_path => $parameters->{blast_genome_lib_path},
		
		# Member classes 
		blast_obj             => $parameters->{blast_obj},
		views_obj             => $parameters->{views_obj},
		
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# SECTION: Retrovirus pages
############################################################################

#***************************************************************************
# Subroutine:  make_retrovirus_record_html 
# Description: make the retrovirus fossil record HTML 
#***************************************************************************
sub make_retrovirus_record_html {

	my ($self, $site_ref) = @_;

	# Load the retroviral reference sequences
	my %refseqs;
	my @refseqs;
	my $refseqlib_obj = RefSeqLibrary->new();
	my $libpath = $self->{refseq_lib_path};
	#my $libpath = "test/WEBSITE/db_tiny/refseq/VIRUS/";
	
	print "\n\t loading GLUE reference sequence library ' $libpath'";
	$refseqlib_obj->load_refseq_library($libpath, \@refseqs, \%refseqs);	

	# Get the virus data
	my %by_taxonomy;
	my %flat;
	my @ranks = qw [ virus_family virus_subfamily virus_supertribe 
                     virus_tribe virus_genus virus_subgroup name ];
	$refseqlib_obj->order_by_taxonomy(\@refseqs, \@ranks, \%by_taxonomy, \%flat);
	#$devtools->print_hash(\%by_taxonomy); die;

	# Set up for colecting statistics on the family
	my %retroviridae_data;
	$self->setup_rv_fossil_stats(\%retroviridae_data);

	# Get host taxonomy information
	print "\n\t Getting host taxonomic information";
	my $host_views = Host_Views->new($self);
	$host_views->setup_host_taxonomy(\%flat); 
	#$devtools->print_hash(\%flat); #die;

	# Iterate through the viruses, writing the page for each 
	my @viruses = sort keys %flat;
	foreach my $virus (@viruses) {
		
		my $refseq = $flat{$virus};
		unless ($refseq)  { die; }
		print "\n\t Processing '$virus'";
		
		# Get host lineage data  
		print "\n\t\t Getting '$virus' host lineage data";
		$self->stats_by_host_group($refseq, \%retroviridae_data);
		
		# Get statistics for this retroviral reference sequence
		print "\n\t\t Getting '$virus' statistics";
		$self->stats_by_virus_group($refseq, \%retroviridae_data);
		
		print "\n\t\t Writing fossil record page for '$virus'";
		my $erv_views = ERV_Views->new($self);
		#$erv_views->write_retrovirus_fossil_record_page($site_ref, $refseq);
		
		print "\n\t\t Writing Refseq reference page for '$virus'";
		$self->write_retrovirus_refseq_page($site_ref, $refseq);
	}
	
	# Write RV fossil index page
	$self->make_retrovirus_index($site_ref, \%flat, \%by_taxonomy);

	# Write the overview page
	$self->make_rv_fossil_overview($site_ref, \%retroviridae_data, \%flat);
	
	# Create the html
	my $views_obj = $self->{views_obj};
	$views_obj->make_paleo_website_html($site_ref);

	# Get site pages
	my $pages_array_ref = $site_ref->{pages_array};
	unless ($pages_array_ref) { die; }
	print "\n\n\t # Writing HTML\n";
	foreach my $page_ref (@$pages_array_ref) {
		my $page_name = $page_ref->{name};
		print "\n\t # Making $page_name HTML";
		my @html;
		$writer->assemble_page($page_ref, \@html);
		#$devtools->print_array(\@html);
		$writer->write_page($page_ref, \@html);
	}
}

#***************************************************************************
# Subroutine:  write_retrovirus_refseq_page 
# Description: write the reference sequence description page for an RV 
#***************************************************************************
sub write_retrovirus_refseq_page {

	my ($self, $site_ref, $refseq, ) = @_;

	# Output status and get refseq data 
	my $refseq_name = $refseq->{name};
	
	# Copy to the flat directory
	my $path = $refseq->{file_path};
	my $flat_dir = $self->{refseq_use_path};
	unless ($flat_dir) { die; }
	my $command = "cp $path $flat_dir";	
	system $command;
	
	# Create the main panel (left) for the refseq HTML
	my @refseq_main;
	$self->create_rv_refseq_view($refseq, \@refseq_main);

	# Write the main block (left panel) to the refseq HTML directory
	my $page_path = $db_pages_path . $refseq_name . '.html'; 
	$fileio->write_output_file($page_path, \@refseq_main);
	#print "\n\t PATH $page_path:";

	# Define the page
	my $full_name   = $refseq->{full_name};
	my $genus       = $refseq->{virus_genus};
	my %page;
	$page{title}  = "Paleovirology online: $full_name ($refseq_name)";
	$page{title} .= " reference sequence";
	my $meta_keywords = "Paleovirology, virus fossil record, $refseq_name";
	$meta_keywords .= ", $genus, GLUE, reference sequence";
	$page{meta_keywords} = $meta_keywords;
	my $description = "Reference sequence description page for $refseq_name";
	$description .= ", including genome features and putative protein sequences";
	$page{description} = $description;
	my $views_obj = $self->{views_obj};
	$views_obj->define_paleo_page($site_ref, \%page, $refseq_name, 3, 'fossil');
}

#***************************************************************************
# Subroutine:  make_retrovirus_index	
# Description: Fossil record HTML: Retroviridae index page
#***************************************************************************
sub make_retrovirus_index {

	my ($self, $site_ref, $flat_ref, $by_taxonomy_ref) = @_;

    # Create main index page
	my @main;
	my $indent  = ' ' x 10;
	my $version = '2013-05-25'; # TODO: set dynamically not hard-coded
	push (@main, "\n$indent <h3>Retrovirus reference sequence library (version $version)</h3>");
	push (@main, "\n$indent <div class='separator'></div>");
	push (@main, "\n$indent Click on virus names to view reference sequence details.");
	push (@main, "\n$indent For an explanation of the fossil record pages, click");
	push (@main, "\n$indent <a href=\"http://saturn.adarc.org/vglue_sandbox/site/html/RV_refseq_definitions.html\"><u>here</u></a>.<br><br>");
	#$devtools->print_hash($by_taxonomy_ref);

	my $retroviridae = $by_taxonomy_ref->{Retroviridae};
	my @subfamilies = sort keys %$retroviridae;
	foreach my $subfamily (@subfamilies) {
		push (@main, "\n$indent <h3> Subfamily $subfamily</h3>");
		push (@main, "\n$indent <div class='divider'></div>");
		my $supertribes_ref = $retroviridae->{$subfamily};
		my @supertribes = sort keys %$supertribes_ref;
		foreach my $supertribe (@supertribes) {
		
			my $f_supertribe = $supertribe;
			if ($supertribe eq 'Mintakaretrovirinae') { $f_supertribe .= ' (class II)';  }
			if ($supertribe eq 'Alnilamretrovirinae') { $f_supertribe .= ' (class I)';   }
			if ($supertribe eq 'Alnitakretrovirinae') { $f_supertribe .= ' (class III)'; }
			push (@main, "\n$indent <u><h4> Supertribe $f_supertribe</h4></u><br>");
			my $tribes_ref = $supertribes_ref->{$supertribe};
			my @tribes = sort keys %$tribes_ref;

			foreach my $tribe (@tribes) {
				# Write tribe table
				my $f_tribe = $tribe;
				my $genus_ref = $tribes_ref->{$tribe};
				my @genera = sort keys %$genus_ref;
				push (@main, "\n$indent <u><h4> Tribe $f_tribe</h4></u><br>");
                my $genus_count = 0;
				foreach my $genus (@genera) {
					$genus_count++;
					print "\n\t Creating index table for genus '$genus' in tribe '$tribe'";
					my @table;
					my $table_open = $writer->open_table($indent, 'text'); # Open summary table
					push (@table, $table_open);

					# column widths
					my $column_tags =  "<col width=\"40\" />
										<col width=\"25\" />
										<col width=\"25\" />";

#										<col width=\"15\" />
#										<col width=\"15\" />
#										<col width=\"15\" />";
					push (@table, $column_tags);

					my @head_row_data;
					push (@head_row_data, 'Name');
					push (@head_row_data, "Subgroup");
					#push (@head_row_data, "Provirus");
					#push (@head_row_data, "Fragment");
					#push (@head_row_data, "Solo LTR");
					#push (@head_row_data, "Status");
					push (@head_row_data, 'Accession #');
					my $head_row = $writer->create_header_row(\@head_row_data); 
					#push (@table, "\n$indent <br><h4> Genus $genus</h4>");
					push (@table, "\n$indent <h4> Genus $genus</h4>");
					push (@table, "\n$indent $head_row");


					my $subgroup_ref = $genus_ref->{$genus};
					my @subgroups = sort keys %$subgroup_ref;
					my $added_to_table = 0;
					foreach my $subgroup (@subgroups) {
						my $viruses_ref = $subgroup_ref->{$subgroup};
						foreach my $virus (@$viruses_ref) {
							print "\n\t Virus $virus: genus '$genus' in tribe '$tribe' and subgroup '$subgroup'";
							#print " '$virus'";
							my $refseq = $flat_ref->{$virus};
							my $state = $refseq->{genome_state};
							my $path  = $refseq->{path};
							if ($state) {
								if ($state eq 'complete') {
									print "\n\t State = '$state; for refseq under $path";
									next;
								}
							}
							else {
								sleep 1;
								print "\n\t No genome state for '$virus'";
							}
							# Format the name
							my $f_virus = $virus; 
							if ($state eq 'exogenous') {  $f_virus = "<u>$virus</u>"; }
							my $f_name  = "<a href=\"http://saturn.adarc.org/$site_name";
							$f_name .= '/site/html/fossil-record/';
							$f_name .= "$virus.html" . '">' . "$f_virus</a>";
							
							my @row;
							push (@row, $f_name);
							#push (@row, $genus);
							if ($subgroup eq 'Unclassified') { $subgroup = '-'; }
							push (@row, $subgroup);
							
							# Get the fossil record data for this reference sequence
							my %record_data;
							#$self->get_rv_record_data(\%record_data, $refseq);
							#my $solo           = $record_data{solo};
							#my $internal_noltr = $record_data{internal_noltr};
							#my $provirus       = $record_data{provirus};
							my $solo           = "Temp";
							my $internal_noltr = "Temp";
							my $provirus       = "Temp";
							#push (@row, $provirus);
							#push (@row, $internal_noltr);
							#push (@row, $solo);
							#push (@row, 'open');
							my $accession = $refseq->{accession};
							if ($accession eq 'NULL') {
                            	$accession = '-'
							}
                            push (@row, $accession);
								
							my $row = $writer->create_row(\@row); 
							push (@table, "\n$indent $row"); 
							$added_to_table++;
						}
						# NAH
						#push (@table, "\n$indent <div class='divider'></div>");
					}


					# Close table
					my $table_close = $writer->close_table($indent); 
					push (@table, $table_close);
					if ($added_to_table) { 
						if ($added_to_table eq 1 and $genus_count eq 1) {
						}
						push (@main, @table); 
					}
				}
			}
		}
	}
	
	push (@main, "\n$indent <div class='separator'></div>");

	# Write the main block (left panel) to the refseq HTML directory
	my $page_name = 'retroviruses';
	my $page_path = $db_pages_path . "$page_name.html"; 
	$fileio->write_output_file($page_path, \@main);
	
	# Define the page
	my %page;
	$page{title} = 'Paleovirology online: Retrovirus fossil record';
	my $views_obj = $self->{views_obj};
	$views_obj->define_paleo_page($site_ref, \%page, $page_name, 2, 'fossil');

}

#***************************************************************************
# Subroutine:  get_rv_record_data
# Description: 
#***************************************************************************
sub get_rv_record_data {

	my ($self, $record_data, $refseq) = @_;

	# Get refseq data
	my $name = $refseq->{name};
	
	# Get the right screening database 
	my $db_name = $self->select_screening_db_using_taxonomy($refseq);
	print "\n\t Database $db_name";
	
	# Instantiations
	my $seq_obj    = Sequence->new();
	my $db = DB->new();	
	my @html;
	$db->load_screening_db($db_name);
	my $loci_table = $db->{loci_table};
	
	# 1: numbers and relative proportions of each loci class in
	#$self->get_consolidate_erv_counts_by_refseq($db, \@html);	
	# 2: overall numbers of class I, II and III
	# 3: overall numbers of loci for each refseq
	# LOCI versus EXTRACTED
	# 1: High copy number matches to ERVs that are not in cognate set 
		
	print "\n\t ## Getting data $name";
	my $field = 'COUNT(*)';
	
	# Solo LTRs	
	my @total;
	my $total_where = "WHERE assigned_to = '$name'";
	$loci_table->select_value($field, \@total, $total_where);	
	my $total_count = shift @total;	


	# Solo LTRs	
	my @solo;
	my $solo_where = "WHERE assigned_to = '$name' 
					  AND   genome_structure  = 'LTR'";
	$loci_table->select_value($field, \@solo, $solo_where);	
	my $solo_count = shift @solo;	
	unless ($solo_count) { $solo_count = '0'; }
	my $solo_prop;
	if ($total_count) {
		$solo_prop = $solo_count / $total_count;
	}

	# Internals with LTRs	
	my @provirus;
	my $provirus_where  = "WHERE assigned_to = '$name'"; 
	   $provirus_where .= "AND   genome_structure  != 'LTR'
					  AND   genome_structure LIKE '\%LTR\%'";
	$loci_table->select_value($field, \@provirus, $provirus_where);	
	my $provirus_count = shift @provirus;	
	unless ($provirus_count) { $provirus_count = '0'; }
	my $provirus_prop;
	if ($total_count) {
		$provirus_prop = $provirus_count / $total_count;
	}
		
	# Internals without LTRs	
	my @internals;
	my $internals_where = "WHERE assigned_to = '$name' 
					  AND   genome_structure NOT LIKE '\%LTR\%'";
	$loci_table->select_value($field, \@internals, $internals_where);	
	my $internals_count = shift @internals;
		
	unless ($internals_count) { $internals_count = '0'; }
	my $internals_prop;
	if ($total_count) {
		$internals_prop = $internals_count / $total_count;
	}

	print "\t total ($total_count)";
	print "\t solo ($solo_count)\tprovirus ($provirus_count)\tinternal no LTR ($internals_count)";

	unless ($total_count) { $total_count = '0'; }
	unless ($solo_count)      { $solo_count      = 0; }
	unless ($internals_count) { $internals_count = 0; }
	unless ($provirus_count)  { $provirus_count  = 0; }

	# get data
	$record_data->{solo}           = $solo_count;
    $record_data->{internal_noltr} = $internals_count;
    $record_data->{provirus}       = $provirus_count;

}

#***************************************************************************
# Subroutine:  create_rv_refseq_view	
# Description: Create the main (left) panel content for the refseq HTML
#***************************************************************************
sub create_rv_refseq_view {

	my ($self, $refseq, $html_ref) = @_;

	# Get refseq data
	my $indent = ' ' x 10;
	my $refseq_name  = $refseq->{name};
	unless ($refseq_name) { die; }

	my %orfs;
	$refseq->get_translated_orfs(\%orfs);
	$refseq->{translated_orfs} = \%orfs;

	# Write the taxonomy section
    $self->write_retrovirus_header($refseq, $html_ref);

	# Write the features table
    $self->write_feature_table($refseq, $html_ref);

	# Make genome-specific gene feature table
	push (@$html_ref, "$indent<h3>Genome features</h3>");
	push (@$html_ref, "$indent<div class=divider></div>");
	$self->make_genome_feature_table($refseq, $html_ref);

	# Dinucleotide composition
	push (@$html_ref, "$indent<h3>Dinucleotide composition</h3>");
	push (@$html_ref, "$indent<div class=divider></div>");
	$self->make_dinuc_composition_table($refseq, $html_ref);

	# Write the gene sequences table
    $self->write_gene_sequences_table($refseq, $html_ref);
	push (@$html_ref, "\n$indent <div class='separator'></div>");
	#$devtools->print_array($html_ref); die;
}

############################################################################
# SECTION: Fossil record HTML: Fossil record overview statistics
############################################################################

#***************************************************************************
# Subroutine:  make_rv_fossil_overview 
# Description: 
#***************************************************************************
sub make_rv_fossil_overview {

	my ($self, $site_ref, $stats_ref, $flat_refseqs) = @_;

    # Create main index page
	my @main;
	my $indent = ' ' x 10;
	push (@main, "\n$indent <h3>Retrovirus fossil record: Overview</h3>");
	push (@main, "\n$indent <div class='separator'></div>");
	push (@main, "\n$indent For an explanation of the retrovirus reference sequence pages, click");
	push (@main, "\n$indent <a href=\"http://saturn.adarc.org/vglue_sandbox/site/html/RV_refseq_definitions.html\"><u>here</u></a>.<br><br>");

	# Total ERVs by class
	my @table1;
	$self->total_ervs_by_tribe(\@table1, $stats_ref);
	push (@main, @table1);
	
	# Total ERVs by genus
	my @table2;
	$self->total_ervs_by_genus(\@table2, $stats_ref);
	push (@main, @table2);

	# By host taxonomy
	my @species; 
	my @refseq_names = keys %$flat_refseqs;
	foreach my $name (@refseq_names) {
		my $refseq = $flat_refseqs->{$name};
		my $species = $refseq->{host_sci_name};
		unless ($species)  { die; }
		push (@species, $species);
	}
	my %sorted;
	my $host_views = Host_Views->new($self);
	$host_views->sort_hosts_by_taxonomy(\@species, \%sorted);
	my @table3;
	push (@main, "\n$indent <br><h3> Retrovirus by host taxonomy</h3>");
	push (@main, "\n$indent <br><div class='separator'></div>");
	my $by_superorder = $sorted{by_superorder_mammals};
	push (@table3, "\n$indent <br><h4> Class Mammalia (Superorders)</h4>");
	$self->rvs_by_host_taxonomy(\@table3, $stats_ref, $by_superorder);
	$by_superorder = $sorted{by_superorder_birds};
	push (@table3, "\n$indent <br><h4> Class Aves (Superorders)</h4>");
	$self->rvs_by_host_taxonomy(\@table3, $stats_ref, $by_superorder);

	# BY ORDER
	my $by_order = $sorted{by_order_mammals};
	push (@table3, "\n$indent <br><h4> Class Mammalia (Orders)</h4>");
	$self->rvs_by_host_taxonomy2(\@table3, $stats_ref, $by_order);
	$by_order = $sorted{by_order_birds};
	push (@table3, "\n$indent <br><h4> Class Aves (Orders)</h4>");
	$self->rvs_by_host_taxonomy2(\@table3, $stats_ref, $by_order);
	push (@main, @table3);
	
	# Write the main block (left panel) to the refseq HTML directory
	push (@main, "\n$indent <br><div class='separator'></div>");
	my $page_name = 'viral_fossil_record';
	my $page_path = $db_pages_path . "$page_name.html"; 
	$fileio->write_output_file($page_path, \@main);
	
	# Define the page
	my %page;
	$page{title} = 'Paleovirology online: Virus fossil record overview';
	my $views_obj = $self->{views_obj};
	$views_obj->define_paleo_page($site_ref, \%page, $page_name, 2, 'fossil');
}

#***************************************************************************
# Subroutine:  write_retrovirus_header 
# Description: 
#***************************************************************************
sub write_retrovirus_header {

	my ($self, $refseq, $html_ref) = @_;

	# Set variables / get refseq data
	my $indent = ' ' x 10;
	my $record_path     = $refseq->{record_path};
	my $name            = $refseq->{name};
	my $full_name       = $refseq->{full_name};
	my $genus           = $refseq->{virus_genus};
	my $subgroup        = $refseq->{virus_subgroup};
	my $virus_subfamily = $refseq->{virus_subfamily};
	my $virus_tribe     = $refseq->{virus_tribe};
	my $virus_class     = $refseq->{virus_class};
	my $state           = $refseq->{virus_state};
	my $accession       = $refseq->{accession};
	my $host_superclass = $refseq->{host_superclass};
	my $host_class      = $refseq->{host_class};
	my $host_superorder = $refseq->{host_superorder};
	my $host_order      = $refseq->{host_order};
	my $host_suborder   = $refseq->{host_suborder};
	my $host_family     = $refseq->{host_family};
	my $host_genus      = $refseq->{host_genus};
	my $host_species    = $refseq->{host_sci_name};
	my $host_common     = $refseq->{host_common_name};
	my $host_taxonomy   = $refseq->{taxonomy_string};
	unless ($name and $full_name) { 
		die;  # We should have all these variables
	}

	push (@$html_ref, "\n$indent<h3>$name reference sequence</h3>");
	push (@$html_ref, "\n$indent<div class=separator></div>");
	push (@$html_ref, "\n$indent For an explanation of the retrovirus reference sequence pages, click");
	push (@$html_ref, "\n$indent <a href=\"http://saturn.adarc.org/vglue_sandbox/site/html/RV_refseq_definitions.html\"><u>here</u></a>.<br><br>");
	
    push (@$html_ref, "\n$indent<h3>Virus taxonomy</h3>");
	push (@$html_ref, "\n$indent<div class=divider></div>");
	my $table_open = $writer->open_table($indent, 'text'); # Open summary table
	push (@$html_ref, $table_open);

	# Write taxonomy etc
	my @name_row;
	push (@name_row, "Virus name");
	push (@name_row, "$name ($full_name)");
	my $name_row = $writer->create_row(\@name_row); 
	push (@$html_ref, "\n$indent $name_row");
	
	my @family_row;
	push (@family_row, "Virus family");
	push (@family_row, "<i>Retroviridae</i>");
	my $family_row = $writer->create_row(\@family_row); 
	push (@$html_ref, "\n$indent $family_row");

	my @sub_family_row;
	unless ($virus_subfamily) {die "\n\t Subfamily undefined\n\n"; }
	push (@sub_family_row, "Virus subfamily");
	push (@sub_family_row, "<i>$virus_subfamily</i>");
	my $sub_family_row = $writer->create_row(\@sub_family_row); 
	push (@$html_ref, "\n$indent $sub_family_row");

	my @tribe_row;
	unless ($virus_tribe) {die "\n\t Tribe undefined\n\n"; }
	push (@tribe_row, "Tribe");
	push (@tribe_row, "$virus_tribe");
	my $tribe_row = $writer->create_row(\@tribe_row); 
	push (@$html_ref, "\n$indent $tribe_row");

	my @vgenus_row;
	push (@vgenus_row, "Virus genus");
	push (@vgenus_row, "$genus");
	my $vgenus_row = $writer->create_row(\@vgenus_row); 
	push (@$html_ref, "\n$indent $vgenus_row");

	if ($subgroup eq 'Unclassified') { $subgroup = '-'; }
	my @subgroup_row;
	push (@subgroup_row, "Virus subgroup");
	push (@subgroup_row, "$subgroup");
	my $subgroup_row = $writer->create_row(\@subgroup_row); 
	push (@$html_ref, "\n$indent $subgroup_row");

	# Format the accession
	my @accession_row;
	unless ($accession eq 'NULL') {
		my $f_accession  = "<a href='http://www.ncbi.nlm.nih.gov/nuccore/";
		$f_accession .= "$accession' target='_blank'>$accession</a>";
		push (@accession_row, "Genbank accession");
		push (@accession_row, "$f_accession");
		my $accession_row = $writer->create_row(\@accession_row); 
		push (@$html_ref, "\n$indent $accession_row");
	}
	
	# Close table
	my $table_close = $writer->close_table($indent); 
	push (@$html_ref, $table_close);

	push (@$html_ref, "\n$indent<h3>Host taxonomy</h3>");
	push (@$html_ref, "\n$indent<div class=divider></div>");
	$table_open = $writer->open_table($indent, 'text'); # Open summary table
	push (@$html_ref, $table_open);

	if ($host_superclass) {
		my @host_superclass_row;
		push (@host_superclass_row, "Host superclass");
		push (@host_superclass_row, "$host_superclass");
		my $host_superclass_row = $writer->create_row(\@host_superclass_row); 
		push (@$html_ref, "\n$indent $host_superclass_row");
	}
	
	my @host_class_row;
	push (@host_class_row, "Host class");
	if ($host_class) {
		if ($host_class eq 'Aves') {
			$host_class .= ' (birds)';
		}
		elsif ($host_class eq 'Actinopterygii') {
			$host_class .= ' (ray-finned fishes)';
		}
		push (@host_class_row, "$host_class");
		my $host_class_row = $writer->create_row(\@host_class_row); 
		push (@$html_ref, "\n$indent $host_class_row");
	}
	else {
		#$devtools->print_hash($refseq);
	}
	
	if ($host_superorder) {
		my @host_superorder_row;
		push (@host_superorder_row, "Host superorder");
		if ($host_superorder eq 'Euarchontoglires') {
			$host_superorder .= ' (Supraprimates)';
		}
		push (@host_superorder_row, "$host_superorder");
		my $host_superorder_row = $writer->create_row(\@host_superorder_row); 
		push (@$html_ref, "\n$indent $host_superorder_row");
	}
	
	my @host_order_row;
	push (@host_order_row, "Host order");
	unless ($host_order) { $host_order = '-'; }
	push (@host_order_row, "$host_order");
	my $host_order_row = $writer->create_row(\@host_order_row); 
	push (@$html_ref, "\n$indent $host_order_row");

	if ($host_suborder) {	
		my @host_suborder_row;
		push (@host_suborder_row, "Host suborder");
		push (@host_suborder_row, "$host_suborder");
		my $host_suborder_row = $writer->create_row(\@host_suborder_row); 
		push (@$html_ref, "\n$indent $host_suborder_row");
	}

	my @host_family_row;
	push (@host_family_row, "Host family");
	push (@host_family_row, "$host_family");
	my $host_family_row = $writer->create_row(\@host_family_row); 
	push (@$html_ref, "\n$indent $host_family_row");

	# This is fairly redundant
	#my @host_genus_row;
	#push (@host_genus_row, "Host genus");
	#push (@host_genus_row, "$host_genus");
	#my $host_genus_row = $writer->create_row(\@host_genus_row); 
	#push (@$html_ref, "\n$indent $host_genus_row");
	
	my @host_species_row;
	push (@host_species_row, "Host species");
	$host_species =~ s/_/ /g;
	push (@host_species_row, "<i>$host_species</i> ($host_common)");
	my $host_species_row = $writer->create_row(\@host_species_row); 
	push (@$html_ref, "\n$indent $host_species_row");

	# Close table
	$table_close = $writer->close_table($indent); 
	push (@$html_ref, $table_close);
	#push (@$html_ref, "\n$indent<div class=divider></div>");

	if ($record_path) {
		my $f_fossil_record  = "Click <a href='$record_path'>here</a>";
		$f_fossil_record .= " to view the fossil record for this virus<br>";
		push (@$html_ref, $f_fossil_record);
	}

	# Link to raw reference sequence in GLUE format
	my $raw_link  = "\n$indent Click <a href='http://saturn.adarc.org/";
	   $raw_link .= "$site_name/db/refseq_flat/$name'>here</a>";
       $raw_link .= "\n$indent to view the raw reference for $name<br><br>";
	push (@$html_ref, $raw_link);
}

#***************************************************************************
# Subroutine:  write_feature_table 
# Description: 
#***************************************************************************
sub write_feature_table {

	my ($self, $refseq, $html_ref) = @_;

	# Set variables / get refseq data
	my $indent = ' ' x 10;
	my $name         = $refseq->{name};
	my $genes_ref    = $refseq->{genes};
	#$devtools->print_array($genes_ref); die;
	my $orf_hash_ref = $refseq->{translated_orfs};
	unless ($name ) { 
		die;  # We should have 
	}

	# Write genes table
	push (@$html_ref, "\n$indent<h3>Genes & sequence features</h3>");
	push (@$html_ref, "\n$indent<div class=divider></div>");
	my $table_open = $writer->open_table($indent, 'data'); # Open summary table
	push (@$html_ref, $table_open);

	#my @head_row_data = qw [ property value ];
	#my @amino_set1 = qw [ A G I L P V ];
	#my @amino_set2 = qw [ F Y W ];
	#my @amino_set3 = qw [ D E N Q S T ];
	#my @amino_set4 = qw [ C M ];
	#my @amino_set5 = qw [ R K H ];
	my @amino_set1 = qw [ K ];
	my @amino_set2 = qw [ R ];
	my @amino_set3 = qw [ H ];
	my @amino_set4 = qw [ P ];
	my @amino_set5 = qw [ C ];

	#my @head_row_data = qw [ Gene Coordinates Length Stops ];
	#my @head_row_data = qw [ Gene Coordinates Length AGILPV FYW DENQST R K H CM  * ? ];
	my @head_row_data = qw [ Gene Coordinates Length K R H P C * ];
	my $head_row = $writer->create_header_row(\@head_row_data); 
	push (@$html_ref, "\n$indent $head_row");
	foreach my $gene_ref (@$genes_ref) {
		
		#$devtools->print_hash($gene_ref);die;
		my @row;
		my $gene_name   = $gene_ref->{name};
		my $gene_start  = $gene_ref->{start};
		my $gene_stop   = $gene_ref->{stop};
		my $coordinates = $refseq->get_gene_coordinate_string($gene_name);
		unless ($coordinates) { die; }
		my $length = $gene_stop - $gene_start + 1;
		my $f_gene_name = lc $gene_name;
		push (@row, "<i>$f_gene_name</i>");
		push (@row, $coordinates);
		push (@row, $length);
		my $peptide_seq = $orf_hash_ref->{$gene_name};
		my %results;
		$seq_obj->get_aa_composition($peptide_seq, \%results);
		my %aminos_by_type;
		#$aminos_by_type{AGILPV}  =0;
		#$aminos_by_type{FYW}     =0;
		#$aminos_by_type{DENQST}  =0;
		#$aminos_by_type{CM}      =0;
		#$aminos_by_type{'?'}       =0;
		$aminos_by_type{K}       =0;
		$aminos_by_type{R}       =0;
		$aminos_by_type{H}       =0;
		$aminos_by_type{C}       =0;
		$aminos_by_type{P}       =0;
		foreach my $amino (@amino_set1) {
			my $amino_data_ref = $results{$amino};
			my $prop = $amino_data_ref->{f_percentage};
			#my $prop = $amino_data_ref->{f_proportion};
			#$aminos_by_type{AGILPV} = $aminos_by_type{AGILPV}  + $prop;
			$aminos_by_type{$amino} = $aminos_by_type{$amino} + $prop;
		}
		foreach my $amino (@amino_set2) {
			my $amino_data_ref = $results{$amino};
			my $prop = $amino_data_ref->{f_percentage};
			#my $prop = $amino_data_ref->{f_proportion};
			#$aminos_by_type{FYW} = $aminos_by_type{FYW} + $prop;
			$aminos_by_type{$amino} = $aminos_by_type{$amino} + $prop;
		}
		foreach my $amino (@amino_set3) {
			my $amino_data_ref = $results{$amino};
			my $prop = $amino_data_ref->{f_percentage};
			#my $prop = $amino_data_ref->{f_proportion};
			#$aminos_by_type{DENQST} = $aminos_by_type{DENQST} + $prop;
			$aminos_by_type{$amino} = $aminos_by_type{$amino} + $prop;
		}
		foreach my $amino (@amino_set4) {
			my $amino_data_ref = $results{$amino};
			my $prop = $amino_data_ref->{f_percentage};
			#my $prop = $amino_data_ref->{f_proportion};
			#$aminos_by_type{CM} = $aminos_by_type{CM} + $prop;
			$aminos_by_type{$amino} = $aminos_by_type{$amino} + $prop;
		}
		foreach my $amino (@amino_set5) {
			my $amino_data_ref = $results{$amino};
			my $prop = $amino_data_ref->{f_percentage};
			#my $prop = $amino_data_ref->{f_proportion};
			$aminos_by_type{$amino} = $aminos_by_type{$amino} + $prop;
		}
		#push (@row, $aminos_by_type{AGILPV});
		#push (@row, $aminos_by_type{FYW});
		#push (@row, $aminos_by_type{DENQST});
		push (@row, $aminos_by_type{K});
		push (@row, $aminos_by_type{R});
		push (@row, $aminos_by_type{H});
		push (@row, $aminos_by_type{P});
		push (@row, $aminos_by_type{C});
		#push (@row, $aminos_by_type{CM});

		my @peptide = split ('', $peptide_seq);
		pop @peptide;
		my $num_stops = 0;
		foreach my $char (@peptide) {
			if ($char eq '*') { $num_stops++; }
		}
		unless ($num_stops) { $num_stops = '-'; }
		push (@row, $num_stops);

		#my $missing_or_unrecognised = $aminos_by_type{'?'};
		#unless ($missing_or_unrecognised) { 
		#	$missing_or_unrecognised = '-';
		#}
		#push (@row, $missing_or_unrecognised);

        my $gene_row = $writer->create_row(\@row); 
		push (@$html_ref, "\n$indent $gene_row");

		# Get amino composition breakdown
		#push (@$html_ref, "$indent<br>ORF amino composition table<br>");
		#push (@$html_ref, "$indent<br><div class=divider></div>");
		#$self->make_aa_composition_table($refseq, $html_ref);
	}
	
	# Close table
	my $table_close = $writer->close_table($indent); 
	push (@$html_ref, $table_close);
	#push (@$html_ref, "\n$indent<div class=divider></div>");
}

#***************************************************************************
# Subroutine:  make_genome_feature_table
# Description: 1 value per virus
#***************************************************************************
sub make_genome_feature_table {

	my ($self, $refseq, $html_ref) = @_;

	my $seq_obj = Sequence->new();
	my $indent = ' ' x 10;
	
	# Open summary table
	my $table_open = $writer->open_table($indent, 'text');
	push (@$html_ref, $table_open);
	my @columns;
	push (@columns, "Feature name");
	push (@columns, "Description");
	push (@columns, "Other information");
	my $head_row = $writer->create_header_row(\@columns); 
	push (@$html_ref, "\n$indent $head_row");

	# TODO: fix this
	my $orf_hash_ref = $refseq->{translated_orfs};
	my $gag_orf = $orf_hash_ref->{gag};
	unless ($gag_orf) { $gag_orf = $orf_hash_ref->{Gag}; }
	unless ($gag_orf) { $gag_orf = $orf_hash_ref->{GAG}; }
	unless ($gag_orf) { $gag_orf = $orf_hash_ref->{NC}; }
	my $pol_orf = $orf_hash_ref->{pol};
	unless ($pol_orf) { $pol_orf = $orf_hash_ref->{Pol}; }
	unless ($pol_orf) { $pol_orf = $orf_hash_ref->{POL}; }
	unless ($pol_orf) { $pol_orf = $orf_hash_ref->{RT}; }
	my $pro_orf = $orf_hash_ref->{pro};
	my $separate_pro = 'true';
	unless ($pro_orf) { $pro_orf = $orf_hash_ref->{Pro}; }
	unless ($pro_orf) { $pro_orf = $orf_hash_ref->{PRO}; }
	unless ($pro_orf) { $pro_orf = $orf_hash_ref->{PR}; }
	unless ($pro_orf) { $pro_orf = $orf_hash_ref->{pr}; }
	unless ($pro_orf) { 
		$separate_pro = undef;
		$pro_orf = $pol_orf;
	}
	my $env_orf = $orf_hash_ref->{env};
	unless ($env_orf) { $env_orf = $orf_hash_ref->{Env}; }
	unless ($env_orf) { $env_orf = $orf_hash_ref->{ENV}; }
	unless ($env_orf) { $env_orf = $orf_hash_ref->{SU}; }

	my %coordinate_based;
	$self->get_coordinate_based_data($refseq, \%coordinate_based);

	#### LTR / leader region
	if ($coordinate_based{ltr_length}) {
		# Length of 5' leader (use coordinates)
		my @ltr_row;
		push (@ltr_row, "Long terminal repeat length");
		if ($coordinate_based{ltr_length}) {
			my $ltr_length = $coordinate_based{ltr_length};
			push (@ltr_row, $ltr_length);
			push (@ltr_row, '');
			my $row = $writer->create_row(\@ltr_row); 
			push (@$html_ref, "\n$indent $row");
		}
	}
	if ($coordinate_based{leader_length}) {

		# Length of 5' leader (use coordinates)
		my @leader_row;
		push (@leader_row, "Leader length");
		if ($coordinate_based{leader_length}) {
			my $lea_length = $coordinate_based{leader_length};
			push (@leader_row, $lea_length);
			push (@leader_row, '');
		}
		else {
			push (@leader_row, "-");
			push (@leader_row, '');
		}
		my $row = $writer->create_row(\@leader_row); 
		push (@$html_ref, "\n$indent $row");
	
		# PBS binding specificity  (use existing routine)
		my @pbs_row;
		push (@pbs_row, "Primer binding site (PBS)");
		#my $pbs = $self->get_pbs_specificity($refseq);
		push (@pbs_row, "-");
		push (@pbs_row, "");
		$row = $writer->create_row(\@pbs_row); 
		push (@$html_ref, "\n$indent $row");
	}

	# Gag
	if ($gag_orf) {
 		# Late domains
		$self->look_for_late_domain($gag_orf, $html_ref);

		# Zinc knuckle
		my @zinc_knuckle_row;
		push (@zinc_knuckle_row, "Zinc knuckles");
		my $zinc = $self->look_for_zinc_knuckle($gag_orf);
		if ($zinc) {
			push (@zinc_knuckle_row, $zinc); 
		}
		else { push (@zinc_knuckle_row, "none identified"); }
		push (@zinc_knuckle_row, '');
		my $row = $writer->create_row(\@zinc_knuckle_row); 
		push (@$html_ref, "\n$indent $row");

		# Major homology region (regex) D
		my @mhr_row;
		push (@mhr_row, "Major homology region");
		my $mhr = $self->look_for_major_homology_region($gag_orf);
		if ($mhr) {
			push (@mhr_row, 'yes'); 
			push (@mhr_row, $mhr); 
		}
		else { 
			push (@mhr_row, 'not identified'); 
			push (@mhr_row, ''); 
		}
		$row = $writer->create_row(\@mhr_row); 
		push (@$html_ref, "\n$indent $row");
	
		if ($pro_orf) {
			# Gag-pro/pol expression strategy (frameshift - use Gag and Pol coordinates)
			my @gagpro_row;
			push (@gagpro_row, "gag-pro junction");
			my $gagpro = $coordinate_based{gagpro};
			unless ($gagpro) { $gagpro = '-'; }	
			push (@gagpro_row, $gagpro);
			push (@gagpro_row, '');
			$row = $writer->create_row(\@gagpro_row); 
			push (@$html_ref, "\n$indent $row");
		}

		# Myristoylation (motif?)

	}

	if ($pro_orf) {

		# Protease active site (motif search) D
		my @pro_row;
		push (@pro_row, "Protease active site");
		my $pro_active = $self->look_for_pro_domain($pro_orf);
		if ($pro_active) {
			push (@pro_row, $pro_active); 
		}
		else { 
			push (@pro_row, '-'); 
		}
		push (@pro_row, ''); 
		my $row = $writer->create_row(\@pro_row); 
		push (@$html_ref, "\n$indent $row");
	
		# pro/pol expression strategy (frameshift - use Gag and Pol coordinates)
		my @propol_row;
		push (@propol_row, "pro-pol junction");
		my $propol;
		if ($separate_pro) {
			$propol = $coordinate_based{propol};
			unless ($propol) { $propol = '-'; }
			#unless ($propol) { die; }	
		}
		else {
			$propol = 'fused';
		}
		push (@propol_row, $propol);
		push (@propol_row, ''); 
		$row = $writer->create_row(\@propol_row); 
		push (@$html_ref, "\n$indent $row");
	}

	if ($pol_orf) {
		# RT active site (YMDD, YVDD, YIDD) (motif search) D
		my @rt_row;
		push (@rt_row, "RT active site");
		if ($pol_orf) {
			#print "\n\t looking for RT active";
			my $rt_active = $self->look_for_rt_domain($pol_orf);
			if ($rt_active) {
				push (@rt_row, $rt_active); 
			}
			else { 
				push (@rt_row, '-'); 
			}
			push (@rt_row, ''); 
			my $row = $writer->create_row(\@rt_row); 
			push (@$html_ref, "\n$indent $row");
		}
	}

	#### Envelope 
	if ($env_orf) {
		#	- immunosuppressive domain (regex?) 
		my @csk17_row;
		push (@csk17_row, "Env immunosuppressive region");
		if ($env_orf) {
			my $csk17_active = $self->look_for_csk17_domain($env_orf);
			if ($csk17_active) {
				push (@csk17_row, 'CSK17'); 
				push (@csk17_row, $csk17_active); 
			}
			else { 
				push (@csk17_row, '-'); 
				push (@csk17_row, ''); 
			}
			my $row = $writer->create_row(\@csk17_row); 
			push (@$html_ref, "\n$indent $row");
		}

		#	- 'env lineage' (BLAST against library)
		$self->get_env_type($refseq, $html_ref);

		#	- also a syncytin? (need to consult external file)  Rob find 
		#	- also a in its native species (need to consult external file) Rob find out
		# NLS REV-like (lysine & argine rich patches of sequence) (motif search (special case)
		$self->look_for_nls($refseq, $html_ref);
	}

	#### Accessory proteins
	# Close table
	my $table_close = $writer->close_table($indent); 
	push (@$html_ref, $table_close);
	#push (@$html_ref, "\n$indent<div class=divider></div>");
}

#***************************************************************************
# Subroutine:  make_dinuc_composition_table 
# Description: 
#***************************************************************************
sub make_dinuc_composition_table {

	my ($self, $refseq, $html_ref) = @_;

	# Write genes table
	my $seq_obj = Sequence->new();
	my $indent = ' ' x 10;
	my $table_open = $writer->open_table($indent, 'data'); # Open summary table
	push (@$html_ref, $table_open);
	#my @head_row_data = qw [ property value ];
	my @dinucs = qw [ AA AG AC AT GA GG GC GT CA CG CC CT TA TG TC TT ];
	my @header_row = @dinucs;
	unshift (@header_row, 'Region');
	my $head_row = $writer->create_header_row(\@header_row); 
	push (@$html_ref, "\n$indent $head_row");

	# Iterative through genes
	my %freqs;
	my %props;
	my $sequence = $refseq->{sequence};
	$seq_obj->get_dinucleotide_composition($sequence, \%freqs, \%props);
	my @row;
	unshift (@row, 'Genome');
	foreach my $dinuc (@dinucs) {
		my $prop = $props{$dinuc};
		my $f_prop = sprintf("%.2f", $prop);
		unless ($prop) { $f_prop = '-'}
		push (@row, $f_prop);
	}
	my $row = $writer->create_row(\@row); 
	push (@$html_ref, "\n$indent $row");
	
	# Iterative through genes
	my %genes;
	$refseq->get_orfs(\%genes);

	my $genes_ref = $refseq->{genes};
	foreach my $gene_ref (@$genes_ref) {
		#print "\n\t GENE $gene";	
		my $gene = $gene_ref->{name};
		my $gene_sequence = $genes{$gene};
		$seq_obj->get_dinucleotide_composition($gene_sequence, \%freqs, \%props);
		my @row2;
		my $f_gene_name = lc $gene;
		unshift (@row2, "<i>$f_gene_name</i>");
		foreach my $dinuc (@dinucs) {
			my $prop = $props{$dinuc};
			my $f_prop = sprintf("%.2f", $prop);
			unless ($prop) { $f_prop = '-'}
			push (@row2, $f_prop);
		}
		my $row = $writer->create_row(\@row2); 
		push (@$html_ref, "\n$indent $row");
	}

	# Close table
	my $table_close = $writer->close_table($indent); 
	push (@$html_ref, $table_close);
}

#***************************************************************************
# Subroutine:  write_gene_sequences_table 
# Description: 
#***************************************************************************
sub write_gene_sequences_table {

	my ($self, $refseq, $html_ref) = @_;

	# Set variables / get refseq data
	my $indent = ' ' x 10;
	my $name      = $refseq->{name};
	my $genes_ref = $refseq->{genes};
	my $sequence     = $refseq->{sequence};
	unless ($name ) { 
		die;  # We should have 
	}

	push (@$html_ref, "\n$indent<h4>Protein Sequences</h4>");
	push (@$html_ref, "\n$indent<div class=divider></div>");

	# Write translated ORFs
	my $orf_hash_ref = $refseq->{translated_orfs};
	unless ($orf_hash_ref) { die; }
	foreach my $gene_ref (@$genes_ref) {
		
		my $orf_name = $gene_ref->{name};
		my %orf_data;
		my $orf_seq = $orf_hash_ref->{$orf_name};
		my $gene_name = $name . '_' . $orf_name;
		my @data;
		my @seq;
		my @orf_seq = split ('',$orf_seq);
		my $f_seq;
		my $count= 0;
		foreach my $char (@orf_seq) {
			$count++;
			$f_seq .= $char;
			if ($count % 76 eq 0) {
				#$f_seq .= "\n$indent<br>";
			    $f_seq .= "\n";
			}
		}
		push (@$html_ref, "\n$indent <pre>");
		push (@$html_ref, "\n>$gene_name");
		push (@$html_ref, "\n$f_seq");
		push (@$html_ref, "\n$indent </pre>");
	}
	
	# Show sequence
	my @sequence = split('', $sequence);
	my $f_seq = '';
	my $count = 0;
	foreach my $char (@sequence) {
		$count++;
		$f_seq .= $char;
		if ($count % 76 eq 0) {
			#$f_seq .= "\n$indent";
			$f_seq .= "\n";
		}
	}
	#push (@$html_ref, "\n$indent <br><br>");
	push (@$html_ref, "\n$indent<h4>Genome Sequences</h4>");
	push (@$html_ref, "\n$indent<div class=divider></div>");
	push (@$html_ref, "\n$indent <pre>");
	#push (@$html_ref, "\n$indent $name<br>\n$indent $f_seq\n$indent<br><br>");
	#push (@$html_ref, "\n$sequence");
	push (@$html_ref, "\n>$name");
	push (@$html_ref, "\n$f_seq");
	push (@$html_ref, "\n$indent </pre>");

}

#***************************************************************************
# Subroutine:  setup_rv_fossil_stats 
# Description: 
#***************************************************************************
sub setup_rv_fossil_stats {

	my ($self, $stats_ref) = @_;

	# Create data structure 
	my @virus_names;
	my %refseqs;
	my %nt_orfs;
	my %aa_orfs;
	my %utrs;
	my @genomes;
	my %by_species;
	my %by_genus;
	my %by_host_order;
	my %by_host_superorder;
	my %by_subgroup;
	my %by_tribe;
	$stats_ref->{'by_species'}     = \%by_species;
	$stats_ref->{'by_genus'}       = \%by_genus;
	$stats_ref->{'by_subgroup'}    = \%by_subgroup;
	$stats_ref->{'aa_orfs'}        = \%aa_orfs;
	$stats_ref->{'nt_orfs'}        = \%nt_orfs;
	$stats_ref->{'utrs'}           = \%utrs;
	$stats_ref->{'genomes'}        = \@genomes;
	$stats_ref->{'viruses'}        = \@virus_names;
	$stats_ref->{'refseqs'}        = \%refseqs;
	$stats_ref->{'by_host_class'}  = \%by_species;
	$stats_ref->{'by_host_order'}  = \%by_host_order;
	$stats_ref->{'by_host_superorder'}  = \%by_host_superorder;
	$stats_ref->{'by_subgroup'}    = \%by_subgroup;
	$stats_ref->{'by_tribe'}       = \%by_tribe;
	
	# Set up fields 
	my @type_fields  = qw [ complete fragment ];
	my @state_fields = qw [ endogenous exogenous  ];
	# TO DO: need this?
}

#***************************************************************************
# Subroutine:  stats_by_virus_group
# Description:  
#***************************************************************************
sub stats_by_virus_group {

	my ($self, $refseq, $stats_ref) = @_;
	
	# Get refseq data
	my $refseq_name  = $refseq->{name};
	my $tribe        = $refseq->{virus_tribe};
	my $genus        = $refseq->{virus_genus};
	my $subgroup     = $refseq->{virus_subgroup};
	my $type         = $refseq->{sequence_type};
	my $state        = $refseq->{genome_state};
	my $coverage     = $refseq->{genome_coverage};
	unless ($type) { 
		print "\n\t\t.1 Refseq $refseq_name genome type undefined";
		$type = 'complete';
		$refseq->{genome_type} = $type;
	}
	unless ($state) { 
		print "\n\t\t.1 Refseq $refseq_name genome state undefined";
		$state = 'undefined';
	}
	#$genus =~ s/ //g;
	
    # Record by tribe
	my $stats_tribe    = $stats_ref->{'by_tribe'};
	unless ($stats_tribe) { die; }
	$self->stats_by_two_classifiers($type, $tribe, $stats_tribe);
	$self->stats_by_two_classifiers($state, $tribe, $stats_tribe);

	# Record by genus
	my $stats_genus = $stats_ref->{'by_genus'};
	unless ($stats_genus) { die; }
	$self->stats_by_two_classifiers($type, $genus, $stats_genus);
	$self->stats_by_two_classifiers($state, $genus, $stats_genus);
	
	# Record stats by genus subgroup and ge
	my $stats_subgroup = $stats_ref->{'by_subgroup'};
	#$self->stats_by_two_classifiers($type, $subgroup, $stats_subgroup);
	#$self->stats_by_two_classifiers($state, $subgroup, $stats_subgroup);
	#$devtools->print_hash($stats_ref); die;
}

#***************************************************************************
# Subroutine:  stats_by_host_group
# Description: statistics by host group
#***************************************************************************
sub stats_by_host_group {

	my ($self, $refseq, $stats_ref) = @_;

	# Get refseq data
	my $refseq_name     = $refseq->{name};
	my $virus_tribe     = $refseq->{virus_tribe};
	my $genus           = $refseq->{virus_genus};
	my $subgroup        = $refseq->{virus_subgroup};
	my $type            = $refseq->{sequence_type};
	my $state           = $refseq->{genome_state};
	my $coverage        = $refseq->{genome_coverage};
	my $host_order      = $refseq->{host_order};
	my $host_superorder = $refseq->{host_superorder};
	unless ($type) { 
		print "\n\t\t.1 Refseq $refseq_name genome type undefined";
		$type = 'complete';
		$refseq->{genome_type} = $type;
	}
	unless ($state) { 
		print "\n\t\t.1 Refseq $refseq_name genome state undefined";
		$state = 'undefined';
	}

    # Record by host species
	my $stats_species    = $stats_ref->{'by_species'};
	#$self->stats_by_two_classifiers($virus_class, $host_species, $stats_species);
	#$self->stats_by_two_classifiers($virus_genus, $host_species, $stats_species);

    # Record by host superorder
	my $stats = $stats_ref->{'by_host_superorder'};
	unless ($stats) { die; }
	if ($host_superorder) {
		$self->stats_by_two_classifiers($virus_tribe, $host_superorder, $stats);
	}
	#$devtools->print_hash($stats);die;

    # Record by host superorder
	$stats = $stats_ref->{'by_host_order'};
	unless ($stats) { die; }
	if ($host_order) {
		$self->stats_by_two_classifiers($virus_tribe, $host_order, $stats);
	}
}

#***************************************************************************
# Subroutine:  stats_by_two_classifiers 
# Description:  
#***************************************************************************
sub stats_by_two_classifiers {

	my ($self, $key, $stratifier, $stats_ref) = @_;

	#print "\n\t Class key $key ($stratifier)";
	if ($stats_ref->{$stratifier}) {
		my $strat_ref = $stats_ref->{$stratifier};
		if ($strat_ref->{$key}) {
			$strat_ref->{$key}++;
		}
		else {
			$strat_ref->{$key} = 1;
		}
	}
	else {
		my %stratifier;	
		my %stats;
		$stats{$key} = 1;
		$stats_ref->{$stratifier} = \%stats;
	}
}



############################################################################
# SECTION: WRITING RETROVIRUS REFSEQ LIBRARY STATS PAGE 
############################################################################

#***************************************************************************
# Subroutine:  total_ervs_by_tribe	
# Description: 
#***************************************************************************
sub total_ervs_by_tribe {

	my ($self, $table_ref, $stats_ref) = @_;
	
	# Set up fields 
	my $stats_tribe = $stats_ref->{'by_tribe'};
	my @tribes = sort keys %$stats_tribe;
	my %totals;
	foreach my $tribe (@tribes) { # Initialise
		$totals{$tribe} = 0;
	}
	#my @type_fields  = qw [ complete fragment ];
	my @type_fields  = qw [ locus consensus ancestral ];

	# Open summary table
	push (@$table_ref, "\n$indent <br><h3> Retroviruses by tribe</h3>");
	my $table_open = $writer->open_table($indent, 'text'); 
	push (@$table_ref, $table_open);
	my $column_tags =  "<col width=\"40\" />
                        <col width=\"15\" />
                        <col width=\"15\" />
                        <col width=\"15\" />
                        <col width=\"15\" />";
	push (@$table_ref, $column_tags);
	my @columns  = qw [ Tribe total locus consensus ancestral ];
	my $head_row = $writer->create_header_row(\@columns); 
	push (@$table_ref, "\n$indent $head_row"); 

	# Stats by tribe	
	foreach my $tribe (@tribes) { 
		my @row;
		my $tribe_stats = $stats_tribe->{$tribe};
		#$devtools->print_hash($class_stats);
		my $tribe_total = 0;
		foreach my $field (@type_fields) {
			my $value =	$tribe_stats->{$field};
			unless ($value) { $value = '0'; }
			print "\n\t Tribe $tribe: $field $value";
			$tribe_total = $tribe_total + $value;
			$totals{$tribe} = $totals{$tribe} + $value;
			$totals{$field} = $totals{$field} + $value;
			unless ($value) { $value = '-'; }
			push (@row, "$value");
		}
		unshift (@row, $tribe_total);
		unshift (@row, $tribe);
		my $row = $writer->create_row(\@row); 
		push (@$table_ref, "\n$indent $row"); 
	}
	my @total_row;
	my $sum_total;
	foreach my $field (@type_fields) {
		my $total = $totals{$field};
		push (@total_row, $total);
		$sum_total = $sum_total + $total;
	}
	unshift (@total_row, "<u>$sum_total</u>");
	unshift (@total_row, "Totals:");
	my $row = $writer->create_row(\@total_row); 
	push (@$table_ref, "\n$indent $row"); 

	# Close table
	my $table_close = $writer->close_table($indent); 
	push (@$table_ref, $table_close);
}

#***************************************************************************
# Subroutine:  total_ervs_by_genus	
# Description: 
#***************************************************************************
sub total_ervs_by_genus {

	my ($self, $table_ref, $stats_ref) = @_;
	
	# Open summary table
	push (@$table_ref, "\n$indent <br><h3> Retroviruses by genus</h3>");
	#my @type_fields  = qw [ complete fragment ];
	my @type_fields  = qw [ locus consensus ancestral ];
	my $table_open = $writer->open_table($indent, 'text'); 
	push (@$table_ref, $table_open);
	my $column_tags =  "<col width=\"40\" />
                        <col width=\"15\" />
                        <col width=\"15\" />
                        <col width=\"15\" />
                        <col width=\"15\" />";
	push (@$table_ref, $column_tags);
	my @columns  = qw [ Genus total locus consensus ancestral ];
	my $head_row = $writer->create_header_row(\@columns); 
	push (@$table_ref, "\n$indent $head_row"); 

	# Set up totals
	my %totals;
	foreach my $field (@type_fields) { $totals{$field} = 0; }

	# Stats by genus	
	my $stats_genus  = $stats_ref->{'by_genus'};
	unless ($stats_genus)  { die; }
	my @sort_genera = sort keys %$stats_genus;
	foreach my $genus (@sort_genera) { 
		my @row;
		my $genus_stats = $stats_genus->{$genus};
		#$devtools->print_hash($class_stats);
		my $total = 0;
		foreach my $field (@type_fields) {
			my $value =	$genus_stats->{$field};
			unless ($value) { $value = '0'; }
			$total = $total + $value;
			print "\n\t Genus $genus: $field $value";
			$totals{$field} = $totals{$field} + $value;
			unless ($value) { $value = '-'; }
			push (@row, $value);
		}
		unshift (@row, $total);
		unshift (@row, "$genus");
		my $row = $writer->create_row(\@row); 
		push (@$table_ref, $row); 
	}

	# Total table
	my @total_row;
	my $sum_total;
	foreach my $field (@type_fields) {
		my $total = $totals{$field};
		push (@total_row, $total);
		$sum_total = $sum_total + $total;
	}
	unshift (@total_row, "<u>$sum_total</u>");
	unshift (@total_row, "Totals:");
	my $row = $writer->create_row(\@total_row); 
	push (@$table_ref, "\n$indent $row");
 
	# Close table
	my $table_close = $writer->close_table($indent); 
	push (@$table_ref, $table_close);
}

#***************************************************************************
# Subroutine:  rvs_by_host_taxonomy
# Description: 
#***************************************************************************
sub rvs_by_host_taxonomy {

	my ($self, $table_ref, $stats_ref, $superorder_ref) = @_;

	# Stats by genus	
	#$devtools->print_hash($stats_ref); die;
	my $stats_tribe       = $stats_ref->{'by_tribe'};
	my $stats_order       = $stats_ref->{'by_host_order'};
	my $stats_superorder  = $stats_ref->{'by_host_superorder'};
	unless ($stats_tribe)  { die; }

	# Open summary table
	my $table_open = $writer->open_table($indent, 'text'); 
	push (@$table_ref, $table_open);

	my $column_tags =  "<col width=\"30\" />
                        <col width=\"10\" />
                        <col width=\"10\" />
                        <col width=\"10\" />
                        <col width=\"10\" />
                        <col width=\"10\" />
                        <col width=\"10\" />
                        <col width=\"10\" />";
	push (@$table_ref, $column_tags);

	my @tribes = sort keys %$stats_tribe;
	my @head_row_data = @tribes;
	unshift (@head_row_data, "Superorder");
	my $head_row = $writer->create_header_row(\@head_row_data); 
	push (@$table_ref, $head_row);

	# Create genomes table
	my $indent = ' ' x 10;
	my @names = sort keys %$superorder_ref;
	foreach my $name (@names) {

		# Create row
		my @row;
		push (@row, $name);
		my $superorder_ref = $stats_superorder->{$name};
		foreach my $tribe (@tribes) {    
			#$devtools->print_hash($by_superorder);
			my $value = $superorder_ref->{$tribe};
			unless ($name and $tribe and $value) {
				print "\n\t ## Superorder $name, tribe '$tribe': '$value'";
			}
			else {
				print "\n\t ## No data for superorder '$name'";
			}
			unless ($value)  { $value = '-'; }
			push (@row, $value);
			#my $value = $tribe_stats->{$field};
			#foreach my $field (@type_fields) {
			#}
		}
		my $row = $writer->create_row(\@row); 
		push (@$table_ref, $row); 

	}
	# Close table
	my $table_close = $writer->close_table($indent); 
	push (@$table_ref, $table_close);
}

#***************************************************************************
# Subroutine:  rvs_by_host_taxonomy
# Description: 
#***************************************************************************
sub rvs_by_host_taxonomy2 {

	my ($self, $table_ref, $stats_ref, $order_ref) = @_;

	# Stats by genus	
	#$devtools->print_hash($stats_ref); die;
	my $stats_tribe       = $stats_ref->{'by_tribe'};
	my $stats_order       = $stats_ref->{'by_host_order'};
	my $stats_superorder  = $stats_ref->{'by_host_superorder'};
	unless ($stats_tribe)  { die; }

	# Open summary table
	my $table_open = $writer->open_table($indent, 'text'); 
	push (@$table_ref, $table_open);

	my $column_tags =  "<col width=\"30\" />
                        <col width=\"10\" />
                        <col width=\"10\" />
                        <col width=\"10\" />
                        <col width=\"10\" />
                        <col width=\"10\" />
                        <col width=\"10\" />
                        <col width=\"10\" />";
	push (@$table_ref, $column_tags);



	my @tribes = sort keys %$stats_tribe;
	my @head_row_data = @tribes;
	unshift (@head_row_data, "Order");
	my $head_row = $writer->create_header_row(\@head_row_data); 
	push (@$table_ref, $head_row);

	# Create genomes table
	my $indent = ' ' x 10;
	my @names = sort keys %$order_ref;
	foreach my $name (@names) {

		print "\n\t ## Order '$name'";
		# Create row
		my @row;
		push (@row, $name);
		my $order_ref = $stats_order->{$name};
		foreach my $tribe (@tribes) {    
			#$devtools->print_hash($by_superorder);
			my $value = $order_ref->{$tribe};
			unless ($value)  { $value = '-'; }
			print "\n\t ## Order $name, tribe '$tribe': '$value'";
			push (@row, $value);
			#my $value = $tribe_stats->{$field};
			#foreach my $field (@type_fields) {
			#}
		}
		my $row = $writer->create_row(\@row); 
		push (@$table_ref, $row); 

	}
	# Close table
	my $table_close = $writer->close_table($indent); 
	push (@$table_ref, $table_close);

}

############################################################################
# SECTION: look for retrovirus specific features
############################################################################

#***************************************************************************
# Subroutine: 	get_coordinate_based_data
# Description:  
#***************************************************************************
sub get_coordinate_based_data {

	my ($self, $refseq, $data_ref) = @_;

	# Get UTRs 	
	my %utrs;
	$refseq->get_utrs(\%utrs);

	my $ltr = $utrs{LTR};
	if ($ltr) {	
		my $ltr_len = length $ltr;
		$data_ref->{ltr_length} = $ltr_len;

	}
	my $lea = $utrs{LEA};
	if ($lea) {	
		my $lea_len = length $lea;
		$data_ref->{leader_length} = $lea_len;
	}

	my $genes_ref = $refseq->{genes};
	my %frames;
	foreach my $gene_ref (@$genes_ref) {

		my $start = $gene_ref->{start};
		my $name = $gene_ref->{name};
		my $remainder = $start % 3;
		my $frame;
		if ($remainder eq 1) {
			$frame = 2;
		}
		elsif ($remainder eq 2) {
			$frame = 3;
		}
		else {
			$frame = 1;
		}
		my $f_name = lc $name;
		#print "\n\t\t\t ~~ gene $f_name is in frame $frame";
		$frames{$f_name} = $frame;
	}
	my $gag_frame = $frames{gag};
	my $pol_frame = $frames{pol};
	my $pro_frame = $frames{pro};
	

	# GAG-PRO
	if ($gag_frame and $pro_frame) {
		if ($gag_frame eq $pro_frame) {
			$data_ref->{gagpro} = 'readthrough';
		}
		elsif ($gag_frame < $pro_frame) {
			if ($pro_frame - $gag_frame eq 1) {
				$data_ref->{gagpro} = '-1 FS';
			}
			elsif ($pro_frame - $gag_frame eq 2) {
				$data_ref->{gagpro} = '-2 FS';
			}
		}
		elsif ($pro_frame < $gag_frame) {
			if ($gag_frame - $pro_frame eq 1) {
				$data_ref->{gagpro} = '-1 FS';
			}
			elsif ($gag_frame - $pro_frame eq 2) {
				$data_ref->{gagpro} = '-2 FS';
			}

		}
		else { die; }
	}
	# GAG-POL
	elsif ($gag_frame and $pol_frame) {
		if ($gag_frame eq $pol_frame) {
			$data_ref->{gagpro} = 'readthrough';
		}
		elsif ($gag_frame < $pol_frame) {
			if ($pol_frame - $gag_frame eq 1) {
				$data_ref->{gagpro} = '-1 FS';
			}
			elsif ($pol_frame - $gag_frame eq 2) {
				$data_ref->{gagpro} = '-2 FS';
			}
		}
		elsif ($pol_frame < $gag_frame) {
			if ($gag_frame - $pol_frame eq 1) {
				$data_ref->{gagpro} = '-1 FS';
			}
			elsif ($gag_frame - $pol_frame eq 2) {
				$data_ref->{gagpro} = '-2 FS';
			}
		}
	}

	# PRO POL
	if ($pol_frame and $pro_frame) {
		if ($pol_frame eq $pro_frame) {
			$data_ref->{propol} = '';
		}
		elsif ($pro_frame < $pol_frame) {
			if ($pol_frame - $pro_frame eq 1) {
				$data_ref->{propol} = '-1 FS';
			}
			elsif ($pol_frame - $pro_frame eq 2) {
				$data_ref->{propol} = '-2 FS';
			}
		}
		elsif ($pol_frame < $pro_frame) {
			if ($pro_frame - $pol_frame eq 1) {
				$data_ref->{propol} = '-1 FS';
			}
			elsif ($pro_frame - $pol_frame eq 2) {
				$data_ref->{propol} = '-2 FS';
			}
		}
	}

}

#***************************************************************************
# Subroutine: 	look_for_late_domain
# Description:  
#***************************************************************************
sub look_for_late_domain {

	my ($self, $orf, $html_ref) = @_;
	
	# Make a file for BLAST
	my $indent = ' ' x 10;

	my $ptap = 'P[TS]AP';
	my $ppxy = 'PP\wY';
	my $ypxl = 'YP\wL';

	# Late domains 
	my $number_ptap =()= $orf =~ /$ptap/gi;
	my $number_ppxy =()= $orf =~ /$ppxy/gi;
	my $number_ypxl =()= $orf =~ /$ypxl/gi;
	my @details;
	my $found = undef;
	if ($number_ptap) {
		push (@details, "PTAP ($number_ptap)");
		$found = "yes";
	}
	if ($number_ppxy) {
		push (@details, "PPxY ($number_ppxy)");
		$found = "yes";
	}
	if ($number_ypxl) {
		push (@details, "YPxL ($number_ypxl)");
		$found = "yes";
	}

	# Details
	my @late_domain_row;
	push (@late_domain_row, "Gag late domain");
	my $details = join(',', @details);	
	if ($found) {
		push (@late_domain_row, $found);
		push (@late_domain_row, $details);
	}
	else {
		push (@late_domain_row, 'no');
		push (@late_domain_row, '');
	}
	my $row = $writer->create_row(\@late_domain_row); 
	push (@$html_ref, "\n$indent $row");

}

#***************************************************************************
# Subroutine: 	look_for_major_homology_region
# Description:  
#***************************************************************************
sub look_for_major_homology_region {

	my ($self, $orf) = @_;
	my $mhr_regex = '(I|V)\wQG\w\wE(P|S|D)(F|Y|P)\w\w(F|Y)\w\wR(L|F)';
	my $match = $self->find_peptide_motif($orf, $mhr_regex);
	if ($match) {
		return $match;
	}
	# Other information
}

#***************************************************************************
# Subroutine: 	look_for_zin_knuckle
# Description:  
#***************************************************************************
sub look_for_zinc_knuckle {

	my ($self, $orf) = @_;
	my $zinc_regex = 'C\w\wC\w\w\w\wH\w\w\w\wC';
	#my $match = $self->find_peptide_motif($orf, $zinc_regex);
	my $number =()= $orf =~ /$zinc_regex/gi;
	return $number;
}

#***************************************************************************
# Subroutine: 	look_for_pro_domain
# Description:  
#***************************************************************************
sub look_for_pro_domain {

	my ($self, $orf) = @_;
	my $pro_regex = 'DTGA';
	if ($self->find_peptide_motif($orf, $pro_regex)) {
		return 'D<u>T</u>GA';
	}
	$pro_regex = 'DSGA';
	if ($self->find_peptide_motif($orf, $pro_regex)) {
		return 'D<u>S</u>GA';
	}
	$pro_regex = 'DLGV';
	if ($self->find_peptide_motif($orf, $pro_regex)) {
		return 'D<u>L</u>G<u>V</u>';
	}
	$pro_regex = 'DTGS';
	if ($self->find_peptide_motif($orf, $pro_regex)) {
		return $pro_regex;
	}
	$pro_regex = 'DVGA';
	if ($self->find_peptide_motif($orf, $pro_regex)) {
		return 'D<u>V</u>GA';
	}
}

#***************************************************************************
# Subroutine: 	look_for_rt_domain
# Description:  
#***************************************************************************
sub look_for_rt_domain {

	my ($self, $orf) = @_;

	my $rt_regex = 'YMDD';
	if ($self->find_peptide_motif($orf, $rt_regex)) {
		return 'Y<u>M</u>DD';
	}
	$rt_regex = 'YVDD';
	if ($self->find_peptide_motif($orf, $rt_regex)) {
		return 'Y<u>V</u>DD';
	}
	$rt_regex = 'YIDD';
	if ($self->find_peptide_motif($orf, $rt_regex)) {
		return 'Y<u>I</u>DD';
	}
	# Other information
}

#***************************************************************************
# Subroutine: 	look_for_csk17_domain
# Description:  
#***************************************************************************
sub look_for_csk17_domain {

	my ($self, $orf) = @_;
	my $csk17_regex = '(G|A|D|M)(L|P|I)(D|H|N)(\w{2,5})(L|F|P|S|V|T|I)(\w{10,14})(G|D|K|E|S)(G|E|S|R|H)';
	#my $csk17_regex = '(L|Y|F|P|W|A)(Q|N|E)N(\w{6,6})(G|A|D|M)(L|P|I)(D|H|N)(\w{2,5})(L|F|P|S|V|T|I)(\w{10,14})(G|D|K|E|S)(G|E|S|R|H)';
	my $match = $self->find_peptide_motif($orf, $csk17_regex);
	if ($match) {
		return $match;
	}
	# Other information
}

#***************************************************************************
# Subroutine: 	get_env_type
# Description:  
#***************************************************************************
sub get_env_type {

	my ($self, $refseq, $html_ref) = @_;
	
	# Get the env gene
	my $orf_hash_ref = $refseq->{translated_orfs};
	my $env_orf = $orf_hash_ref->{env};
	unless ($env_orf) { $env_orf = $orf_hash_ref->{Env}; }
	unless ($env_orf) { $env_orf = $orf_hash_ref->{ENV}; }
	unless ($env_orf) { $env_orf = $orf_hash_ref->{SU}; }

	# Make a file for BLAST
	my $name = $refseq->{name};
	my $fasta      = ">$name env\n$env_orf";
	my $result_path = './site/tmp/'; # TODO remove hard-coded path
	my $query_file = $result_path . $name . '.fas';
	$fileio->write_text_to_file($query_file, $fasta);
	my $result_file = $result_path . $name . '.blast_result';
		
	# BLASTx against the env library
	my $lib_path = 'db/blast/envelope_groups.fasta';
	my $blast_alg = 'blastp';
	#print "\n\t ##### ENV BLAST for $name";
	my $blast_obj = $self->{blast_obj};
	unless ($blast_obj) { die; }
	$blast_obj->blast($blast_alg, $lib_path, $query_file, $result_file);
	my @results;
	$blast_obj->parse_tab_format_step_one($result_file, \@results);
	# Get the best match from this file
	my $top_match = shift @results;
	my $env = $top_match->{scaffold};
	my $bitscore = $top_match->{bit_score};
	#$devtools->print_hash($top_match); die;
	my $env_type1;
	my $env_type2;
	unless ($env) {
		$env_type1 = '-';
		$env_type2 = '-';
	}
	else {
		my @env = split(/\|/, $env);
		$env_type1 = $env[0];
		$env_type2 = $env[1];
		$env_type2 .= " ($bitscore)";
	}

	my $indent = ' ' x 10;
	my @env_type_row;
	push (@env_type_row, "<i>env</i> type");
	push (@env_type_row, $env_type1);
	push (@env_type_row, $env_type2);
	my $row = $writer->create_row(\@env_type_row); 
	push (@$html_ref, "\n$indent $row");
}

#***************************************************************************
# Subroutine: 	look_for_nls
# Description:  look for nuclear localization signal 
#***************************************************************************
sub look_for_nls {

	my ($self, $refseq, $html_ref) = @_;
	
	# Get the env gene
	my $orf_hash_ref = $refseq->{translated_orfs};
	my $env_orf = $orf_hash_ref->{env};
	unless ($env_orf) { $env_orf = $orf_hash_ref->{Env}; }
	unless ($env_orf) { $env_orf = $orf_hash_ref->{ENV}; }
	unless ($env_orf) { $env_orf = $orf_hash_ref->{SU}; }
	
	my @nls_row;
	push (@nls_row, "Env NLS (predicted)");
	my $nls = '[KR][KR][KR]';
	my $number =()= $env_orf =~ /$nls/gi;
	unless ($number) {
		$number = 'no';
	}	
	push (@nls_row, $number);
	push (@nls_row, '');
	my $row = $writer->create_row(\@nls_row); 
	my $indent = ' ' x 10;
	push (@$html_ref, "\n$indent $row");

}

#***************************************************************************
# Subroutine: 	
# Description: 
#***************************************************************************
sub find_peptide_motif {

	my ($self, $orf, $motif) = @_;
	if ($orf =~ /($motif)/) {
		return $1;
	}
}


#***************************************************************************
# Subroutine: 	
# Description: 
#***************************************************************************
sub merge_libraries {

	my ($self, $lib_path1, $lib_path2) = @_;

	my $parser_obj    = RefSeqParser->new();
	my @files1;
	my @files2;
	$fileio->read_directory_to_array($lib_path1, \@files1);	
	$fileio->read_directory_to_array($lib_path2, \@files2);	
	my $outdir = 'merged/';

	foreach my $file (@files1) {

		print "\n\t DOING $file";

		# Create reference sequence object
		my %params1;
		my $path1 = $lib_path1 . $file;
		$parser_obj->parse_refseq_flatfile($path1, \%params1);
		my $refseq1 = RefSeq->new(\%params1);

		# Create reference sequence object
		my %params2;
		my $path2 = $lib_path2 . $file;
		my $exists = $fileio->check_file_exists($path2);
		unless ($exists) {
			print "\n\t    SKIPPING $file";
			next;
		}

		$parser_obj->parse_refseq_flatfile($path2, \%params2);
		my $refseq2 = RefSeq->new(\%params2);
	
		$refseq1->{features} = $refseq2->{features};
		my $tribe = $refseq1->{virus_tribe};
		my $full_path = $outdir . $tribe . '/';
		$refseq1->write_self_to_text($full_path);
	} 
}


############################################################################
# CREATING RV BLAST LIBRARIES FROM REFERENCE LIBRARY
############################################################################

#***************************************************************************
# Subroutine:  create_blast_libraries
# Description: 
#***************************************************************************
sub create_blast_libraries {

	my ($self) = @_;

	# Load the retroviral reference sequences
	my %refseqs;
	my @refseqs;
	my $refseqlib_obj = RefSeqLibrary->new();
	my $libpath = $self->{refseq_lib_path};
	print "\n\t loading GLUE reference sequence library ' $libpath'";
	$refseqlib_obj->load_refseq_library($libpath, \@refseqs, \%refseqs);	

	# Get the virus data
	my %by_taxonomy;
	my %flat;
	my @ranks = qw [ virus_family virus_subfamily virus_supertribe 
                     virus_tribe virus_genus virus_subgroup name ];
	$refseqlib_obj->order_by_taxonomy(\@refseqs, \@ranks, \%by_taxonomy, \%flat);
	#$devtools->print_hash(\%by_taxonomy); die;

	# Write out reference sequence group and summary files
	my @all_genomes;
	my @all_leaders;
	my @all_ltrs;
	my @all_orfs_aa;
	my @all_orfs_nt;
	my @all_envs_aa;
	my $blast_db_path = $self->{blast_db_path};
	
	# Iterate through the viruses 
	my @viruses = sort keys %flat;
	foreach my $virus (@viruses) {
	
		my $refseq = $flat{$virus};
		print "\n\t Adding data for '$virus'";
		
		# Get complete genome
		print "\n\t\t Adding '$virus' genome";
		my $genome = $refseq->{sequence};
		my $g_fasta = ">$virus\n$genome\n";
		push (@all_genomes, $g_fasta);
		
		# Get UTRs 	
		my %utrs;
		$refseq->get_utrs(\%utrs);

		my $ltr = $utrs{LTR};
		if ($ltr) {	
			print "\n\t\t Adding '$virus' LTR";
			my $ltr_fasta = ">$virus" . "_ltr\n$ltr\n";
			push (@all_ltrs, $ltr_fasta);
		}
		my $lea = $utrs{LEA};
		if ($lea) {	
			print "\n\t\t Adding '$virus' LEA";
			my $lea_fasta = ">$virus" . "_lea\n$lea\n";
			push (@all_leaders, $lea_fasta);
		}

		# Get ORFs	
		my %genes;
		$refseq->get_orfs(\%genes);
		my @genes = keys %genes;
		my %orfs;
		$refseq->get_translated_orfs(\%orfs);
		foreach my $gene_name (@genes) {
			#$devtools->print_hash($gene_ref);
			print "\n\t\t Adding '$virus' gene '$gene_name'";
			my $gene_seq = $genes{$gene_name};
			unless ($gene_name and $gene_seq) { die; }
			my $nt_fasta = ">$virus" . "_$gene_name\n$gene_seq\n";
			push (@all_orfs_nt, $nt_fasta);

			my $protein_seq = $orfs{$gene_name};
			unless ($protein_seq) { die; }
			my $aa_fasta = ">$virus" . "_$gene_name\n$protein_seq\n";
			push (@all_orfs_aa, $aa_fasta);
			
			# Get envs only
			my $lc_gene_name = lc $gene_name;
			if ($lc_gene_name eq 'env') {
				push (@all_envs_aa, $aa_fasta);
			}
		}
	}
	#$devtools->print_array(\@all_orfs_nt);die;
	
	# Write entire library ORFs as NT
	my $blast_nt_orf_lib_path = $self->{blast_nt_orf_lib_path};
	$fileio->write_output_file($blast_nt_orf_lib_path, \@all_orfs_nt);

	# Format for BLAST
	my $formatdb_command = "./bin/blast/makeblastdb -in $blast_nt_orf_lib_path -dbtype nucl > /dev/null";
	print "\n\t #  $formatdb_command";
	system $formatdb_command;
	
	# Write entire library ORFs as AA
	my $blast_orf_lib_path = $self->{blast_orf_lib_path};
	$fileio->write_output_file($blast_orf_lib_path, \@all_orfs_aa);
	# Format for BLAST
	$formatdb_command = "./bin/blast/makeblastdb -in $blast_orf_lib_path > /dev/null";
	print "\n\t #  $formatdb_command";
	system $formatdb_command;

	# Write env ORFs as AA
	my $blast_env_lib_path = $self->{blast_env_lib_path};
	unless ($blast_env_lib_path) { die; }
	$fileio->write_output_file($blast_env_lib_path, \@all_envs_aa);
	# Format for BLAST
	$formatdb_command = "./bin/blast/makeblastdb -in $blast_env_lib_path > /dev/null";
	print "\n\t #  $formatdb_command";
	system $formatdb_command;

	
	# Write entire library UTRs
	my @all_utrs = ( @all_ltrs, @all_leaders);
	my $blast_utr_lib_path = $self->{blast_utr_lib_path};
	$fileio->write_output_file($blast_utr_lib_path, \@all_utrs);
	
	# Format for BLAST
	$formatdb_command = "./bin/blast/makeblastdb -in $blast_utr_lib_path -dbtype nucl > /dev/null";
	print "\n\t #  $formatdb_command";
	system $formatdb_command;
	
	# Write entire library complete genomes
	my $blast_genome_lib_path = $self->{blast_genome_lib_path};
	$fileio->write_output_file($blast_genome_lib_path, \@all_genomes);
	# Format for BLAST
	$formatdb_command = "./bin/blast/makeblastdb -in $blast_genome_lib_path -dbtype nucl > /dev/null";
	print "\n\t #  $formatdb_command";
	system $formatdb_command;
}

#***************************************************************************
# Subroutine:  validate refseq
# Description: validate refseq
#***************************************************************************
sub validate_retrovirus_refseq {

	my ($self, $refseq, $stats_ref) = @_;

	# Get refseq data
	my $name        = $refseq->{name};
	my $full_name   = $refseq->{full_name};
	my $type        = $refseq->{genome_type};
	my $state       = $refseq->{virus_state};
	my $group       = $refseq->{virus_group};
	my $subgroup    = $refseq->{virus_subgroup};
	my $genus       = $refseq->{virus_group};
	my $class       = $refseq->{virus_class};
	my $host_family = $refseq->{host_family};
	my $host_order  = $refseq->{host_order};
	my $host_sci_name    = $refseq->{host_sci_name};
	my $host_superorder  = $refseq->{host_superorder};
	my $host_common_name = $refseq->{host_common_name};
	$class    =~ s/\s+//g; # Remove whitespace
	$group    =~ s/\s+//g; # Remove whitespace
	$subgroup =~ s/\s+//g; # Remove whitespace
	$name     =~ s/\s+//g; # Remove whitespace

	# Subfamily
	my $subfamily;
	if ($genus eq 'Spumavirus') { 
		$subfamily = "Spumaretrovirinae";	
	}
	else {
		$subfamily = "Orthoretrovirinae";	
	}

	# Supertribe
	my $supertribe;
    if    ($class eq 'I')   { $supertribe = 'Alnilamretrovirinae'; }
    elsif ($class eq 'III') { 
			$supertribe = 'Alnitakretrovirinae'; 
			$subfamily = "Spumaretrovirinae";	
	}
    elsif ($class eq 'II')  { $supertribe = 'Mintakaretrovirinae'; }
	else { 
	     my $class = $refseq->{virus_class};
         $devtools->print_hash($refseq);
	     die "\n\t no class for $name get '$class'"; 
    }
	$refseq->{virus_supertribe} = $supertribe;
	$refseq->{virus_subfamily} = $subfamily;


	# Tribe
	# Lenti/Spuma/Class I-II-III
	my $tribe = 'Unclassified';
	my $new_genus = 'Unclassified';
	my $new_subgroup;
	if ($genus) {
		if ($genus eq 'Alpharetrovirus') {
			$tribe = 'Alpharetrovirus';
			$new_genus = 'Alpharetrovirus';
		}
		elsif ($genus eq 'Betaretrovirus') {
			$tribe = 'Mintaka';
			$new_genus = 'Betaretrovirus';
		}
		elsif ($genus eq 'Gammaretrovirus') {
			$tribe = 'Alnilam';
			$new_genus = 'Gammaretrovirus';
		}
		elsif ($genus eq 'Deltaretrovirus') {
			$tribe = 'Delta';
			$new_genus = 'Deltaretrovirus';
		}
		elsif ($genus eq 'Epsilonretrovirus') {
			$tribe = 'Epsilon';
			$new_genus = 'Epsilonretrovirus';
		}
		elsif ($genus eq 'Lentivirus') {
			$tribe = 'Lenti';
			$new_genus = 'Lentivirus';
		}
		elsif ($genus eq 'Spumavirus') {
			$tribe = 'Spuma';
			$new_genus = 'Spumavirus';
		}
		elsif ($genus eq 'Unclassified') {
			if ($class eq 'I') {
				$tribe = 'Alnilam';
			}
			elsif ($class eq 'II') {
				$tribe = 'Mintaka';
			}
			elsif ($class eq 'III') {
				$tribe = 'Alnitak';
			}
		}
		else {
			if ($class eq 'I') {
				$tribe = 'Alnilam';
			}
			elsif ($class eq 'II') {
				$tribe = 'Mintaka';
			}
			elsif ($class eq 'III') {
				$tribe = 'Alnitak';
			}
			else {
				$tribe = 'Unclassified';
			}
		}
	}
	$refseq->{virus_tribe} = $tribe;
	$refseq->{virus_genus} = $new_genus;

	unless ($state) { $state = 'exogenous' }
	$refseq->{genome_state} = $state;


	# Show the data
	unless ($name) { die; }
	unless ($state) { 
		my $question = "\n\t Set state for $name";
		#my $set_state = $console->ask_question($question);
		#$refseq->{state} = $set_state;
	}
	unless ($full_name) {
		my $question = "\n\t Full n for '$name' not set ";
		$question = "\n\t - set now as: ";
		#my $set_full_name = $console->ask_question($question);
		#$refseq->{full_name} = $set_full_name;
	}
	unless ($host_sci_name) {
		my $question = "\n\t Set host scientific name for $name: ";
		#my $host_sci_name = $console->ask_question($question);
		#$refseq->{host_sci_name} = $host_sci_name;
	}
	unless ($type) { 
		my $question = "\n\t Set genome type for $name: ";
		#my $set_type = $console->ask_question($question);
		#$refseq->{type} = $set_type;
    }


	# Deal with the genes
	my $genes_ref = $refseq->{genes};
	foreach my $gene (@$genes_ref) {
		#$devtools->print_hash($gene);
		my $gene_name = $gene->{name};
		my $full_name = $gene->{full_name};
		my $lc_gene_name = lc $gene_name;
		$gene->{name} = $lc_gene_name;
		print "\n\t ORF short $lc_gene_name full $full_name";
	}

	# Set genome type
	unless ($type) {
		$type = 'complete';
	}
	$refseq->{genome_coverage} = $type;

	if ($name =~ /HERV/ or $name eq 'RELIK' or $name eq 'SloEFV') {
		$refseq->{sequence_type} = 'consensus';

	}
	else {
		$refseq->{sequence_type} = 'locus';
	}
}


############################################################################
# DEPRECATED 
############################################################################

#***************************************************************************
# Subroutine:  make_aa_composition_table 
# Description: 
#***************************************************************************
sub make_aa_composition_table {

	my ($self, $refseq, $html_ref) = @_;

	# Write genes table
	my $indent = ' ' x 10;
	my $table_open = $writer->open_table($indent, 'data'); # Open summary table
	push (@$html_ref, $table_open);
	#my @head_row_data = qw [ property value ];
	my @aminos = qw [ A G I L P V
                      F Y W
                      D E N Q R H S T K
                      C M ];
	my @header_row = @aminos;
	unshift (@header_row, 'Gene');
	my $head_row = $writer->create_header_row(\@header_row); 
	push (@$html_ref, "\n$indent $head_row");

	# Iterate through genes
	my $orf_hash_ref = $refseq->{translated_orfs};
	my @gene_names = sort keys %$orf_hash_ref;
	foreach my $gene (@gene_names) {	
		#print "\n\t GENE $gene";	
		my $gene_sequence = $orf_hash_ref->{$gene};
		my %results;
		$seq_obj->get_aa_composition($gene_sequence, \%results);
		my @row;
		push (@row, $gene);
		foreach my $amino (@aminos) {
			my $amino_data_ref = $results{$amino};
			my $prop = $amino_data_ref->{f_percentage};
			$devtools->print_hash($amino_data_ref); die;
			#my $prop = $amino_data_ref->{f_proportion};
			unless ($prop) { $prop = '0'}
			push (@row, $prop);
		}
		my $gene_row = $writer->create_row(\@row); 
		push (@$html_ref, "\n$indent $gene_row");
	}

	# Close table
	my $table_close = $writer->close_table($indent); 
	push (@$html_ref, $table_close);
}

############################################################################
# EOF 
############################################################################
