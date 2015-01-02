#!/usr/bin/perl -w
############################################################################
# Module:      ERV_Views.pm
# Description: fxns for profiling reference sequence libraries and 
#              screening databases
# History:     November 2011: Created by Robert Gifford 
############################################################################
package ERV_Views;

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
		
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# SECTION: Top level handler functions
############################################################################

#***************************************************************************
# Subroutine:  write_retrovirus_fossil_record_page
# Description: Fossil record HTML: Retroviridae index page
#***************************************************************************
sub write_retrovirus_fossil_record_page {

	my ($self, $site_ref, $refseq) = @_;

	# Get reference sequence data
	my $name      = $refseq->{name};
	unless ($name eq 'HERV-T') { return; }

	# Get the correct database
	print "\n\t\t selecting the database for '$name'";
	my $host_views = Host_Views->new($self);	
	my $db_name = $host_views->select_screening_db_using_taxonomy($refseq);
	print "\n\t\t Got db $db_name";
	my $db = DB->new();	
	$db->load_screening_db($db_name);
	my $loci_table = $db->{loci_table};

	# Create the main panel (left) for the refseq HTML
	my @fossil_main;
	#$self->create_fossil_view($refseq, \@fossil_main);
	# Get refseq data
	my $indent = ' ' x 10;
	my $refseq_name  = $refseq->{name};
	unless ($refseq_name) { die; }
	# Write the taxonomy section
    #$self->write_retrovirus_header($refseq, $html_ref);
	#push (@$html_ref, "\n$indent <div class='separator'></div>");
	#$devtools->print_array($html_ref); die;

	# Write the main block (left panel) to the refseq HTML directory
	my $page_path = $db_pages_path . $name . '.html'; 
	$fileio->write_output_file($page_path, \@fossil_main);
	print "\n\t PATH $page_path:";
	my %counts;
	$self->get_erv_locus_counts($db, $refseq_name, \%counts);
	die;

	# Define the page
	my $full_name   = $refseq->{full_name};
	my $genus       = $refseq->{virus_genus};
	my %page;
	$page{title}  = "Paleovirology online: $name fossil record";
	my $meta_keywords = "Paleovirology, virus fossil record, $name";
	$meta_keywords .= ", $genus, fossil record";
	$page{meta_keywords} = $meta_keywords;
	my $description = "Fossil record description page for $name";
	$page{description} = $description;
	my $views_obj = $self->{views_obj};
	$views_obj->define_paleo_page($site_ref, \%page, $name, 3, 'fossil');

	# Store the path
	$refseq->{record_path} = $page_path;

}

#***************************************************************************
# Subroutine:  get_erv_locus_counts 
# Description: 
#***************************************************************************
sub get_erv_locus_counts {  

	my ($self, $db, $reference, $counts_ref) = @_;
	
	# Get objects and paths 
	my $extracted_table     = $db->{extracted_table};	
	my $loci_table          = $db->{loci_table};
	my $db_name             = $db->{db_name};
	
	print "\n\t ## Getting data $reference";
	my $field = 'COUNT(*)';
	
	# Solo LTRs	
	my @total;
	my $total_where = "WHERE assigned_to = '$reference'";
	$loci_table->select_value($field, \@total, $total_where);	
	my $total_count = shift @total;	
	unless ($total_count) { $total_count = '0'; }

	# Solo LTRs	
	my @solo;
	my $solo_where = "WHERE assigned_to = '$reference' 
					  AND   genome_structure  = 'LTR'";
	$loci_table->select_value($field, \@solo, $solo_where);	
	my $solo_count = shift @solo;	
	unless ($solo_count) { $solo_count = '0'; }
	my $solo_prop = $solo_count / $total_count;

	# Internals with LTRs	
	my @provirus;
	my $provirus_where  = "WHERE assigned_to = '$reference'"; 
	   $provirus_where .= "AND   genome_structure  != 'LTR'
					  AND   genome_structure LIKE '\%LTR\%'";
	$loci_table->select_value($field, \@provirus, $provirus_where);	
	my $provirus_count = shift @provirus;	
	unless ($provirus_count) { $provirus_count = '0'; }
	my $provirus_prop = $provirus_count / $total_count;
		
	# Internals without LTRs	
	my @internals;
	my $internals_where = "WHERE assigned_to = '$reference' 
					  AND   genome_structure NOT LIKE '\%LTR\%'";
	$loci_table->select_value($field, \@internals, $internals_where);	
	my $internals_count = shift @internals;	
	unless ($internals_count) { $internals_count = '0'; }
	my $internals_prop = $internals_count / $total_count;
	
	print "\t total ($total_count)";
	print "\t solo ($solo_count)\tprovirus ($provirus_count)\tinternal no LTR ($internals_count)";
	
	$counts_ref->{total}      = $total_count;
	$counts_ref->{solo_ltrs}  = $solo_count;
	$counts_ref->{proviruses} = $provirus_count;
	$counts_ref->{internals}  = $internals_count;
	$counts_ref->{solo_ltr_prop}  = $solo_prop;
	$counts_ref->{provirus_prop}  = $provirus_prop;
	$counts_ref->{internals_prop} = $internals_prop;

}

############################################################################
# DEVELOPMENT: PLOTS
############################################################################

#***************************************************************************
# Subroutine:  
# Description: 
#***************************************************************************
sub create_erv_plots {

	my ($self, $db_name, $references_ref, $counts_ref ) = @_;

	print "\n\n\t Writing total loci PLOT";
	my $output_path = 'site/plots/';	
	#my $output_path = $self->{output_path};
	my $total_ref  = $counts_ref->{total};
	my $data       = [ $references_ref, $total_ref ]; 
	my $graph_name = $output_path . "/$db_name" . "_total_loci.jpg";
	$self->plot_by_refseq_data($data, $db_name, $graph_name, 'marine');

	print "\n\n\t Writing solo LTR PLOT";
	my $ltr_provirus_count  = $counts_ref->{proviruses};
	$data = [ $references_ref, $ltr_provirus_count ]; 
	$graph_name = $output_path . "/$db_name" . "_proviruses.jpg";
	$self->plot_by_refseq_data($data, $db_name, $graph_name, 'lorange');

	print "\n\n\t Writing solo LTR PLOT";
	my $ltr_solo_count  = $counts_ref->{solo_ltrs};
	$data = [ $references_ref, $ltr_solo_count ]; 
	$graph_name = $output_path . "/$db_name" . "_solo_LTRs.jpg";
	$self->plot_by_refseq_data($data, $db_name, $graph_name, 'dpink');
	
	print "\n\n\t Writing internals with no LTRs PLOT";
	my $internal_count  = $counts_ref->{internals};
	$data = [ $references_ref, $internal_count ]; 
	$graph_name = $output_path . "/$db_name" . "_internals-without.jpg";
	$self->plot_by_refseq_data($data, $db_name, $graph_name, 'dgreen');

	print "\n\n\t Writing proportions plot";
	my $ltr_solo_prop = $counts_ref->{solo_ltrs_prop}; 
	my $provirus_prop = $counts_ref->{provirus_prop}; 
	my $internals_prop = $counts_ref->{internals_prop}; 
	$data = [ $references_ref, $ltr_solo_prop, $provirus_prop, $internals_prop ]; 
	$graph_name = $output_path . "/$db_name" . "_proportions.jpg";
	my @colors = qw [ dpink lorange dgreen ];
	$self->plot_by_refseq_proportions($data, $db_name, $graph_name, \@colors);

}

#***************************************************************************
# Subroutine:  
# Description: 
#***************************************************************************
sub plot_by_refseq_proportions {

	my ($self, $data, $db_name, $graph_name, $colors_ref) = @_;
	
	my $db_data = GD::Graph::Data->new($data) 
       or die GD::Graph::Data->error;
	unless ($db_name) { die; }
	my $graph = GD::Graph::bars->new(1000, 500);
	
	$graph->set(
		x_label           => 'ERV name',
		x_labels_vertical => 1,
		y_label           => 'count',
		title             => "$db_name proportions",
		transparent       => 0,
		cumulate          => 2,
		dclrs             => $colors_ref
	)
	or warn $graph->error;
    $graph->set_legend('Solo LTRs', 'Provirus', 'Internal'); 

	# Plots and print the graphs out
	$graph->plot($data) or die $graph->error; 
	my $output_path = 'site/plots';
	open(GRAPH,">$graph_name") || die "Cannot open $graph_name\n";
	print GRAPH $graph->gd->png();	
	close GRAPH;
}

#***************************************************************************
# Subroutine:  
# Description: 
#***************************************************************************
sub plot_by_refseq_data {

	my ($self, $data, $db_name, $graph_name, $color) = @_;
	
	my $db_data = GD::Graph::Data->new($data) 
       or die GD::Graph::Data->error;
	unless ($db_name) { die; }
	my $graph = GD::Graph::bars->new(800, 500);
	my @colors;
	push (@colors, $color);
	$graph->set(
		x_label           => 'ERV name',
		x_labels_vertical => 1,
		y_label           => 'count',
		title             => "$db_name data",
		transparent       => 0,
		dclrs             => \@colors 
	)
	or warn $graph->error;

	# Plots and print the graphs out
	$graph->plot($data) or die $graph->error; 
	my $output_path = 'site/plots';
	open(GRAPH,">$graph_name") || die "Cannot open $graph_name\n";
	print GRAPH $graph->gd->png();	
	close GRAPH;
}

#***************************************************************************
# Subroutine:  draw_num_by_sim_score_plots 
# Description: Create data views using charts
#***************************************************************************
sub draw_num_by_sim_score_plots {  

	my ($self, $db_name) = @_;
	
	# Get objects and paths 
	my $bar_output_path     = $self->{bar_output_path};	
	my $db = DB->new();	
	$db->load_screening_db($db_name);
	my $extracted_table     = $db->{extracted_table};	
	my $genome_chunks_table = $db->{genome_chunks_table};
	#$devtools->print_hash($db_ref);
	
	my @references;
	my @counts;
	my @fields = qw [ assigned_to ];
	$extracted_table->select_distinct(\@fields, \@references);	
	foreach my $reference (@references) {

		# set params
		my $field = 'COUNT(*)';
		my $where = "WHERE assigned_to = '$reference' 
		             AND   bit_score > 100 
					 ORDER BY assigned_to ";
		$extracted_table->select_value($field, \@counts, $where);	
	}
		
	# Set up the bar graph
	my $graph = GD::Graph::bars->new(800, 500)
		or die "\nCan't create graph\n\n";
	$graph->set( title          => "$db_name summary",
				x_label         => "reference",
				x_labels_vertical  => 1,
				y_label         => "number of hits",
				show_values     => 1,
				);
		
	# plot the graph
	my $gd = $graph->plot([ \@references, \@counts ])
		or die "\nCan't plot graph\n\n";

	# print the graph out
	my $graph_name = $db_name . "_hit_summary.jpg";
	my $output_path = $bar_output_path . $graph_name;
	open(GRAPH,">$output_path") || die "Cannot open $output_path\n";
	print GRAPH $graph->gd->jpeg(100);	

	my $key = $db_name;
	$self->{$key} = $output_path;
}

############################################################################
# EOF 
############################################################################
