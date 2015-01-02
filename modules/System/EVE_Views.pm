#!/usr/bin/perl -w
############################################################################
# Module:      EVE_Views.pm
# Description: fxns for profiling reference sequence libraries and 
#              screening databases
# History:     November 2011: Created by Robert Gifford 
############################################################################
package EVE_Views;

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
# SECTION: EVE pages
############################################################################

#***************************************************************************
# Subroutine:  make_eve_fossil_record_html 
# Description: make the EVE fossil record HTML 
#***************************************************************************
sub make_eve_fossil_record_html {

	my ($self, $site_ref) = @_;
	return;

	# Load the EVE reference sequences
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
	my @ranks = qw [ virus_family virus_genus virus_subgroup name ];
	$refseqlib_obj->order_by_taxonomy(\@refseqs, \@ranks, \%by_taxonomy, \%flat);

	# Set up for colecting statistics on the family
	my %eve_data;
	$self->setup_rv_fossil_stats(\%eve_data);

	# Get host taxonomy information
	print "\n\t Getting host taxonomic information";
	$self->setup_host_taxonomy(\%flat); 

	# Iterate through the viruses, writing the page for each 
	my @viruses = sort keys %flat;
	foreach my $virus (@viruses) {
		
		my $refseq = $flat{$virus};
		unless ($refseq)  { die; }
		print "\n\t Processing '$virus'";
		
		# Get host lineage data  
		print "\n\t\t Getting '$virus' host lineage data";
		$self->stats_by_host_group($refseq, \%eve_data);
		
		# Get statistics for this EVE reference sequence
		print "\n\t\t Getting '$virus' statistics";
		$self->stats_by_virus_group($refseq, \%eve_data);
		
		print "\n\t\t Writing Refseq reference page for '$virus'";
		$self->write_eve_refseq_page($site_ref, $refseq);
	}
	
	# Write RV fossil index page
	$self->make_eve_fossil_index($site_ref, \%flat, \%by_taxonomy);

	# Write the overview page
	$self->make_eve_fossil_overview($site_ref, \%eve_data, \%flat);
	
	# Create the html
	$self->make_paleo_website_html($site_ref);

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

############################################################################
# EOF 
############################################################################
