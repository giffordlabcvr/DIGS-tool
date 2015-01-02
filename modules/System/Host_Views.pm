#!/usr/bin/perl -w
############################################################################
# Module:      Host_Views.pm
# Description: fxns for profiling reference sequence libraries and 
#              screening databases
# History:     November 2011: Created by Robert Gifford 
############################################################################
package Host_Views;

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
# SECTION: GETTING HOST TAXONOMY FROM THE NCBI TAXONOMY DB
############################################################################

#***************************************************************************
# Subroutine:  setup_host_taxonomy 
# Description: 
#***************************************************************************
sub setup_host_taxonomy {

	my ($self, $retroviridae_ref) = @_;

	# Instantiate DB obj
	my $db_obj  = DB->new($self);

	# Create data structure 
	$db_obj->load_ncbi_taxonomy();

	# Get tables
	my $names_table = $db_obj->{names_table}; 
	my $nodes_table = $db_obj->{nodes_table}; 

	# Iterate through the viruses
	my @viruses = sort keys %$retroviridae_ref;
	foreach my $virus (@viruses) {

		print "\n\t # Getting taxonomy for '$virus'";
		my $refseq = $retroviridae_ref->{$virus};
		my $refseq_name   = $refseq->{name};
		my $host_species  = $refseq->{host_sci_name};
		$host_species =~ s/_/ /g;
		my @fields = qw [ taxonid name uniquename class ];
		my $where = "WHERE name = '$host_species'";
		my %data; 
		$names_table->select_row(\@fields, \%data, $where);
		my %taxonomy;
		my $taxonid = $data{taxonid};
		#$devtools->print_hash(\%data); die;
		unless ($taxonid) { die "\n\t No taxon_id found for species $host_species" ; }
		my $limit = 'class';
		my $i = 0;
		$self->retrieve_taxonomy($db_obj, $taxonid, \%taxonomy, $limit, $i);
		#$devtools->print_hash(\%taxonomy);
		# die;
		$refseq->{host_superclass}  = $taxonomy{superclass};
		$refseq->{host_class}       = $taxonomy{class};
		my $order = $taxonomy{order};
		$refseq->{host_order}       = $order;
		$refseq->{host_suborder}    = $taxonomy{suborder};
		$refseq->{host_family}      = $taxonomy{family};
		$refseq->{host_superorder}  = $taxonomy{superorder};
		$refseq->{host_genus}       = $taxonomy{genus};
		$refseq->{taxonomy_string}  = $taxonomy{taxonomy_string};
	}
}

#***************************************************************************
# Subroutine:  retrieve_taxonomy 
# Description:  
#***************************************************************************
sub retrieve_taxonomy {

	my ($self, $db_obj, $taxonid, $taxonomy, $limit_rank, $i) = @_;

	$i++;
	if ($i > 20) { return; } # Prevent deep recursion

	#print "\n\t Getting data for $taxonid";	
	my $names_table = $db_obj->{names_table}; 
	my $nodes_table = $db_obj->{nodes_table}; 
	#sleep 1;

	my %node_data;
	my @node_fields = qw [ parenttaxonid rank ]; 
	my $where  = "WHERE taxonid = '$taxonid' ";
	$nodes_table->select_row(\@node_fields, \%node_data, $where);
	my $rank            = $node_data{rank};
	my $parenttaxonid   = $node_data{parenttaxonid};
	unless ($parenttaxonid and $rank) {
		print "\n\t got nothing for id $taxonid";
		$devtools->print_hash($taxonomy); die;
	}

	my %data;
	my @fields = qw [ taxonid name uniquename class ];
	$where = "WHERE taxonid = '$taxonid' AND class = 'scientific name' ";
	$names_table->select_row(\@fields, \%data, $where);
	my $name = $data{name};


	unless ($taxonid and $name) {
		$devtools->print_hash($taxonomy); die;
		return;
	}
	else {

		unless ($rank eq 'no rank') {
			$taxonomy->{$rank} = $name;
		}

		my $taxa_string = $taxonomy->{taxonomy_string};
		unless ($taxa_string) {
			$taxa_string = $name;
		}
		else {
			$name .= ":$taxa_string";
		}
		$taxonomy->{taxonomy_string} = $name;

		#if ($rank eq $limit_rank) {
		if ($rank eq 'domain') {
			return;
		}
		elsif ($parenttaxonid) {
			#print "\n id $taxonid = $name ($rank), parent = $parenttaxonid)";
			$self->retrieve_taxonomy($db_obj, $parenttaxonid, $taxonomy, $limit_rank, $i);
		}
	}
}

############################################################################
# SECTION: Make host DB HTML
############################################################################

#***************************************************************************
# Subroutine:  make_host_db_html
# Description: 
#***************************************************************************
sub make_host_db_html {

	my ($self, $site_ref) = @_;

	# Create the HTML
	my @main;
	my $indent = ' ' x 10;
	push (@main, "\n$indent <h3>Genome database summary</h3>");
	push (@main, "\n$indent <div class='separator'></div>");
	
	# Get the hosts organsims sorted into groups
	my $db_obj = DB->new();
	$db_obj->load_genome_database();
	my $table = $db_obj->{genome_chunks_table};
	my @organisms;	
	my %sorted;
	$table->select_distinct_single_field('organism', \@organisms);	
	$self->sort_hosts_by_taxonomy(\@organisms, \%sorted);

	push (@main, "\n$indent <br><br>");
	push (@main, "\n$indent <h3>Class Mammalia</h3>");
	push (@main, "\n$indent <div class='separator'></div>");
	my $superorder_mammals = $sorted{by_superorder_mammals};
	#$devtools->print_hash($superorder_mammals); die;
	my @table1;
	$self->create_host_genome_table($table, $superorder_mammals, \@table1, 'Superorder');
	push (@main, @table1); 
	
	push (@main, "\n$indent <h3>Class Aves</h3>");
	push (@main, "\n$indent <div class='separator'></div>");
	my $superorder_birds = $sorted{by_superorder_birds};
	#$devtools->print_hash($super_order_birds); die;
	my @table2;
	$self->create_host_genome_table($table, $superorder_birds, \@table2, 'Superorder');
	push (@main, @table2); 

	push (@main, "\n$indent <div class='separator'></div>");

	# Write the main block (left panel) to the refseq HTML directory
	my $page_name = 'genomes';
	my $page_path = $site_path . "main/screening/$page_name.html"; 
	$fileio->write_output_file($page_path, \@main);
	
	# Define the page
	my %page;
	$page{title} = 'Paleovirology online: Retrovirus fossil record';
	$self->define_paleo_page($site_ref, \%page, $page_name, 2, 'fossil');

	print "\n\n\t Recreated genome DB HTML! under $page_path\n\n\n";
}

#***************************************************************************
# Subroutine:  create_host_genome_table
# Description:  
#***************************************************************************
sub create_host_genome_table {

	my ($self, $genome_chunks_table, $hosts, $table_ref, $rank) = @_;

	# Create genomes table
	my $indent = ' ' x 10;
	my @names = sort keys %$hosts;
	foreach my $name (@names) {
	
		push (@$table_ref, "\n$indent <h4>$rank $name</h4>");
		push (@$table_ref, "\n$indent <div class='divider'></div>");
		my $table_open = $writer->open_table($indent, 'text'); 
		push (@$table_ref, $table_open);

		# Open summary table
		my $column_tags =  "<col width=\"40\" />
							<col width=\"20\" />
							<col width=\"20\" />
							<col width=\"20\" />";
		push (@$table_ref, $column_tags);

		my @head_row_data = qw [ Organism ];
		push (@head_row_data , 'type');
		push (@head_row_data , '# scaffolds');
		push (@head_row_data , '# bases');
		my $head_row = $writer->create_header_row(\@head_row_data); 
		push (@$table_ref, $head_row);
		
		my $species_ref = $hosts->{$name};
		my @alphabetical = sort @$species_ref;
		foreach my $organism (@alphabetical) {    

			# Connect to the database and select rows from the 'chunks' table
			my @fields = qw [ organism version ];
			my @genome_versions;
			my $where = " WHERE organism = '$organism'";
			$genome_chunks_table->select_distinct(\@fields, \@genome_versions, $where);	
			foreach my $version_ref (@genome_versions) {	
			
				my $version  = $version_ref->{version};
				print "\n\t ## $name organism '$organism' ($version)";
				my @data;
				@fields = qw [ number_scaffolds total_bases source_type ];
				my $where  = " WHERE Organism = '$organism' AND Version = '$version' ";
				$genome_chunks_table->select_rows(\@fields, \@data, $where);	
				my $num = scalar @data;
				#print "\n\t number $num";

				$organism =~ s/_/ /g;
				my $total_bases = 0;
				my $total_scaffolds = 0;
				my $type;
				foreach my $chunk_ref (@data) {
					my $scaffolds   = $chunk_ref->{number_scaffolds};
					my $bases       = $chunk_ref->{total_bases};
					$type           = $chunk_ref->{source_type};
					$type =~ s/_/ /g;
					$total_bases    = $total_bases + $bases;	
					$total_scaffolds = $total_scaffolds + $scaffolds;
				}
				
				my @new_row = ( "<i>$organism</i>", $type, $total_scaffolds, $total_bases );
				my $row = $writer->create_row(\@new_row); 
				push (@$table_ref, $row);
			}
		}
		# Close table
		my $table_close = $writer->close_table($indent); 
		push (@$table_ref, $table_close);
		push (@$table_ref, "\n$indent <br>\n\n");
	}
}

############################################################################
# SORTing by HOST taxonomy
############################################################################

#***************************************************************************
# Subroutine:  sort_hosts_by_taxonomy 
# Description:  
#***************************************************************************
sub sort_hosts_by_taxonomy {

	my ($self, $host_list_ref, $sorted_ref) = @_;
	
	# Instantiate DB obj
	my $db_obj  = DB->new($self);

	# Create data structure 
	$db_obj->load_ncbi_taxonomy();

	# Get tables
	my $names_table = $db_obj->{names_table}; 
	my $nodes_table = $db_obj->{nodes_table}; 

	# Iterate through the species
	my %taxonomy;
	foreach my $organism (@$host_list_ref) {	
		my $latin_binomial = $organism;
		$latin_binomial =~ s/_/ /g;	
		my $where = "WHERE name = '$latin_binomial'";
		my %data; 
		my @fields = qw [ taxonid name uniquename class ];
		$names_table->select_row(\@fields, \%data, $where);
		my %taxonomy_info;
		my $taxonid = $data{taxonid};
		unless ($taxonid) { die "\n\t No taxon_id found for species '$latin_binomial'" ; }
		my $limit = 'class';
		my $i = 0;
		$self->retrieve_taxonomy($db_obj, $taxonid, \%taxonomy_info, $limit, $i);
		#$devtools->print_hash(\%taxonomy_info);
		$taxonomy{$organism} = \%taxonomy_info; 
	}

	# Sort into classes
	my %by_class;
	$self->sort_by_rank(\@$host_list_ref, \%taxonomy, 'class', \%by_class);
	my $birds_ref   = $by_class{Aves};
	my $mammals_ref = $by_class{Mammalia};
	print "\n\t BY CLASS";
	#$devtools->print_hash(\%by_class); 

	# Sort into superorders 
	my %by_superorder_mammals;
	$self->sort_by_rank($mammals_ref, \%taxonomy, 'superorder', \%by_superorder_mammals);
	#$devtools->print_hash(\%by_superorder_mammals); exit;
	#$by_superorder_mammals{Marsupials} = $by_superorder_mammals{undefined};
	#delete $by_superorder_mammals{undefined};
	my %by_superorder_birds;
	$self->sort_by_rank($birds_ref, \%taxonomy, 'superorder', \%by_superorder_birds);
	#$devtools->print_hash(\%by_superorder_mammals); die;

	# Sort into orders 
	my %by_order_mammals;
	$self->sort_by_rank($mammals_ref, \%taxonomy, 'order', \%by_order_mammals);
	my %by_order_birds;
	$self->sort_by_rank($birds_ref, \%taxonomy, 'order', \%by_order_birds);

	# Sort into suborders 
	
	# Store data
	$sorted_ref->{by_class} = \%by_class;
	$sorted_ref->{by_superorder_mammals} = \%by_superorder_mammals;
	$sorted_ref->{by_superorder_birds} = \%by_superorder_birds;
	$sorted_ref->{by_order_mammals} = \%by_order_mammals;
	$sorted_ref->{by_order_birds} = \%by_order_birds;

}

#***************************************************************************
# Subroutine:  sort_by_rank
# Description: 
#***************************************************************************
sub sort_by_rank {  

	my ($self, $species_list, $taxonomy_ref, $rank, $sorted) = @_;
	
	# Iterate through the species
	foreach my $organism (@$species_list) {	
	
		my $species_taxonomy = $taxonomy_ref->{$organism};
		unless ($species_taxonomy) { die; }
		my $rank_name = $species_taxonomy->{$rank};
		print "\n\t GOT $rank $rank_name for $organism";
		unless ($rank_name) { 
			$rank_name = 'undefined';
		}
		if ($sorted->{$rank_name}) {
			my $array_ref = $sorted->{$rank_name};
            push (@$array_ref, $organism);
		}	
		else { 
			my @array;
            push (@array, $organism);
			$sorted->{$rank_name} = \@array;
		}
	}
	#$devtools->print_hash($sorted); 
}

############################################################################
# VIEW OPTIONS 
############################################################################

#***************************************************************************
# Subroutine:  make_screendb_html
# Description: handler for tools to show views on screening DBs
#***************************************************************************
sub make_screen_db_html {

	my ($self) = @_;

	# Connect to the screening DB database (keeps track of screening DBs)
	my $server   = 'localhost';
    my $username = 'root';
    my $password = 'palb0x';
	my $dbh = DBI->connect("dbi:mysql:ScreeningDB:$server", $username, $password);
	unless ($dbh) {	die "\n\t # Couldn't connect to ScreeningDB database\n\n"; }
	
	# Instantiations
	my $seq_obj    = Sequence->new();

	# write the header
	my @html;
	
	# Execute the query
    my $query = "SELECT DB_name, Purpose from DB_detail"; 
	my $sth = $dbh->prepare($query);
    unless ($sth->execute()) { print $query; exit;}	
    my %data;

	# Write genes table
	my $width = 550;
	my $table_open = $writer->open_table($width); # Open summary table
	push (@html, $table_open);
	my @head_row_data = qw [ Database Description ];
	my $head_row = $writer->create_header_row(\@head_row_data); 
	push (@html, $head_row);
	
	while (my $row = $sth->fetchrow_arrayref()) {
        my $db_name = @$row[0];
        my $purpose = @$row[1];
		my $db_page = $db_name . '_screendb.html';
		#$self->write_screening_db_html($db_name, $db_page);
		my $f_db_name = "<a href='$db_page'>$db_name</a>"; 
		my @table_row;
		push (@table_row, $f_db_name);
		push (@table_row, $purpose);
		print "\n\t$db_name: $purpose\n";
		#print "\n\t$string";
		#push (@access, $string);
		#$devtools->print_hash($gene_ref); die;
		my $db_row = $writer->create_row(\@table_row); 
		push (@html, $db_row); 
    }
	#my $path  = $site_path . 'screening_dbs.html';
	#print "\n\t $path \n\n";
	#$fileio->write_output_file($path, \@html);
	
	# Load the screening DB
	# Show summary table
	# Show retreive fe	
	# Select DISTINCT assigned to
	# Select DISTINCT assigned to groups
	# Select DISTINCT assigned to subgroups
}

#***************************************************************************
# Subroutine:  select_screening_db_using_taxonomy
# Description: $refseq: a reference sequence object (WITH HOST TAXONOMY DATA)
#***************************************************************************
sub select_screening_db_using_taxonomy {

	my ($self, $refseq) = @_;

	my $pause;
	# Get refseq data
	my $name       = $refseq->{name};	
	my $class      = $refseq->{host_class};	
	my $superorder = $refseq->{host_superorder};	
	my $order      = $refseq->{host_order};	
	my $suborder   = $refseq->{host_suborder};	
	my $family     = $refseq->{host_family};	
	my $species    = $refseq->{host_sci_name};	
	$species =~ s/_/ /g;

	# Class to DB 
	my %class_to_db;
	$class_to_db{Aves}           = 'Aves';
	$class_to_db{Actinopterygii} = 'Fish';
	$class_to_db{Amphibia}       = 'Amphibia';
	
	# Superorder to DB 
	my %superorder_to_db;
	$superorder_to_db{Afrotheria}       = 'Afrotheria';
	$superorder_to_db{Xenarthra}        = 'Xenarthra';
	$superorder_to_db{Lepidosauria}     = 'Reptiles';
	$superorder_to_db{Acanthopterygii}  = 'Fish';

	# Order to DB 
	my %order_to_db;
	$order_to_db{Insectivora}    = 'Insectivora';
	$order_to_db{Carnivora}      = 'Carnivora';
	$order_to_db{Chiroptera}     = 'Chiroptera';
	$order_to_db{Lagomorpha}     = 'Lagomorpha';
	$order_to_db{Rodentia}       = 'Rodents';
	$order_to_db{Scandentia}     = 'Scandentia';
	$order_to_db{Cetacea}        = 'Cetartiodactyla';
	$order_to_db{Diprotodontia}  = 'Australidelphia';
	$order_to_db{Dasyuromorphia} = 'Australidelphia';
	$order_to_db{Didelphimorphia}= 'Ameridelphia';
	$order_to_db{Perissodactyla} = 'Perissodactyla';
	$order_to_db{Monotremata}    = 'Monotremata';
	 
	# Suborder to DB 
	my %suborder_to_db;
	$suborder_to_db{Platyrrhini}    = 'Platyrrhini';
	$suborder_to_db{Strepsirrhini}  = 'Lemuriformes';
	$suborder_to_db{Ruminantia}     = 'Cetartiodactyla';

	# Family to DB 
	my %family_to_db;
	$family_to_db{Suidae}       	= 'Cetartiodactyla';
	$family_to_db{Hominidae}        = 'Catarrhini';
	$family_to_db{Cercopithecidae}  = 'Catarrhini';
	$family_to_db{Tarsiidae}        = 'Tarsidae';
	$family_to_db{Pitheciidae}      = 'Platyrrhini';
	$family_to_db{Atelidae}         = 'Platyrrhini';
	$family_to_db{Hylobatidae}      = 'Catarrhini';
	$family_to_db{Cebidae}          = 'Platyrrhini';
	
	my $db_name;
	if ($class) {
		print "\n\t class    $class";
		if ($class_to_db{$class}) {
			$db_name = $class_to_db{$class};
		}
	}
	if ($superorder) {
		if ($superorder_to_db{$superorder}) {
			$db_name = $superorder_to_db{$superorder};
		}
	}
	if ($order) {
		print "\n\t order    $order";
		if ($order_to_db{$order}) {
			$db_name = $order_to_db{$order};
		}
	}
	if ($suborder) {
		print "\n\t suborder $suborder";
		if ($suborder_to_db{$suborder}) {
			$db_name = $suborder_to_db{$suborder};
		}
	}
	if ($family) {
		print "\n\t family $family";
		if ($family_to_db{$family}) {
			$db_name = $family_to_db{$family};
		}
	}

	unless ($db_name) {
		print "\n\n\n\\t No database for '$species'\n\n";
	}

	return $db_name;
}

############################################################################
# EOF 
############################################################################
