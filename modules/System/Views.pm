#!/usr/bin/perl -w
############################################################################
# Module:      Views.pm
# Description: fxns for profiling reference sequence libraries and 
#              screening databases
# History:     November 2011: Created by Robert Gifford 
############################################################################
package Views;

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
use Program::System::RV_Views;   # Retrovirus reference sequences views
use Program::System::ERV_Views;  # Fossil record data views
use Program::System::EVE_Views;  # Non-retroviral EVE views
use Program::System::Host_Views; # Host species data + views

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
# Subroutine:  run_view_fxn 
# Description: top level handler 
#***************************************************************************
sub run_view_function {

	my ($self, $glue_mode) = @_;
	
	# Set the usage statement
   	my $USAGE .= "\n\t  v = Views.pm website management\n";
   	$USAGE    .= "\n\t     1 = make static site HTML";
   	$USAGE    .= "\n\t     2 = make RV fossil record HTML"; 
   	$USAGE    .= "\n\t     3 = genome screening HTML";
   	$USAGE    .= "\n\t     4 = refresh entire site"; 
   	$USAGE    .= "\n\t     5 = recreate BLAST libraries";
	$USAGE    .= "\n\n";

	# Load data structure describing basic website
	# TO DO: Move from here, into fxns, dpoesnt always need 'load site'
	$self->show_title();
	my %site;
	my @pages;
	my %pages;
	$site{pages_hash}  = \%pages;
	$site{pages_array} = \@pages;

	# Hand of to function
	if ($glue_mode eq 1) { 
		# Write the static pages only 
		$self->load_paleo_site(\%site);
		$self->make_paleo_website_html(\%site);
	}
	# Single process calls
	elsif ($glue_mode eq 2) { # Refresh screening DB pages
		print "\n\t # Making reference genome HTML";
		my $rv_views = RV_Views->new($self);
		$rv_views->{views_obj} = $self;
		$rv_views->make_retrovirus_record_html(\%site); 
		#$self->make_eve_fossil_record_html(\%site); 
	}
	elsif ($glue_mode eq 3) { 
		$self->make_screening_html(\%site);
	}
	elsif ($glue_mode eq 4) { 
		$self->make_entire_web_site(\%site);
	}
	elsif ($glue_mode eq 5) { 
    	# Create the BLAST libraries 
		$self->create_blast_libraries();
	}
	else { die $USAGE; }	
}

############################################################################
# SECTION: Second level handler functions
############################################################################

#***************************************************************************
# Subroutine:  make_paleo_website_html 
# Description: make the (non fossil record) pages of the paleo website
#***************************************************************************
sub make_paleo_website_html {

	my ($self, $site_ref) = @_;

	# Get site pages
	my $pages_array_ref = $site_ref->{pages_array};
	unless ($pages_array_ref) { 
		# Something went wrong if no pages were loaded
		die "\n\t No pages were loaded in the 'site' data structure\n\n";
	}

	# Write the pages 
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
# Subroutine:  make_screening_html 
# Description: refresh genome screening page(s)  - # views on Screening DBs
#***************************************************************************
sub make_screening_html {

	my ($self, $site_ref) = @_;
	
	# Create html describing host DB
	my $macroevolution;
	print "\n\t # Making screening DB HTML";
	$self->make_host_db_html($site_ref);
	print "\n\t # Making server genomes HTML";
	$self->make_screendb_html($site_ref);
}

#***************************************************************************
# Subroutine:  make_entire_web_site 
# Description: refresh paleovirology website
#***************************************************************************
sub make_entire_web_site {

	my ($self, $site_ref) = @_;
	
	print "\n\t # Making web site HTML\n";
	$self->make_paleo_website_html($site_ref);
	
	print "\n\t # Making reference genome HTML";
	my $rv_views = RV_Views->new($self);
	$rv_views->{views_obj} = $self;
	$rv_views->make_retrovirus_record_html($site_ref); 
	
	# Make server HTML
	print "\n\t # Making server genomes HTML";
	$self->make_host_db_html($site_ref);
	
	# Make screening DB HTML
	print "\n\t # Making screening DB HTML";
	$self->make_screen_db_html($site_ref);
}

############################################################################
# SECTION: Load the static pages and research/news 'streams' from the site
############################################################################

#***************************************************************************
# Subroutine:  load_paleo_site
# Description: load paleovirology site
#***************************************************************************
sub load_paleo_site {

	my ($self, $site_ref) = @_;

	# Store the site pages for PALEO
	my @paleo    = qw [ about background teaching ];
	my @lab      = qw [ gifford niewiadomska ];
	my @fossil   = qw [ viral_fossil_record non_retroviral_eves
                        definitions RV_refseq_definitions 
                        ERV_download  ];
	my @tool     = qw [ tools_overview fasta_tool blast_rvs 
                        glue_tool glue_msa_tool refseq_define 
                        glue_tool_release_notes fasta_tool_release_notes
                        blast_tool_release_notes refseq_tool_release_notes
					 ];
	my @screening = qw [ screening genomes ];
	my @external  = qw [ external 404 ];

	# Define the main, stable pages
	$self->define_paleo_pages(\@paleo, $site_ref,  2, 'paleo');
	$self->define_paleo_pages(\@lab, $site_ref, 2, 'lab');
	$self->define_paleo_pages(\@fossil, $site_ref, 2, 'fossil');
	$self->define_paleo_pages(\@tool, $site_ref, 2, 'tool');
	$self->define_paleo_pages(\@screening, $site_ref, 2, 'screening');
	$self->define_paleo_pages(\@external, $site_ref, 2, 'external');

	# Define the streams (i.e. posts in sequence)
	my @post_files;
	my @posts;
	my $research_post_dir_path = $site_path . 'main/posts/';
	$fileio->read_directory_to_array($research_post_dir_path, \@post_files);
	foreach my $file (@post_files) {
		my @file_bits = split (/\./, $file);
		my $file_type = pop @file_bits;
		my $file_stem = join('.', @file_bits);
		my %page;
		my $file_path = $research_post_dir_path . $file;
		my @file;
		$fileio->read_input_file($file_path, \@file);
		foreach my $line (@file) {
            # Look for title set token
            if (($line =~ /<!--/) and ($line =~ /TITLE/)) {
				my @line = split (/'/, $line);
				my $title = $line[1];
				#print "\n\t Setting post $file_stem title '$title'";
				$page{title} = "Paleovirology online: $title";
			}

            # Look for keywords set token
            if (($line =~ /<!--/) and ($line =~ /TAGS/)) {
				my @line = split (/'/, $line);
				my $keywords = $line[1];
				#print "\n\t Setting post $file_stem keywords '$keywords'";
				$page{meta_keywords} = $keywords;
			}
		}
		$self->define_paleo_page($site_ref, \%page, $file_stem, 3, 'research');
		push (@posts, $file_stem);
	} 

	# Create the update summary stream
	my %stream;
	$stream{page_names}  = \@posts;
	$stream{main_path}   = $site_path . 'main/posts/'; 
	$stream{outpath}     = $site_path . 'main/paleo/'; 
	$stream{stream_path} = $site_path . 'html'; 
	$stream{make_index}  = 'true';
	$stream{post_limit}  = 4;
	$stream{stream_name} = 'paleovirology';
	$self->create_stream($site_ref, \%stream, 'paleo');
	
	# Create the research summary stream
	my %stream2;
	$stream2{page_names}  = \@posts;
	$stream2{main_path}   = $site_path . 'main/posts/'; 
	$stream2{outpath}     = $site_path . 'main/lab/'; 
	$stream2{stream_path} = $site_path . 'html'; 
	$stream2{post_limit}  = 3;
	$stream2{stream_name} = 'research';
	$stream2{include_token} = 'RESEARCH POST';
	#$self->create_stream($site_ref, \%stream2, 'research', 2);
}

#***************************************************************************
# Subroutine:  create_stream
# Description: create a stream of trucated posts linked to more complete doc
#***************************************************************************
sub create_stream {

	my ($self, $site_ref, $stream_ref, $style) = @_;

	my $indent = ' ' x 10;
	my $stream_pages_ref = $stream_ref->{page_names};
	my $main_path     = $stream_ref->{main_path};
	my $stream_path   = $stream_ref->{stream_path};
	my $outpath       = $stream_ref->{outpath};
	my $post_limit    = $stream_ref->{post_limit};
	my $stream_name   = $stream_ref->{stream_name};
	my $make_index    = $stream_ref->{make_index};
	my $include_token = $stream_ref->{include_token};
	unless ($main_path and $post_limit and $stream_name and $outpath) { die; }

	my %stream;
	my @current_stream;
	my $i = 0; # Post count
	my $j = 1; # Stream count
	my $total_pages = scalar @$stream_pages_ref;
	my @sorted = reverse sort @$stream_pages_ref;
	my %stream_to_main_path;
	my %stream_to_final_path;
	my %stream_to_content;
	foreach my $stream_page (@sorted) {
		
		$i++;	
		#print "\n\t Adding post '$stream_page' to stream $j";
		
		# Read the page
		my $post_path = $main_path . $stream_page . '.html';
		my @post;
		$fileio->read_input_file($post_path, \@post);		
		
		# Put the truncation in (if there is one)
		my @truncated;
		my $truncated      = undef;
		my $found_token    = undef;
		my $mainpost_flag  = undef;
		my $added_stream   = undef;
		foreach my $line (@post) {

			if (($line =~ /<!--/) and ($line =~ /TOP-POST/)) {
				$mainpost_flag = 1;
			}

            # Look for include/exclude token
			if ($include_token) {
				if (($line =~ /<!--/) and ($line =~ /$include_token/)) {
					$found_token = 'true';
				}
			}

            # Look for truncate token
			if (($line =~ /<!--/) and ($line =~ /TRUNCATE/)) {
				$truncated = 'true';
				last; 
			}
			else { push (@truncated, $line); }
		}
					
		# Add the truncated post
		if ($include_token) {
			if ($found_token) {
				$added_stream = 'true';
				push(@current_stream, @truncated);
			}
		}
		else {
			$added_stream = 'true';
			push(@current_stream, @truncated);
		}
		# Add the link to the long post
		if ($truncated and $added_stream) {
			my $absolute = $url . "site/html/posts/";
			my $link = '           ';
			$link   .= '<br><br><a href="' .  $absolute . $stream_page . '.html">';
			$link   .= '[Read more...]';
			$link   .= " </a><br>\n";
			push (@current_stream, $link);
		}
	
		if ($added_stream) {
			# Add the divider
        	my $separator = "         <br>\n";
			push (@current_stream, $separator);
			push (@current_stream, $separator);
			
			# See if we should start a new stream
			if ($i % $post_limit == 0) {
				$stream{$j} = \@current_stream;
				my $main_path;
				if ($j eq 1 and $make_index) {
					$stream_to_final_path{$j} = $url . 'index.html';
					$main_path = $outpath . "index.html";
				}
				else {
					my $name = $stream_name . '_' . $j;
					$stream_to_final_path{$j} = $url . "$stream_path/$name.html";
					$main_path = $outpath . "$name.html";
				}
				#print "\n\t\t WRITING stream $j to '$main_path";
				$stream_to_main_path{$j} = $main_path; # Track streams to paths
				my @stream_fin = @current_stream;
				$stream_to_content{$j} = \@stream_fin; # Track streams to content
        	    @current_stream = ();
				$j++;
			}
		}
		else { $i = $i - 1; }
	}

	# get any stragglers
	my $num_current = scalar @current_stream;
	if ($num_current > 1) {
		$stream{$j} = \@current_stream;
		my $main_path;
		if ($j eq 1 and $make_index) {
			$stream_to_final_path{$j} = $url . 'index.html';
			$main_path = $outpath . "index.html";
		}
		else {
			my $name = $stream_name . '_' . $j;
			$stream_to_final_path{$j} = "$url" . "$stream_path/$name.html";
			$main_path = $outpath . "$name.html";
		}
		#print "\n\t\t WRITING stream $j to '$stream_path";
		my @stream_fin = @current_stream;
		$stream_to_main_path{$j} = $main_path; # Track streams to paths
		$stream_to_content{$j}    = \@stream_fin; # Track streams to content
        @current_stream = ();
	}
	#$devtools->print_hash(\%stream_to_content); die;

	# Write the streams with the end linking section
	my @streams = sort by_number keys %stream_to_main_path;
	foreach my $stream (@streams) {
		my $stream_path = $stream_to_main_path{$stream};
		my $content_ref = $stream_to_content{$stream};
		my $link=$self->create_page_footer_link(\@streams, \%stream_to_final_path, $stream);	
		push (@$content_ref, "$link<br>");
		push (@$content_ref, "$indent<div class=\"separator\"></div>");
		$fileio->write_output_file($stream_path, $content_ref);
	}

	# Create stream pages
	if ($j < 1) { die; }
	my $stream_count = 0;
	foreach my $stream (@streams) {
		#print "\n\t Setting post $stream";
		my %page;
		my $stream_path = $stream_to_main_path{$stream};
		$stream_count++;
		if ($stream_count eq 1 and $make_index) {
			#print "\n\t Creating index stream";
			$self->define_paleo_page($site_ref, \%page, 'index', 1, $style);
		}
		else {
			my $name = $stream_name . '_' . $stream_count;
			#print "\n\t Creating stream $name";
			$self->define_paleo_page($site_ref, \%page, $name, 2, $style);
		}
	} 
}

#***************************************************************************
# Subroutine:  create_page_footer_link
# Description: 
#***************************************************************************
sub create_page_footer_link {

	my ($self, $array_ref, $hash_ref, $stream) = @_;

	# Create the page foot stream links
	my $indent = ' ' x 10;
	my $link = "\n$indent  <br><div class=\"centered\">";
	$link   .= "\n$indent  Page:&nbsp;";
	foreach my $stream_num (@$array_ref) {
		my $stream_path = $hash_ref->{$stream_num};
		if ($stream_num eq $stream) {
			$link .= "\n$indent  <u>$stream_num</u>&nbsp;";
		}
		else {
			$link .= "\n$indent  <a href=\"$stream_path\">$stream_num</a>&nbsp;";
		}
	}	
	$link   .= "\n$indent</div>\n\n";
	return $link;
}

############################################################################
# SECTION: Define a web page in paleovirology site format 
############################################################################

#***************************************************************************
# Subroutine:  define paleo pages
# Description: just a handler loop for the 'define_paleo_page' sub 
#***************************************************************************
sub define_paleo_pages {

	my ($self, $pages_ref, $site_ref, $level, $style) = @_;
	
	# Iterate through the pages
	#print "\n\t # Defining $style pages (level $level)";
	foreach my $page (@$pages_ref) {
		# Define the page
		my %page;
		$self->define_paleo_page($site_ref, \%page,  $page, $level, $style);
	}
}

#***************************************************************************
# Subroutine:  define paleo page
# Description: 
#***************************************************************************
sub define_paleo_page {

	my ($self, $site_ref, $page_ref, $page, $level, $style) = @_;
	
	my $ssi_path = $site_path . 'ssi/';
	
	#print "\n\t ####### DOING PAGE $page";
	#	$devtools->print_hash($titles_ref);
	#	$devtools->print_hash($page_ref);
	# Set the universal values	
	$page_ref->{name}          = $page;
	$page_ref->{doctype}       = $ssi_path . 'doctype.html';
	$page_ref->{meta_robots}   = $ssi_path . 'meta-robots.html';
	$page_ref->{meta_content}  = $ssi_path . 'meta-content.html';
	$page_ref->{meta_author}   = $ssi_path . 'meta-author.html';
	$page_ref->{meta_owner}    = $ssi_path . 'meta-owner.html';
	$page_ref->{urchin}        = $ssi_path . 'urchin.html';
	$page_ref->{footer}        = $ssi_path . 'footer.html';
	$page_ref->{top_nav}       = $ssi_path . 'top_nav.html';
	$page_ref->{right}         = $ssi_path . 'sidebar.html';
	$page_ref->{css}           = $ssi_path . 'meta-css.html';
	unless ($page_ref->{meta_keywords}) {
		my $keywords;
		if ($style eq 'paleo') {
			$keywords = "paleovirology,paleovirus,endogenous,ERV,EVE,ERVs,";
			$keywords   .= "virology,retrovirus,virus,fossil,genomics,";
			$keywords   .= "evolution,evolutionary,history,genomes";
		}
		if ($style eq 'fossil') {
			$keywords = "virus,fossil,genes,genomes,retrovirus,fossil record,";
			$keywords   .= "class,family,genus,species,host";
		}
		if ($style eq 'lab' or $style eq 'screening') {
			$keywords = "research,paleovirology,Gifford,bioinformatics,evolution";
		}
		if ($style eq 'tool') {
			$keywords = "bioinformatics,sequence,evolution,FASTA,alignment,tools,virology,";
			$keywords   .= "MSA,online,service,RESTful";
		}
		$page_ref->{meta_keywords} = $keywords;
	}	

	unless ($page_ref->{description}) {
		my $description = "Short articles describing paleovirological research";
		$description .= " in the Gifford lab, and discussing recent advances in";
		$description .= " the field of paleovirology at large.";
		$page_ref->{description}  = $description;
	}

	# Create description, set header image & title
	if ($style eq 'paleo') {
		unless ($page_ref->{title}) {
			$page_ref->{title} = 'The paleovirology website: home of the viral fossil record';
		}
		$page_ref->{head_image}   = $ssi_path . 'headimage_paleo.html';
		if ($level eq 1) {
			$page_ref->{main_path}    = $site_path . 'main/paleo/';
		}
		elsif ($level eq 2) {
			$page_ref->{main_path}    = $site_path . 'main/paleo/';
		}
	}
	elsif ($style eq 'lab') {
		if ($page eq 'gifford') {
			$page_ref->{title} = 'Paleovirology online: Robert J. Gifford';
		}
		else {
			$page_ref->{title} = 'Paleovirology online: Gifford Lab';
		}
		$page_ref->{head_image}   = $ssi_path . 'headimage_lab.html';
		$page_ref->{main_path}    = $site_path . 'main/lab/';
	}
	elsif ($style eq 'fossil') {
		unless ($page_ref->{title}) {
			$page_ref->{title} = 'Paleovirology online: The viral fossil record';
		}
		$page_ref->{head_image}   = $ssi_path . 'headimage_fossil.html';
		$page_ref->{main_path}    = $site_path . 'main/db_pages/';
	}
	elsif ($style eq 'tool') {
		unless ($page_ref->{title}) {
			$page_ref->{title} = 'Paleovirology online: Tools and services';
		}
		$page_ref->{head_image}   = $ssi_path . 'headimage_tool.html';
		$page_ref->{main_path}    = $site_path . 'main/tool/';
	}
	elsif ($style eq 'external') {
		$page_ref->{title} = 'Paleovirology online: External links';
		$page_ref->{head_image}   = $ssi_path . 'headimage_external.html';
		$page_ref->{main_path}    = $site_path . 'main/external/';
	}
	elsif ($style eq 'research') {
		unless ($page_ref->{title}) {
			$page_ref->{title} = 'Paleovirology online: Research from the Gifford Lab';
		}
		$page_ref->{head_image}   = $ssi_path . 'headimage_lab.html';
		if ($level eq 1) {
			$page_ref->{main_path}    = $site_path . 'main/lab/';
			$page_ref->{head_image}   = $ssi_path . 'headimage_paleo.html';
		}
		if ($level eq 2) {
			$page_ref->{main_path}    = $site_path . 'main/lab/';
		}
		elsif ($level eq 3) {
			$page_ref->{main_path}    = $site_path . 'main/posts/';
			$page_ref->{head_image}   = $ssi_path . 'headimage_paleo.html';
		}
	}
	elsif ($style eq 'screening') {
		unless ($page_ref->{title}) {
			$page_ref->{title} = 'Paleovirology online: Genome screening';
		}
		$page_ref->{head_image}   = $ssi_path . 'headimage_screening.html';
		$page_ref->{main_path}    = $site_path . 'main/screening/';
	}
	else { die; }

	# Set outpaths
	if ($level eq 1) {
		$page_ref->{outpath} = "index.html";
	}
	elsif ($level eq 2) {
		$page_ref->{outpath} = $site_path . 'html/' . "$page.html";
	}
	elsif ($level eq 3) {
		if ($style eq 'research') {
			$page_ref->{outpath} = $site_path . 'html/posts/' . "$page.html";
		}
		elsif ($style eq 'fossil') {
			$page_ref->{outpath} = $site_path . 'html/fossil-record/' . "$page.html";
		}
		elsif ($style eq 'tool') {
			my $report_dir = $site_ref->{report_dir};
			unless ($report_dir)  { die; }
			
			$page_ref->{outpath}   = $report_dir . "$page.html";
			$page_ref->{main_path} = $report_dir;
		}
		else { die; }
	}
	else { 
		print "$level"; 
		die; 
	}

	# Set main content ('left panel')
	my $main_path = $page_ref->{main_path};
	unless ($main_path) { die; }
	$page_ref->{left} = $main_path . $page . '.html';

	# Store the page data under its name
	my $pages_hash_ref = $site_ref->{pages_hash};
	if ($pages_hash_ref) {
		$pages_hash_ref->{$page} = $page_ref;
	}
	else {
		my %hash;
		$hash{$page} = $page_ref;
		$site_ref->{pages_hash} = \%hash;
	}

	# Store in an array
	my $pages_array_ref = $site_ref->{pages_array};
	if ($pages_hash_ref) {
		push (@$pages_array_ref, $page_ref);
	}
	else {
		my @array;
		push (@array, $page_ref);
	}

}

############################################################################
# CREATING BLAST LIBRARIES FROM REFERENCE LIBRARY
############################################################################

#***************************************************************************
# Subroutine:  clone_site  
# Description: clone the website 
# TODO
#***************************************************************************
sub clone_site {

	my ($self) = @_;
	
	#print "\n\n\t ### WARNING - you must be root for all changes to take effect\n\n"; 

	# Name of the new site
	my $question = "\tEnter the name of the site";
	my $site = $console->ask_question($question);
	chomp $site;

	# Create the new site directory
	my $command = "mkdir $site";
	system $command;

	# Copy what we need 
	$command = "cp -r bin $site/";
	print "\n\t$command";
	system $command;
	$command = "cp -r db $site/";
	print "\n\t$command";
	system $command;
	$command = "mkdir $site/library";
	print "\n\t$command";
	system $command;
	$command = "cp -r library/Base $site/library/";
	print "\n\t$command";
	system $command;
	$command = "cp -r library/Component $site/library";
	print "\n\t$command";
	system $command;
	$command = "mkdir $site/library/Program";
	print "\n\t$command";
	system $command;
	$command = "cp -r library/Program/GLUE/ $site/library/Program/";
	print "\n\t$command";
	system $command;
	$command = "cp -r library/Program/System/ $site/library/Program/";
	print "\n\t$command";
	system $command;
	$command = "cp -r site $site/";
	print "\n\t$command";
	system $command;
	$command = "rm -rf $site/site/reports/*";
	print "\n\t$command";
	system $command;
	$command = "chown apache $site/site/reports/*";
	print "\n\t$command";
	system $command;

	$command = "cp *.html $site/";
	print "\n\t$command";
	system $command;
	print "\n\t$command";
	$command = "cp system.pl $site/";
	system $command;
	print "\n\t$command";
	$command = "cp glue.cgi $site/";
	print "\n\t$command";
	system $command;

	# Apply the site name to ssi
	my $ssi_dir = "$site/site/ssi/";
	$fileio->replace_phrase_in_all_files($ssi_dir, 'vglue_sandbox', $site);
	# Apply the site name to main posts
	my $post_dir = "$site/site/main/posts/";
	$fileio->replace_phrase_in_all_files($post_dir, 'vglue_sandbox', $site);
	# Apply the site name to main tool dir
	my $tool_dir = "$site/site/main/tool/";
	$fileio->replace_phrase_in_all_files($tool_dir, 'vglue_sandbox', $site);
	# Apply the site name to main fossil dir
	my $fossil_dir = "$site/site/main/fossil/";
	$fileio->replace_phrase_in_all_files($fossil_dir, 'vglue_sandbox', $site);
	# Apply the site name to main lab dir
	my $lab_dir = "$site/site/main/lab/";
	$fileio->replace_phrase_in_all_files($lab_dir, 'vglue_sandbox', $site);
	# Apply the site name to main paleo dir
	my $paleo_dir = "$site/site/main/paleo/";
	$fileio->replace_phrase_in_all_files($paleo_dir, 'vglue_sandbox', $site);

	# Apply the site name to View.pm
	my @views;
	my $view_path = "$site/library/Program/System/Views.pm";
	$fileio->read_input_file($view_path, \@views);
	my @view_switch;
	foreach my $line (@views) {
		$line =~ s/vglue_sandbox/$site/g;
		push (@view_switch, $line);
	}
	$fileio->write_output_file($view_path, \@view_switch);
         
	# Apply the site name to BLAST_Tool.pm
	my $blast_path = "$site/library/Program/GLUE/BLAST_Tool.pm";
	my @blast;
	$fileio->read_input_file($blast_path, \@blast);
	my @switch;
	foreach my $line (@blast) {
		$line =~ s/vglue_sandbox/$site/g;
		push (@switch, $line);
	}
	$fileio->write_output_file($blast_path, \@switch);
         
	# think thats it.
}

############################################################################
# SECTION: Consolidating reference sequence library
############################################################################

#***************************************************************************
# Subroutine:  consolidate library 
# Description: 
#***************************************************************************
sub consolidate_library {

	my ($self) = @_;

	# Load the retroviral reference sequences
	my %refseqs;
	my @refseqs;
	my $refseqlib_obj = RefSeqLibrary->new();
	my $libpath = $self->{refseq_lib_path};
	print "\n\t loading GLUE reference sequence library ' $libpath'";
	$refseqlib_obj->load_refseq_library($libpath, \@refseqs, \%refseqs);	
	
	my %ordered;
	$refseqlib_obj->order_by_family_and_genus(\@refseqs, \%ordered);
	my $family = 'Retroviridae';
	my $retroviridae_ref = $ordered{$family};

	# Iterate through the genera
	my $tmp_dir = './library_consolidate/';
	my @rv_genera = sort keys %$retroviridae_ref;
	foreach my $genus (@rv_genera) {
		
		# Output status and get the refseqs for this genus
		print "\n\t# Processing Genus '$genus'";
		my $genus_ref = $retroviridae_ref->{$genus};
		my @viruses = sort keys %$genus_ref;

		# Iterate through the viruses in this genus
		foreach my $virus (@viruses) {
			# Output status and get refseq data 
			my $refseq = $retroviridae_ref->{$genus}->{$virus};
			my $name = $refseq->{name};
			unless ($name eq $virus) {
				my $error = "\n\t ERROR - filename '$virus' ";
				$error .= "does not match refseq name '$name'\n\n";
			}
			print "\n\t # Validating refseq '$name'";
			$self->validate_retrovirus_refseq($refseq);
			my $tribe = $refseq->{virus_tribe};
			unless ($tribe) { die; }
			my $dir = $tmp_dir . $tribe . '/';
			$refseq->write_self_to_text($dir);
		}
	}
}

############################################################################
# SECTION: Command line title blurb 
############################################################################

#***************************************************************************
# Subroutine:  show_title
# Description: does what it says 
#***************************************************************************
sub show_title {

	$console->refresh();
	my $title       = 'System.pl';
	my $version     = '1.0';
	my $description = 'Manage GLUE-associated website and libraries';
	my $author      = 'Robert J. Gifford';
	my $contact		= '<rgifford@adarc.org>';
	$console->show_about_box($title, $version, $description, $author, $contact);
}

############################################################################
# BASIC FXNS
############################################################################

sub by_number { $a <=> $b }	

############################################################################
# EOF 
############################################################################
