#!/usr/bin/perl -w
############################################################################
# Module:      RefSeqTool.pm
# Description: GLUE-associated FASTA tool
# History:     January 2012: Created by Robert Gifford 
############################################################################
package RefSeqTool;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::BioIO;
use Base::FileIO;
use Base::SeqIO;
use Base::DevTools;
use Base::Console;
use Base::HTML_Utilities;
use Base::Sequence; # For performing basic sequence manipulations

############################################################################
# Globals
############################################################################

# Create base objects
my $fileio     = FileIO->new();
my $seqio      = SeqIO->new();
my $bioio      = BioIO->new();
my $devtools   = DevTools->new();
my $console    = Console->new();

my $header_limit  = 50;
my $default_delimiter =  "\t";  # default for delimited files
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
		output_type          => $parameter_ref->{output_type},
		
		# Member classes
		blast_obj 			 =>	$parameter_ref->{blast_obj},
		
		# Paths
		output_path          => $parameter_ref->{output_path},
		header_path          => $parameter_ref->{header_path},
		refseq_use_path      => $parameter_ref->{refseq_use_path}, 
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# SECTION: Define reference sequence core (online version) 
############################################################################

#***************************************************************************
# Subroutine:  define_refseq
# Description: define a reference sequence
#***************************************************************************
sub define_refseq {

	my ($self, $query) = @_;

	# Get params
	my $stage = 1;
	my %refseq_data;
	$self->get_refseq_form_params($query, \%refseq_data);

	# Create the reference sequence and write it to the report directory
	my $process_id  = $self->{process_id};
	my $output_path = $self->{output_path};
	my $report_dir  = $output_path . $process_id . '/';
	my $mkdir_cmd   = "mkdir $report_dir";
	system $mkdir_cmd;
	
	my $sequence     = $refseq_data{sequence};
	my $metadata_ref = $refseq_data{metadata};
	my $refseq_name  = $metadata_ref->{name};
	unless ($refseq_name) { 
		$devtools->print_hash($metadata_ref); die; 
	}
	my $fasta  = ">$refseq_name\n$sequence\n\n";
	my $raw_path = $report_dir . 'raw.fa';
	$fileio->write_text_to_file($raw_path, $fasta);
	# Convert mac line breaks
	my $command = "perl -pi -e 's/\r/\n/g' $raw_path";
	system $command;
	my @fasta;
	$seqio->read_fasta($raw_path, \@fasta);
	my $seq_ref = shift @fasta;
	$sequence = $seq_ref->{sequence};
	
	# Get the params from the form
	my $seq_type     = $refseq_data{sequence_type};
	my $eve_type     = $refseq_data{eve_type};
	my $genome_type  = $refseq_data{genome_type};
	my $genes_ref    = $refseq_data{genes};
	$refseq_data{sequence} = $sequence;
	my $refseq = RefSeq->new(\%refseq_data);
	$refseq->write_self_to_text($report_dir);
	
	# Create the main block (left panel) for the refseq HTML
	my $views_obj = Views->new($self);
	# First create the form for adding further annotation
	my @refseq_main;
	my %data;
	$data{refseq_name} = $refseq->{name};
	$data{report_dir}  = $report_dir;
	$self->create_form_part(\@refseq_main, \%data);

	# Then add refseq HTML below
	$views_obj->create_refseq_view($refseq, \@refseq_main);
	#$devtools->print_array(\@refseq_main);
	# Write the main block (left panel) to the report directory
	my $page_path = $report_dir . $refseq_name . '.html'; 
	#print "\n\t page path $page_path";
	$fileio->write_output_file($page_path, \@refseq_main);
    # Define the page in the correct format and write it
	my %site;
	$site{report_dir} = $report_dir;
	my %page;
	$page{title} = "Paleovirology online: defining GLUE reference sequence for $refseq_name";
	$views_obj->define_paleo_page(\%site, \%page, $refseq_name,  3, 'tool');
	my $pages_ref = $site{pages_hash};
	my $page_ref = $pages_ref->{$refseq_name};
	my @html;
	#$writer->assemble_page($page_ref, \@html);
	#my $output = join("\n", @html);
	#print $output;
}

#***************************************************************************
# Subroutine:  append_to_refseq
# Description: 
#***************************************************************************
sub append_to_refseq {

	my ($self, $query) = @_;

	# Get params
	my $stage = 1;
	my %form_data;
	$self->get_refseq_append_form_params($query, \%form_data);
	my $metadata_ref = $form_data{metadata};

	# Get the params from the form
	my $feature_type = $form_data{feature_type};
	my $report_dir   = $form_data{report_dir};
	my $append_path  = $form_data{append_path};
	unless ($append_path) { die "NO PATHS FOR APPEND"; }
	#$devtools->print_web_hash($refseq); die;

	# Parse the reference sequence and add new data
	my @refseq;
	my $refseq_parse_obj = RefSeqParser->new();
	my %params; 
	my $parser = RefSeqParser->new();
	$parser->parse_refseq_flatfile($append_path, \%params);
	my $refseq = RefSeq->new(\%params);
	my $refseq_name  = $refseq->{name};
	unless ($refseq_name) { die; }
	$refseq->append_new_data(\%form_data);
	$refseq->write_self_to_text($report_dir);
	#$devtools->print_web_hash($refseq);
	
	# Create the main block (left panel) for the refseq HTML
	my $views_obj = Views->new();

	# First create the form for adding further annotation
	my @refseq_main;
	my %data;
	$data{refseq_name} = $refseq->{name};
	$data{report_dir}  = $report_dir;
	$self->create_form_part(\@refseq_main, \%data);

	# Then add refseq HTML below
	$views_obj->create_refseq_view($refseq, \@refseq_main);
	# Write the main block (left panel) to the report directory
	my $page_path = $report_dir . $refseq_name . '.html'; 
	$fileio->write_output_file($page_path, \@refseq_main);
    # Define the page in the correct format and write it
	my @html;
	my %site;
	$site{report_dir} = $report_dir;
	my %page;
	$page{title} = "Paleovirology online: defining GLUE reference sequence for $refseq_name";
	$views_obj->define_paleo_page(\%site, \%page, $refseq_name,  3, 'tool');
	my $pages_ref = $site{pages_hash};
	my $page_ref = $pages_ref->{$refseq_name};
	#$writer->assemble_page($page_ref, \@html);

	# Write the main block (left panel) to the refseq HTML directory
	$fileio->write_output_file($page_path, \@html);
	my $output = join("\n", @html);
	print $output;
}

#***************************************************************************
# Subroutine:  create_form_part 
# Description: 
#***************************************************************************
sub create_form_part {
	
	my ($self, $output_ref, $hash_ref) = @_;

	my $refseq_name = $hash_ref->{refseq_name};
	my $report_dir  = $hash_ref->{report_dir};
	my $append_path = $report_dir . $refseq_name;

	my $html = '';
	$html .= '<h4>Annotate an existing EVE reference sequence</h4>';
	$html .= '<div class="separator"></div>';
	$html .= '<form method="post" action="http://saturn.adarc.org/vglue_sandbox/glue.cgi" 
				enctype="multipart/form-data" id="resform">
				<p>';
	
	# Reference sequence top line
	$html .= "\n\t\t<INPUT TYPE='hidden' NAME='mode' VALUE='REFSEQ_APPEND'>";
	$html .= "\n\t\t<INPUT TYPE='hidden' NAME='append_path' VALUE='$append_path'>";
	$html .= "\n\t\t<INPUT TYPE='hidden' NAME='report_dir'  VALUE='$report_dir'>";
	$html .= "\n\t\t<INPUT TYPE='hidden' NAME='family' VALUE='Retroviridae'>";
	
	#  Sequence feature
	$html .= '<br><b>Define a new sequence feature</b><br><br>';
	$html .= "Feature name &nbsp; <input type='text' name='feature_name' value='Pol' ";
	$html .= "size='4' maxlength='30' />";
	$html .= '&nbsp;&nbsp;';
	$html .= '<select name="feature_type">';
	$html .= "<option selected value='ORF'> ORF</option>";
	$html .= "<option value='UTR'> UTR</option>";
	$html .= '</select>&nbsp; &nbsp;';
	$html .= "Coordinates &nbsp; <input type='text' name='csv_coordinates' value='4744,7095'";
	$html .= ' size="18" maxleform/ngth="80"/><BR><BR>';	
	$html .= "\n\t\t<div><input type='hidden' name='action' value='ANALYZE'>";
	$html .= "\n\t\t<a href='javascript:document.forms[0].submit();'>";
	$html .= "\n\t\t<img src='http://saturn.adarc.org/vglue_sandbox/site/images/form/btn.png' alt='ANALYZE'/></a>";
	$html .= "\n\t<br></div></form><br>";
	$html .= '<div class="separator"></div>';	
	push (@$output_ref, $html);	

	my $link_text  = "<br> Click <u><a href='$append_path'>here</a></u>";
	$link_text    .= " to view the raw reference sequence <br><br>";
	push (@$output_ref, $link_text);	

}

############################################################################
# SECTION: Functions for getting values from web forms
############################################################################

#***************************************************************************
# Subroutine:  get_refseq_form_params 
# Description: get forms from the start page of the 'refseq define' process
#***************************************************************************
sub get_refseq_form_params {

	my ($self, $query, $data_ref) = @_;

	# Create writer utility obj
	my $writer = HTML_Utilities->new();	
	
	# Get pasted sequence data if its there
	my $referrer = 'vglue_configure';
	my %exons;
	my %feature;
	my %metadata;
	my $raw_data = $query->param('SequenceData');
	if ($raw_data) {
		
		# Create sequence
		my @dataset = split("\n", $raw_data); 
		my @f_dataset;
		foreach my $line (@dataset) { push (@f_dataset, "$line\n"); }
		my $sequence = join ('', @f_dataset);
		$sequence =~ s/\n//g;

		# Initialise CGI output
		my $append_path  = $query->param('append_path');
		my $report_dir   = $query->param('report_dir');
		$data_ref->{report_dir}    = $report_dir;
		$data_ref->{append_path}   = $append_path;
		$data_ref->{sequence}      = $sequence;
		$data_ref->{raw_data}      = \@f_dataset;
		
		# Metadata 
		my $supertribe = $query->param('tribe');
		my $tribe;
		if ($supertribe eq 'Alnilamretrovirinae') {
			$tribe = 'Alnilam';
		}
		elsif ($supertribe eq 'Alnitakretrovirinae') {
			$tribe = 'Alnitak';
		}
		elsif ($supertribe eq 'Mintakaretrovirinae') {
			$tribe = 'Mintaka';
		}
		else { die; }

		$metadata{name}             = $query->param('name');
		$metadata{full_name}        = $query->param('full_name');
		$metadata{family}           = $query->param('family');
		$metadata{virus_subfamily}  = $query->param('subfamily');
		$metadata{virus_supertribe} = $supertribe;
		$metadata{virus_tribe}      = $tribe;
		$metadata{virus_genus}      = $query->param('genus');
		$metadata{virus_subgroup}   = $query->param('subgroup');
		$metadata{host_sci_name}    = $query->param('host_sci_name');
		$metadata{host_common_name} = $query->param('host_common_name');
		$metadata{accession}        = $query->param('accession');
		$metadata{sequence_type}    = $query->param('sequence_type');
		$metadata{genome_type}      = $query->param('genome_type');
		$metadata{genome_coverage}  = $query->param('genome_coverage');
		
		# Features 
		my $feature_type = $query->param('feature_type');
		my $feature_name = $query->param('feature_name');
		my $coordinates  = $query->param('csv_coordinates');
		if ($feature_type and $feature_name and $coordinates) {
			
			my %starts;
			my @starts;
			my @coordinates = split(',', $coordinates);
			my $odd = undef;
			my $start = undef;
			my $coding_start;  # THIS MIGHT BE POINTLESS - need to consolidate
			my $first_start = undef;
			my $stop = undef;
			my $num_exons = 0;
			foreach my $position (@coordinates) {
				
				#print "<BR> $position";
				if ($odd) { 
					$stop = $position;
					$starts{$start} = $position;
					$start = undef;
					$num_exons++;
					$odd = undef;  
				}
				else { 
					$start = $position;
					push (@starts, $start);
					unless ($first_start) {
						$first_start = $start;
					}
					$odd = 'true'; 
					
				}	
			}
			if ($odd) { die; } # Not an even number of coordinates
			
			$feature{stop}         = $stop;
			$feature{coding_stop}  = $stop;
			$feature{starts}       = \@starts;
			$feature{start}        = $first_start;
			$feature{coding_start} = $first_start;
			$feature{type}         = $feature_type;
			$feature{name}         = $feature_name;
			$feature{full_name}    = $feature_name;
			$feature{num_exons}    = $num_exons;
			$feature{exons}        = \%starts;
			#$feature{csv_coordinates}  = $coordinates;
		}
	}
	# print the formatted input page if no sequence data has been received
	else {
		# Read  the VGLUE header
		my @html;
		$fileio->read_input_file('./site/html/refseq_define.html', \@html);
		my $html  = join('', @html);
		print $html;
		exit;
	}

	# Make the input structured the same as it would be from the refseq parser
	my @features;
	push (@features, \%feature);
	$data_ref->{genes}    = \@features;	
	$data_ref->{metadata} = \%metadata;	
}

#***************************************************************************
# Subroutine:  get_refseq_append_form_params
# Description: 
#***************************************************************************
sub get_refseq_append_form_params {

	my ($self, $query, $data_ref) = @_;

	# Create writer utility obj
	my $writer = HTML_Utilities->new();	
	
	# Get pasted sequence data if its there
	my %exons;
	my %feature;
	my %metadata;
	
	# Initialise CGI output
	my $append_path  = $query->param('append_path');
	my $report_dir   = $query->param('report_dir');
	my $genome_type  = $query->param('genome_type');
	my $seq_type     = $query->param('sequence_type');
	$data_ref->{report_dir}    = $report_dir;
	$data_ref->{append_path}   = $append_path;
	
	# Features 
	my $feature_type = $query->param('feature_type');
	my $feature_name = $query->param('feature_name');
	my $coordinates  = $query->param('csv_coordinates');
	if ($feature_type and $feature_name and $coordinates) {
		
		my %starts;
		my @starts;
		my @coordinates = split(',', $coordinates);
		my $odd = undef;
		my $start = undef;
		my $coding_start;  # THIS MIGHT BE POINTLESS - need to consolidate
		my $first_start = undef;
		my $stop = undef;
		my $num_exons = 0;
		foreach my $position (@coordinates) {
			
			#print "<BR> $position";
			if ($odd) { 
				$stop = $position;
				$starts{$start} = $position;
				$start = undef;
				$num_exons++;
				$odd = undef;  
			}
			else { 
				$start = $position;
				push (@starts, $start);
				unless ($first_start) {
					$first_start = $start;
				}
				$odd = 'true'; 
				
			}	
		}
		if ($odd) { die; } # Not an even number of coordinates
		
		$feature{stop}         = $stop;
		$feature{coding_stop}  = $stop;
		$feature{starts}       = \@starts;
		$feature{start}        = $first_start;
		$feature{coding_start} = $first_start;
		$feature{type}         = $feature_type;
		$feature{name}         = $feature_name;
		$feature{full_name}    = $feature_name;
		$feature{num_exons}    = $num_exons;
		$feature{exons}        = \%starts;
		#$feature{csv_coordinates}  = $coordinates;
	}
	else { # print the formatted input page if no sequence data has been received
		# TO DO FIX THIS
		die;
		# Read  the VGLUE header
		my $views_obj = Views->new();
		my $site_ref = $views_obj->{paleo_site};
		my @page;
		#$writer->assemble_page($site_ref, $referrer, \@page);
		print @page;
		exit;
	}

	# Make the input structured the same as it would be from the refseq parser
	my @features;
	push (@features, \%feature);
	$data_ref->{genes}    = \@features;	
	$data_ref->{metadata} = \%metadata;	
}

#***************************************************************************
# Subroutine:  by_number
# Description: by number - for use with perl 'sort'  (cryptic but works) 
#***************************************************************************
sub by_number { $a <=> $b }	

############################################################################
# EOF
############################################################################
