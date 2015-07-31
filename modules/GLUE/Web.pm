#!/usr/bin/perl -w
############################################################################
# Module:      Web.pm
# Description: Module containing functions for creating HTML front ends to 
#              GLUE
# History:     January 2015: Created by Robert Gifford 
############################################################################
package Web;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;
use DBI;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::SeqIO;
use Base::FileIO;
use Base::DevTools;
use Base::Console;

use Interface::MySQLtable;

############################################################################
# Globals
############################################################################

# Create base objects
my $seqio     = SeqIO->new();
my $fileio    = FileIO->new();
my $devtools  = DevTools->new();
my $console   = Console->new();
my $writer    = HTML_Utilities->new();
my $seq_obj   = Sequence->new();

# Global site paths
my $site_name = 'paleo';                       
my $site_path = './site/';                       
my $main_path = './site/main/';                       
my $url       = "http://saturn.adarc.org/paleo/";

1;

#***************************************************************************
# Subroutine:  new
# Description: Parameters
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {
	
		# Paths and constants
		mode                  => $parameter_ref->{mode},
		process_id            => $parameter_ref->{process_id},
		output_type           => $parameter_ref->{output_type},
		
		# Paths and constants	
		refseq_use_path       => $parameter_ref->{refseq_use_path},
		db_list_path          => $parameter_ref->{db_list_path},
		html_path             => $parameter_ref->{html_path},
		html_header           => $parameter_ref->{html_header},
		html_header2          => $parameter_ref->{html_header2},
		html_footer2          => $parameter_ref->{html_footer2},
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# Public member functions
############################################################################

#***************************************************************************
# Subroutine:  retrieve_sequences 
# Description: 
#***************************************************************************
sub retrieve_sequences  {

	my ($self, $query,  $form) = @_;

	# Connection
	my %params;
	$params{server}   = 'localhost';
	$params{username} = 'root';    
	$params{password} = 'blenat2'; 
	my $db_obj  = SequenceDB->new(\%params);

	#$devtools->print_hash($query); die;
	# Set defaults
	my $lim_genotype  = $query->{'limit_by_genotype'};
	my $genotype      = $query->param('Genotype');
	my $lim_host      = $query->{'limit_by_host'};
	my $host          = $query->param('Host');
	my $lim_coverage  = $query->{'limit_by_coverage'};
	my $seq_start     = $query->param('seq_start');
	my $seq_end       = $query->param('seq_end');
	my $lim_country   = $query->{'limit_by_country'};
	my $country       = $query->param('Country');
	my $lim_year      = $query->{'limit_by_year'};
	my $year_rule     = $query->{'year_rule'};
	my $year_low      = $query->{'year_low'};
	my $year_high     = $query->{'year_high'};
	my $head_genotype = $query->{'head_genotype'};
	my $head_country  = $query->{'head_country'};
	my $head_date     = $query->{'head_date'};
	my $head_host     = $query->{'head_host'};
	my $database      = $query->param('database');
	$db_obj->load_sequence_db($database);

	# Create the SQL
	my $select = ' SELECT ';
	my @fields = qw [ Sequence.sequence_id Sequence.isolate_id 
                      Sequence.host Sequence.isolation_date 
                      Sequence.start Sequence.stop
					  Genotype.genotype	
					  Location.location_ID
                      Sequence.sequence
					];

	my $fields = join(",", @fields);
	$select .= $fields;
	my @tables = qw [ Sequence ];
	push (@tables, 'Genotype');
	push (@tables, 'Location');
	my $tables = join(',', @tables);
	$select .= " FROM $tables";

	my $where = " WHERE Sequence.sequence_id = Genotype.sequence_id 
                  AND   Sequence.sequence_id = Location.sequence_id ";
	if ($lim_country) {
		$where .= " AND Location_ID = '$country' AND Location_type = 'country' ";
	}
	if ($lim_genotype) {
		$where .= " AND Genotype = '$genotype' ";
	}
	if ($lim_host) {
		$where .= " AND Host = '$host' ";
	}
	if ($lim_coverage) {
		unless ($seq_start) { die; }
		unless ($seq_end) { die; }
		$where .= " AND Sequence.start >= $seq_start ";
		$where .= " AND Sequence.stop  >= $seq_end ";
	}
	


	$select .= $where;
	#print "\n\n\t\t $select <BR><BR>";

	my $dbh = $db_obj->{dbh};
	my $sth = $dbh->prepare($select);
	unless ($sth->execute()) { print $select; exit; }	
	my $row_count = 0;
	my $lowest_start = undef;
	my $highest_stop = 1;
	my @data;
	while (my $row = $sth->fetchrow_arrayref) {
		
		$row_count++;
		my $i = 0;
		my %row;
		my @header;
		my $sequence;
		my %data;
		foreach my $field (@fields) {
			my $value = @$row[$i];
			$data{$field} = $value;
			unless ($field eq 'Sequence.sequence'
                 or $field eq 'Sequence.start'
                 or $field eq 'Sequence.stop') {
				push (@header, $value);
			}
			if ($field eq 'Sequence.start') {
				unless ($lowest_start) {
					$lowest_start = $value;
				}
				elsif ($value < $lowest_start) {
					$lowest_start = $value;
				}
			}
			if ($field eq 'Sequence.stop') {
				if ($value > $highest_stop) {
					$highest_stop = $value;
				}
			}	
			$i++;
		}
		my $header = join ('|', @header);
		$data{header} = $header;
		push (@data, \%data);
	}
	#print "<BR> LOWEST  $lowest_start";	
	#print "<BR> HIGHEST $highest_stop";	

	my @fasta;
	push (@fasta, "<font='courier'>");
	foreach my $data_ref (@data) {

		my $header   = $data_ref->{'header'};
		my $start    = $data_ref->{'Sequence.start'};
		my $stop     = $data_ref->{'Sequence.stop'};
		my $sequence = $data_ref->{'Sequence.sequence'};
		unless ($sequence)  { die 'RETRIEVE FAILED'; }
		my $lead  = $start - $lowest_start;
		my $trail = $highest_stop - $stop;
		my $lead_pad  = '-' x $lead;
		my $trail_pad = '-' x $trail;
		$sequence = $lead_pad . $sequence . $trail_pad;	
		#print "<BR> Header $header";
		#print "<BR> Lead  pad $lead ($start - $lowest_start)";
		#print "<BR> TRail pad $trail ($highest_stop - $stop)";
		my $f_sequence;
		my @sequence = split ('', $sequence);
		my $i;
		foreach my $char (@sequence) {
			$i++;
			$f_sequence .= $char;
			if ($i eq 100) {
				$i = 0;
				$f_sequence .= "<BR>";
			}
		
		}
		my $fasta = ">$header\n$f_sequence\n";
		push (@fasta, $fasta);

	}
	#die;
	#print " $select<BR><BR>";

	my $sequences = join("\n\n", @fasta);
	if ($sequences) {
		# Print the formatted input page if no sequnce data has been received
		my @html;
		$fileio->read_file('./site/ssi/glue_head2.html', \@html);
		push (@html, "\n<pre>");
		my $html = join("\n", @html);
		print $html; 
		print $sequences; 
		print "\n</pre>";
	}
	else {
		die;
		# Print the formatted input page if no sequnce data has been received
		my @html;
		$fileio->read_file($form, \@html);
		my $html = join("\n", @html);
		print $html; die;
	}

}

#***************************************************************************
# Subroutine:   retrieve_mutations 
# Description: 
#***************************************************************************
sub retrieve_mutations  {

	my ($self, $query, $form) = @_;



	# Print the formatted input page if no sequnce data has been received
	my @html;
	$fileio->read_file('./site/ssi/glue_head2.html', \@html);
	my $html = join("\n", @html);
	print $html; 

	my @mutations;
    push (@mutations, "\n<h3>Mutation frequencies</h3>");
	push (@mutations, "\n<div class=divider></div>");

	# Write genes table
	my %table_set_ref;
	$table_set_ref{indent} = "\n\t\t";	
	$table_set_ref{class} = 'data';
	$table_set_ref{width} = 800;	
	my $table_open = $writer->open_table(\%table_set_ref);
	push (@mutations, $table_open);


	my @header;	
	$header[0] = 'Gene';
	$header[1] = 'Mutation type';
	$header[2] = 'Mutation';
	$header[3] = 'Reference AA';
	$header[4] = 'Position';
	$header[5] = 'AA';
	$header[6] = 'Frequency';
	my $head_row = $writer->create_row(\@header); 
	push (@mutations, "\n\t $head_row");


	my $mutations = join("\n", @mutations);
	print $mutations;

	my $table_close = $writer->close_table(\%table_set_ref);
	push (@mutations, $table_close);


	# Show form again if input is not complete
	#my @html;
#	$fileio->read_file($form, \@html);
	#my $html = join("\n", @html);
	#print $html; die;

}

#***************************************************************************
# Subroutine:  generate_db_retrieve_html
# Description: 
#***************************************************************************
sub generate_db_retrieve_html  {

	my ($self) = @_;

	# Connection
	my %params;
	$params{server}   = 'localhost';
	$params{username} = 'root';    
	$params{password} = 'blenat2'; 
	my $refseq_name; # TODO - how to get this?

	# Get the list of DB names from list of refseqs (note: assumes correspondence)
	my @sequence_dbs;
	my $db_list_path = $self->{db_list_path};
	unless ($db_list_path) { die; }
	$fileio->read_file($db_list_path, \@sequence_dbs);
	#$fileio->read_directory_to_array($db_list_path, \@sequence_dbs);
	#$devtools->print_array(\@sequence_dbs); die;
	
	# For each sequence database
	my %db_data;
	foreach my $seq_db_name (@sequence_dbs) {

		# Load the database
		chomp $seq_db_name;
		my $db_obj = SequenceDB->new(\%params);
		$db_obj->load_sequence_db($seq_db_name);
	
		# Get the parameter ranges for this database and refseq
		my %form_params;
		$form_params{db_name}     = $seq_db_name;
		$form_params{refseq_name} = $seq_db_name ;
		$self->get_db_param_ranges($db_obj, \%form_params);
		#$devtools->print_hash(\%form_params);	die;

		# Write web form page
		my %data;
		$self->write_db_retrieval_html(\%form_params, \%data);
		my $years_ref = $form_params{years};
		my $year_low = shift @$years_ref;
		my $year_high  = pop @$years_ref;
		$data{year_high} = $year_high;
		$data{year_low} = $year_low;
		$data{refseq_full_name} = $form_params{refseq_full_name};
	
		# Get db summary info
		$self->get_summary_info($db_obj, \%data);
		$db_data{$seq_db_name} = \%data;
		#$devtools->print_hash(\%db_data); die;
	}
	
	# Create and write the top-level page to access individual web-forms
	$self->write_database_access_page(\@sequence_dbs, \%db_data);

}

#***************************************************************************
# Subroutine:  get_summary_info
# Description: 
#***************************************************************************
sub get_summary_info {

	my ($self, $db_obj, $data_ref) = @_;

	# Get tables
	my $sequence_table  = $db_obj->{sequence_table};
	my $filter_table    = $db_obj->{filter_table};
	my $genotype_table  = $db_obj->{genotype_table};
	my $location_table  = $db_obj->{location_table};

	# Get sequence count
	my @seq_fields;
	push (@seq_fields, 'COUNT(*)');
	my @seq_count;
	$sequence_table->select_rows(\@seq_fields, \@seq_count);
	my $count_ref = shift @seq_count;
	my $total_seqs = $count_ref->{'COUNT(*)'};
	$data_ref->{total_seqs} = $total_seqs;

	# Get number of distinct genotypes
	my @genotype_fields = qw [ genotype ];
	my @genotypes;
	my $where = " WHERE genotype_method = 'Genbank'  ";
	$genotype_table->select_distinct(\@genotype_fields, \@genotypes, $where);
	#$devtools->print_array(\@genotypes); die;
	my $num_genotypes = scalar @genotypes;
	$data_ref->{total_genotypes} = $num_genotypes;

	# Get number of distinct hosts
	my @host_fields = qw [ host ];
	my @hosts;
	$sequence_table->select_distinct(\@host_fields, \@hosts);
	my $num_hosts = scalar @hosts;
	$data_ref->{total_hosts} = $num_hosts;

	# Get number of distinct countries
	my @country_fields = qw [ location_id ];
	my @countries;
	$where = " WHERE location_type = 'country'  ";
	$location_table->select_distinct(\@country_fields, \@countries, $where);
	my $num_countries = scalar @countries;
	$data_ref->{total_countries} = $num_countries;

	# Get year range


}

#***************************************************************************
# Subroutine:  write_database_access_page
# Description: dynamically create page with links to server Sequence DB forms
#***************************************************************************
sub write_database_access_page {

	my ($self, $databases_ref, $db_data_ref) = @_;

	# Create the page header
	my @html;
	my $html_header = $self->{html_header2};
	$fileio->read_file($html_header, \@html);

	push (@html, "\n\t<div class='left'>");
	push (@html, "\n\t<h4>Databases</h4>");
	push (@html, "\n\t<BR><BR>");
	push (@html, "\n\t\t<TABLE class=data width=840>");


	push (@html, "\n\t\t<tr>\n\t\t");
	push (@html, "\n\t\t <td>");
	my @header;
	push (@header, '<b>Database</b>');
	push (@header, '<b>Virus name</b>');
	push (@header, '<b># Sequences</b>');
	push (@header, '<b># Genotypes</b>');
	push (@header, '<b># Host sp.</b>');
	push (@header, '<b># Countries</b>');
	push (@header, '<b># Year range</b>');
	my $header = join("\n\t\t </td>\n\t\t <td>", @header);	
	push (@html, $header);
	push (@html, "\n\t\t </td>");
	push (@html, "\n\t\t</tr>\n\n");

	foreach my $database (@$databases_ref) {

		push (@html, "\n\t\t<tr>\n\t\t");
		push (@html, "\n\t\t <td>");

		# Get totals
		my $data_ref        = $db_data_ref->{$database};
		#$devtools->print_hash($data_ref); die;

		unless ($data_ref) { die; };	
		my $form_path       = $data_ref->{form_path};
		my $full_name       = $data_ref->{refseq_full_name};
		my $total_seqs      = $data_ref->{total_seqs};
		my $total_genotypes = $data_ref->{total_genotypes};
		my $total_hosts     = $data_ref->{total_hosts};
		my $total_countries = $data_ref->{total_countries};
		my $year_high       = $data_ref->{year_high};
		my $year_low        = $data_ref->{year_low};
		#unless ($year_low and $year_high) { die; }
		
		my @line;
		my $f_database = "<a href=\"$form_path\">$database</a></td><td>$full_name</td>";
		push (@line, $f_database);
		push (@line, $total_seqs);
		push (@line, $total_genotypes);
		push (@line, $total_hosts);
		push (@line, $total_countries);
		my $range = "$year_low-$year_high";
		push (@line, $range);
		my $line = join("\n\t\t </td>\n\t\t <td>", @line);	
		push (@html, $line);

		push (@html, "\n\t\t </td>");
		push (@html, "\n\t\t</tr>\n\n");
	}	
	push (@html, "\n\t\t</TABLE>");

	# Create the page footer
	push (@html, "\n\t<br><br> The above table shows a list of currently active GLUE sequence databases on this server");
	push (@html, "\n\t<br><br> Click the link on the left (in the database name column) to access the data retrieval form for individual databases");
	push (@html, "\n\t</div>");
	my @footer;
	my $html_footer = $self->{html_footer2};
	$fileio->read_file($html_footer, \@footer);
	push (@html, @footer);

	# Write the HTML file
	$fileio->write_file('./site/html/data.html', \@html);
}

#***************************************************************************
# Subroutine:  write_db_retrieval_html
# Description: 
#***************************************************************************
sub write_db_retrieval_html {

	my ($self, $form_params, $data_ref) = @_;

	# Create the page header
	my $html_header = $self->{html_header2};
	my @html;
	$fileio->read_file($html_header, \@html);
	
	# Set the path to HTML
	my $db_name = $form_params->{db_name};
	my $refseq_full_name = $form_params->{refseq_full_name};
	my $html_path = $self->{html_path};
	my $output_file = $db_name . '_retrieve.html';
	my $output_html = $html_path . $output_file;
	
	# Create the web form part
	my @form;
	push (@form, "\n\t<div class='left'>");
	push (@form, "\n\t<h4>Retrieve aligned sequence data from $db_name");
	push (@form, " ($refseq_full_name) database</h4>");
	push (@form, "\n\t<div class=\"divider2\"></div>");
	push (@form, "\n\t<form method=\"post\" ");
	push (@form, "action=\"http://saturn.adarc.org/GLUE-dev/glue.cgi\" ");
	push (@form, "enctype=\"multipart/form-data\" id=\"resform\" target=\"_blank\">\n");
	$self->write_sequence_form(\@form, $form_params);	
	push (@form, "\n\t<br><input type=\"submit\" value=\"Retrieve\">");
	#push (@form, "\n\t<br><input type=\"hidden\" name=\"retrieve_mode\" value=\"sequence\">");
	push (@form, "\n\t<br><input type=\"hidden\" name=\"mode\" value=\"SEQUENCE_RETRIEVE\">");
	push (@form, "\n\t<br><input type=\"hidden\" name=\"database\" value=\"$db_name\">");
	push (@form, "\n\t<br><input type=\"hidden\" name=\"report_style\" value=\"GLUE\">");
	push (@form, "\n\t<br><input type=\"hidden\" name=\"retrieve_form\" value=\"$output_html\">");
	push (@form, "\n\t</form>");
           
	push (@form, "\n\t<br><br>");
	push (@form, "\n\t<h4>Retrieve mutation data from $db_name database</h4>");
	push (@form, "\n\t<div class=\"divider2\"></div>");
	push (@form, "\n\t<form method=\"post\" ");
	push (@form, "action=\"http://saturn.adarc.org/GLUE-dev/glue.cgi\" ");
	push (@form, "enctype=\"multipart/form-data\" id=\"resform\">\n");
	$self->write_mutation_form(\@form, $form_params);	
	push (@form, "\n\t<br><input type=\"submit\" value=\"Retrieve\">");
	#push (@form, "\n\t<br><input type=\"hidden\" name=\"retrieve_mode\" value=\"mutation\">");
	push (@form, "\n\t<input type=\"hidden\" name=\"mode\" value=\"MUTATION_RETRIEVE\">");
	push (@form, "\n\t<input type=\"hidden\" name=\"report_style\" value=\"GLUE\">");
	push (@form, "\n\t<input type=\"hidden\" name=\"retrieve_form\" value=\"$output_html\">");
	push (@form, "\n\t</form>");
	push (@html, @form);

	# Create the page footer
	push (@form, "\n\t</div>");
	my @footer;
	my $html_footer = $self->{html_footer2};
	$fileio->read_file($html_footer, \@footer);
	push (@html, @footer);


	# Write the HTML file
	print "\n\t Writing data retrieval page for DB $db_name to '$output_html'\n";
	$fileio->write_file($output_html, \@html);

	$data_ref->{form_path} = $output_file;
}

#***************************************************************************
# Subroutine:  write mutation form
# Description: 
#***************************************************************************
sub write_mutation_form {

	my ($self, $html_ref, $form_params) = @_;

	my $refseq_len = $form_params->{refseq_len};
	my @form;
	push (@form, "\n\t\t<TABLE class=data>");
	
	# Add genotype select options
	push (@form, "\n\t\t<tr>\n\t\t <td>");
	push (@form, "\n\t\t  Get mutation frequencies in gene");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t <td> <select name=\"Gene\">");
	my $genes = $form_params->{gene_names};
	foreach my $gene (@$genes) {
		my $line = "\n\t\t  <option value=\"$gene\"> $gene</option>";
		push (@form, $line);
	}
	push (@form, "\n\t\t </select>");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t</tr>");

	push (@form, "\n\t\t<tr>\n\t\t <td>");
	push (@form, "\n\t\t  <INPUT TYPE=CHECKBOX NAME=\"limit_by_genotype\" VALUE=\"\" >");
	push (@form, "\n\t\t  Limit by genotype:");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t <td>");
	push (@form, "\n\t\t <select name=\"Genotype\">");
	my $genotypes = $form_params->{genotypes};
	foreach my $genotype (@$genotypes) {
		my $line = "\n\t\t  <option value=\"$genotype\"> $genotype</option>";
		push (@form, $line);
	}
	push (@form, "\n\t\t </select>");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t</tr>\n\n");

	# Add select by country options
	push (@form, "\n\t\t<tr>\n\t\t <td>");
	push (@form, "\n\t\t  <INPUT TYPE=CHECKBOX NAME=\"mut_fiuysuyeru\" VALUE=\"\" >");
	push (@form, "\n\t\t  Limit by country:");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t <td>");
	push (@form, "\n\t\t <select name=\"Country\">");
	my $countries = $form_params->{countries};
	foreach my $country (@$countries) {
		my $line = "\n\t\t  <option value=\"$country\"> $country</option>";
		push (@form, $line);
	}
	push (@form, "\n\t\t </select>");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t</tr>");

	# Add year select options
	push (@form, "\n\t\t<tr>\n\t\t <td>");
	push (@form, "\n\t\t  <INPUT TYPE=CHECKBOX NAMEfiuysuyeru\" VALUE=\"\" >");
	push (@form, "\n\t\t  Limit by year:");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t <td>");
	push (@form, "\n\t\t <select name=\"year_rule\">");
	push (@form, "\n\t\t <option selected value=\"include\"> include</option>");
	push (@form, "\n\t\t <option value=\"exclude\"> exclude</option>");
	push (@form, "\n\t\t </select>");
	push (@form, "\n\t\t &nbsp; sequences in the year ange");
	push (@form, "\n\t\t </td>");

	push (@form, "\n\t\t <td>");
	push (@form, "\n\t\t <select name=\"year_low\">");
	my $years = $form_params->{years};
	foreach my $year (@$years) {
		my $line = "\n\t\t  <option value=\"$year\"> $year</option>";
		push (@form, $line);
	}
	push (@form, "\n\t\t </select>");
	push (@form, "\n\t\t </td>");

	push (@form, "\n\t\t <td>");
	push (@form, "\n\t\t to <select name=\"year_high\">");
	foreach my $year (@$years) {
		my $line = "\n\t\t  <option value=\"$year\"> $year</option>";
		push (@form, $line);
	}
	push (@form, "\n\t\t </select>");
	push (@form, "\n\t\t </td>");

	push (@form, "\n\t\t</tr>\n\n");
	push (@form, "\n\t\t</TABLE>");

	push (@$html_ref, @form);
}

#***************************************************************************
# Subroutine:  write sequence form
# Description: 
#***************************************************************************
sub write_sequence_form {

	my ($self, $html_ref, $form_params) = @_;

	my $refseq_len = $form_params->{refseq_len};
	
	# Add genotype select options
	my @form;
	push (@form, "\n\t\t<TABLE class=data>");
	push (@form, "\n\t\t<tr>");
	push (@form, "\n\t\t <td>");
	push (@form, "\n\t\t  <INPUT TYPE=CHECKBOX NAME=\"limit_by_genotype\" VALUE=\"\" CHECKED=\"checked\">");
	push (@form, "\n\t\t  Limit by genotype:");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t <td>");
	push (@form, "\n\t\t <select name=\"Genotype\">");
	my $genotypes = $form_params->{genotypes};
	foreach my $genotype (@$genotypes) {
		my $line = "\n\t\t  <option value=\"$genotype\"> $genotype</option>";
		push (@form, $line);
	}
	push (@form, "\n\t\t </select>");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t</tr>");

	push (@form, "\n\t\t<tr>");
	push (@form, "\n\t\t <td>");
	push (@form, "\n\t\t  <INPUT TYPE=CHECKBOX NAME=\"limit_by_host\" VALUE=\"\" CHECKED=\"checked\">");
	push (@form, "\n\t\t  Limit by host:");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t <td>");
	push (@form, "\n\t\t <select name=\"Host\">");
	my $hosts = $form_params->{hosts};
	foreach my $host (@$hosts) {
		my $line = "\n\t\t  <option value=\"$host\"> $host</option>";
		push (@form, $line);
	}
	push (@form, "\n\t\t </select>");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t</tr>");



	# Add select by country options
	push (@form, "\n\t\t<tr>\n\t\t <td>");
	push (@form, "\n\t\t <INPUT TYPE=CHECKBOX NAME=\"limit_by_country\" VALUE=\"\">");
	push (@form, "\n\t\t Limit by country");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t <td>");
	push (@form, "\n\t\t <select name=\"Country\">");
	my $countries = $form_params->{countries};
	foreach my $country (@$countries) {
		my $line = "\n\t\t  <option value=\"$country\"> $country</option>";
		push (@form, $line);
	}
	push (@form, "\n\t\t </select>");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t</tr>");

	# Add year select options
	push (@form, "\n\t\t<tr>\n\t\t <td>");
	push (@form, "\n\t\t  <INPUT TYPE=CHECKBOX NAME=\"limit_by_year\" VALUE=\"\" >");
	push (@form, "\n\t\t  Limit by year:");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t <td>");
	push (@form, "\n\t\t <select name=\"year_rule\">");
	push (@form, "\n\t\t <option selected value=\"include\"> include</option>");
	push (@form, "\n\t\t <option value=\"exclude\"> exclude</option>");
	push (@form, "\n\t\t </select>");
	push (@form, "\n\t\t  sequences in the year range");
	push (@form, "\n\t\t <select name=\"year_low\">");
	my $years = $form_params->{years};
	foreach my $year (@$years) {
		my $line = "\n\t\t  <option value=\"$year\"> $year</option>";
		push (@form, $line);
	}
	push (@form, "\n\t\t </select>");
	push (@form, "\n\t\t to <select name=\"year_high\">");
	foreach my $year (@$years) {
		my $line = "\n\t\t  <option value=\"$year\"> $year</option>";
		push (@form, $line);
	}
	push (@form, "\n\t\t </select>");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t</tr>\n\n");


	push (@form, "\n\t\t<tr>");
	push (@form, "\n\t\t <td>");
	push (@form, "\n\t\t <INPUT TYPE=CHECKBOX NAME=\"limit_by_coverage\" VALUE=\"\">");
	push (@form, "\n\t\t Limit by coverage");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t <td>");
	push (@form, "\n\t\t from NT position ");
	push (@form, "\n\t\t <input type=\"text\" name=\"seq_start\" value=\"1\" size=\"4\" maxlength=\"5\" />");
	push (@form, "\n\t\t to");
	push (@form, "\n\t\t <input type=\"text\" name=\"seq_end\" value=\"$refseq_len\" size=\"4\" maxlength=\"5\"/>");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t</tr>\n\n");






	
	push (@form, "\n\t\t <tr>");
	push (@form, "\n\t\t <td>");
	push (@form, "\n\t\t  <INPUT TYPE=CHECKBOX NAME=\"limit_by_mutation\" VALUE=\"\" >");
	push (@form, "\n\t\t  Limit by mutation:");
	push (@form, "\n\t\t </td>");

	push (@form, "\n\t\t <td> in gene <select name=\"Gene\">");
	my $genes = $form_params->{gene_names};
	foreach my $gene (@$genes) {
		my $line = "\n\t\t  <option value=\"$gene\"> $gene</option>";
		push (@form, $line);
	}
	push (@form, "\n\t\t </select>");
	push (@form, ' at position <input type="text" name="mut_position" value="1" size="4" maxlength="5" />');
	push (@form, "\n\t\t to amino acid ");
	push (@form, "\n\t\t <select name=\"Amino\">");
	my @aminos = qw [ A C D E F G H I K L M N P Q R S T V W Y ];
	foreach my $amino (@aminos) {
		my $line = "\n\t\t  <option value=\"$amino\"> $amino</option>";
		push (@form, $line);
	}
	push (@form, "\n\t\t </select>");
	push (@form, "\n\t\t </td>");
	push (@form, "\n\t\t</tr>\n\n");



	push (@form, "\n\t\t</TABLE>");

	push (@$html_ref, @form);
}

#***************************************************************************
# Subroutine:  get_db_param_ranges
# Description: 
#***************************************************************************
sub get_db_param_ranges  {

	my ($self, $db_obj, $data_ref) = @_;

	# get tables
	my $sequence_table  = $db_obj->{sequence_table};
	my $filter_table    = $db_obj->{filter_table};
	my $genotype_table  = $db_obj->{genotype_table};
	my $location_table  = $db_obj->{location_table};

	# Get the range of year dates
	my @dates;
	my @date_fields = qw [ isolation_date ];
	my $where = " WHERE isolation_date != 'NULL'";
	$sequence_table->select_distinct(\@date_fields, \@dates, $where);
	my %years;
	foreach my $row_ref (@dates) {
		my $date = $row_ref->{isolation_date};
		my @data = split('-', $date);
		my $year = pop @data;
		my $length = length $year;
		if ($length eq 2) {
			if ($year < 20) {
				$year = 20 . $year;
			}
			else {
				$year = 19 . $year;
			}
		}
		$years{$year} = 1;
	}
	my @years = sort by_number keys %years;

	# Get the range of genotypes
	my @geno_data;
	my @geno_fields = qw [ genotype ];
	my $order = " ORDER BY LENGTH(genotype), genotype";
	$genotype_table->select_distinct(\@geno_fields, \@geno_data, $order);
	my %genotypes;
	my @genotypes;
	foreach my $row_ref (@geno_data) {
		my $genotype = $row_ref->{genotype};
		unless ($genotypes{$genotype}) {
			$genotypes{$genotype} = 1;
			push (@genotypes, $genotype);
		}
	}

	# Get the range of countries
	my @country_rows;
	my @geo_fields = qw [ location_id ];
	$where = " WHERE Location_type = 'country' ";
	$location_table->select_distinct(\@geo_fields, \@country_rows, $where);
	my %countries;
	foreach my $row_ref (@country_rows) {
		my $country = $row_ref->{location_id};
		$countries{$country} = 1;
	}
	my @countries = sort keys %countries;

	# Get the range of hosts
	my @host_rows;
	my @host_fields = qw [ host ];
	my $host_order = " ORDER BY Host ";
	$sequence_table->select_distinct(\@host_fields, \@host_rows, $host_order);
	my %hosts;
	my @hosts;
	foreach my $row_ref (@host_rows) {
		my $host = $row_ref->{host};
		unless ($hosts{$host}) {
			$hosts{$host} = 1;
			push (@hosts, $host);
		}
	}
	
	# Get the refseq parameters	
	my $refseq_name     = $data_ref->{refseq_name};
	my $refseq_use_path = $self->{refseq_use_path};
	unless ($refseq_use_path) { die; }
	my $refseq_path = $refseq_use_path . $refseq_name;
	my $parser_obj = RefSeqParser->new();
	my %params;
	$parser_obj->parse_refseq_flatfile($refseq_path, \%params);
	my $refseq = RefSeq->new(\%params);
	my $refseq_full_name = $refseq->{full_name};
	unless ($refseq_full_name) { 
		$refseq->describe();
		die; 
	}

	my $sequence = $refseq->{sequence};
	my $refseq_len = length $sequence;	

	my @gene_names;
	my $genes_ref = $refseq->{genes};
	foreach my $gene_ref (@$genes_ref) {
		my $gene_name = $gene_ref->{name};
		push (@gene_names, $gene_name);
	}

	# Get the polymorphism genes
	$data_ref->{countries}  = \@countries;
	$data_ref->{hosts}      = \@hosts;
	$data_ref->{genotypes}  = \@genotypes;
	$data_ref->{years}      = \@years;
	$data_ref->{gene_names} = \@gene_names;
	$data_ref->{refseq_len}       = $refseq_len;
	$data_ref->{refseq_full_name} = $refseq_full_name;

}

#***************************************************************************
# Subroutine:  create_refseq_view	
# Description: Create the main (left) panel content for the refseq HTML
#***************************************************************************
sub create_refseq_view {

	my ($self, $refseq, $html_ref) = @_;

	# Get refseq data
	my $indent = ' ' x 10;
	my $refseq_name  = $refseq->{name};
	unless ($refseq_name) { die; }

	my %orfs;
	$refseq->get_translated_orfs(\%orfs);
	$refseq->{translated_orfs} = \%orfs;

	# Write the taxonomy section
    #$self->write_retrovirus_header($refseq, $html_ref);

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
	push (@$html_ref, "\n$indent <a href=\"http://saturn.adarc.org/paleo/site/html/RV_refseq_definitions.html\"><u>here</u></a>.<br><br>");
	
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
	my %table_set_ref;
	#$table_set_ref{indent} = "\n\t\t";	
	$table_set_ref{class} = 'DATA';
	$table_set_ref{width} = 800;	
	#my $table_open = $writer->open_table(\%table_set_ref);
	#push (@$html_ref, $table_open);

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
# Subroutine:  by number
# Description: by number - for use with perl 'sort'  (cryptic but works) 
#***************************************************************************
sub by_number { $a <=> $b }	

############################################################################
# EOF
############################################################################
