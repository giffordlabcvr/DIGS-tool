#!/usr/bin/perl -w
############################################################################
# Script:      HTML_Utilities 
# Description: Perl module containing utilities for producing HTML
# History:     (Rob Gifford) December 2006: Creation
# History:     (Rob Gifford) December 2013: Updated
############################################################################
package HTML_Utilities;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

use Base::FileIO;

############################################################################
# Globals
############################################################################

# Base objects
my $fileio      = FileIO->new();
my $devtools    = DevTools->new();

# Indent globals
my $indent_size = 2;
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create a new HTML_Utilies object 
# Arguments:   $raw_input: reference to an array containing raw fasta data 
#                          with headers as keys and sequnces as values
#***************************************************************************
sub new {

	my ($invocant, $parameters) = @_;
	my $class = ref($invocant) || $invocant;

	# Member variables
	my $self = {
	

	};

	bless ($self, $class);
	return $self;
}

############################################################################
# Assemble HTML 
############################################################################

#***************************************************************************
# Subroutine:  assemble_page
# Description: create HTML page from static files 
#***************************************************************************
sub assemble_page {

	my ($self, $page_ref, $html_ref, $leave_body_open) = @_;

	# DEBUG
	#$devtools->print_hash($page_ref); die;

	# Add doctype statement
	my $doctype_path = $page_ref->{doctype};
	unless ($doctype_path) { die "\n\t DOCTYPE undefined"; } 
	$self->incorporate_ssi($doctype_path, $html_ref);
	
	# open <html>
	$self->open_tag($html_ref, 'html', $indent_size, 0);

	# Assemble head
	$self->assemble_head($page_ref, $html_ref);
	$self->assemble_body($page_ref, $html_ref, $leave_body_open);
	
	# close </html>
	unless ($leave_body_open) {
		$self->close_tag($html_ref, 'html', $indent_size, 0);
	}
}

#***************************************************************************
# Subroutine:  assemble head
# Description: 
#***************************************************************************
sub assemble_head {

	my ($self, $page_ref, $html_ref) = @_;

	# open <head>
	$self->open_tag($html_ref, 'head', $indent_size, 1);
	
	# Add title
	my $title = $page_ref->{title};
	unless ($title) { die "\n\t title undefined"; } 
	push(@$html_ref, "\n    <title>$title</title>\n\n");

	# Add meta: description
	my $description = $page_ref->{description};
	unless ($description) { die "\n\t meta-description undefined"; } 
	push(@$html_ref, "\n    <meta name=\"description\" content=\"$description\">\n\n");
	#$self->incorporate_ssi($meta_description_path, $html_ref);
	
	# Add meta: keywords
	my $keywords = $page_ref->{meta_keywords};
	if ($keywords) { 
		my $meta_keywords = "\n    ";
		$meta_keywords .= "<meta name=\"keywords\" content=\"$keywords\">\n\n\n";
		push(@$html_ref, $meta_keywords);
	}
	
    # Add meta: content type
	my $meta_content_path = $page_ref->{meta_content};
	unless ($meta_content_path) { die "\n\t Meta content undefined"; } 
	unless ($meta_content_path) { die "\n\t meta_content_path undefined"; } 
	$self->incorporate_ssi($meta_content_path, $html_ref);
	
	# Add meta: robots
	my $meta_robots_path = $page_ref->{meta_robots};
	unless ($meta_robots_path) { die "\n\t meta-robots undefined"; } 
	$self->incorporate_ssi($meta_robots_path, $html_ref);
	
	
	# Add meta: author
	my $meta_author_path = $page_ref->{meta_author};
	if ($meta_author_path) { 
		$self->incorporate_ssi($meta_author_path, $html_ref);
	}
	
	# Add meta: owner
	#my $meta_owner_path = $page_ref->{meta_owner};
	#if ($meta_owner_path) { 
	#	$self->incorporate_ssi($meta_owner_path, $html_ref);
	#}
    
    # Add meta: css
	my $meta_css_path = $page_ref->{css};
	unless ($meta_css_path) { die "\n\t link to css undefined"; } 
	$self->incorporate_ssi($meta_css_path, $html_ref);
	
	# Add google analytics 'urchin' tracking script
	my $urchin_path = $page_ref->{urchin};
	unless ($urchin_path) { die "\n\t link to urchin undefined"; } 
	$self->incorporate_ssi($urchin_path, $html_ref);
	
	# close <head>
	$self->close_tag($html_ref, 'head', $indent_size, 1);
}

#***************************************************************************
# Subroutine:  assemble body
# Description: 
#***************************************************************************
sub assemble_body {

	my ($self, $page_ref, $html_ref, $leave_open) = @_;

	# open <body>
	$self->open_tag($html_ref, 'body', $indent_size, 1);
	
	# open CSS division <wrap>
	$self->open_css_div($html_ref, 'wrap', $indent_size, 2, 'id');
	
	# add top navigation
	my $top_navigation_path = $page_ref->{top_nav};
	unless ($top_navigation_path) { die "\n\t Top navigation undefined"; } 
	$self->incorporate_ssi($top_navigation_path, $html_ref);
	
	# add header-image
	my $headimage_path = $page_ref->{head_image};
	unless ($headimage_path) { die "\n\t Header image undefined"; } 
	$self->incorporate_ssi($headimage_path, $html_ref);
	
	# open CSS division <main>
	# TODO: shouldn't this be ID?
	$self->open_css_div($html_ref, 'main', $indent_size, 3, 'class');
	
	# write right section  (sidebar in paleo)
	my $right_panel_path = $page_ref->{right};
	unless ($right_panel_path) { die "\n\t right panel undefined"; } 
	$self->incorporate_ssi($right_panel_path, $html_ref);
	
	# write left section  (main in paleo)
	$self->open_css_div($html_ref, 'left', $indent_size, 4, 'class');
	my $left_panel_path = $page_ref->{left};
	#print "\n\t incorporating main left panel '$left_panel_path'";
	unless ($left_panel_path) { die "\n\t left panel undefined"; } 
	my $done = $self->incorporate_ssi($left_panel_path, $html_ref);
	unless ($done) { 
		$devtools->print_hash($page_ref);
		#die; 
	}	
	# close CSS division <left>
	unless ($leave_open) {
		$self->close_css_div($html_ref, $indent_size, 4);
	}

	# close CSS division <main>
	unless ($leave_open) {
		$self->close_css_div($html_ref, $indent_size, 3);
	}	

	# close CSS division <wrap>
	unless ($leave_open) {
		$self->close_css_div($html_ref, $indent_size, 2);
	}	

	# write footer
	#my $footer_path = $page_ref->{footer};
	#unless ($footer_path) { die "\n\t Footer undefined"; } 
	#$self->incorporate_ssi($footer_path, $html_ref);

	# close </body>
	unless ($leave_open) {
		$self->close_tag($html_ref, 'body', $indent_size, 1);
	}
}

#***************************************************************************
# Subroutine:  incorporate_ssi 
# Description: 
#***************************************************************************
sub incorporate_ssi {

	my ($self, $ssi_path, $html_ref) = @_;
	my @ssi;
	my $file_io = FileIO->new();
	my $read = $file_io->read_input_file($ssi_path, \@ssi);
	my $ssi  = join('', @ssi);
	push(@$html_ref, "$ssi\n\n");
	return $read;
}

#***************************************************************************
# Subroutine:  open_tag
# Description: 
#***************************************************************************
sub open_tag {

	my ($self, $html_ref, $tag_name, $indent_size, $indent_multiple) = @_;
	my $tag = "<$tag_name>";
	my $i_num  = $indent_size * $indent_multiple; 
	my $indent = ' ' x $i_num;
	my $line = $indent . $tag . "\n\n\n";
	push (@$html_ref, $line);
}

#***************************************************************************
# Subroutine:  close_tag
# Description: 
#***************************************************************************
sub close_tag {

	my ($self, $page_ref, $tag_name, $indent_size, $indent_multiple) = @_;
	my $tag = "</$tag_name>";
	my $i_num  = $indent_size * $indent_multiple; 
	my $indent = ' ' x $i_num;
	my $line = $indent . $tag . "\n\n\n";
	push (@$page_ref, $line);
}

#***************************************************************************
# Subroutine:  open_css_div
# Description: 
#***************************************************************************
sub open_css_div {

	my ($self, $html_ref, $div_name, $indent_size, $indent_multiple, $type) = @_;

	my $div    = '<div ' . $type . '="' . $div_name . '">';
	my $i_num  = $indent_size * $indent_multiple; 
	my $indent = ' ' x $i_num;
	my $line   = $indent . $div . "\n\n\n";
	push (@$html_ref, $line);
}

#***************************************************************************
# Subroutine:  close_css_div
# Description: 
#***************************************************************************
sub close_css_div {

	my ($self, $html_ref, $indent_size, $indent_multiple) = @_;
	my $tag = "</div>";
	my $i_num  = $indent_size * $indent_multiple; 
	my $indent = ' ' x $i_num;
	my $line = $indent . $tag . "\n\n\n";
	push (@$html_ref, $line);
}

#***************************************************************************
# Subroutine:  write_page
# Description: create HTML page from static files 
#***************************************************************************
sub write_page {

	my ($self, $page_ref, $html_ref) = @_;

	# Add doctype statement
	my $doctype_path = $page_ref->{outpath};
	unless ($doctype_path) { die "\n\t Outpath undefined"; } 
	my $written = $fileio->write_output_file($doctype_path, $html_ref);
	unless ($written) {
		$devtools->print_hash($page_ref);
		die;
	}

}

############################################################################
# HTML tables
############################################################################

#***************************************************************************
# Subroutine:  convert text table to html
# Description: 
#***************************************************************************
sub convert_text_table_to_html { 

	my ($self, $text_ref, $html_ref, $width) = @_;

	my $indent = "\t\t";
	my %settings;
	$settings{class}  = 'data';
	$settings{indent} = $indent;
	$settings{width}  = $width;

	my $open_table = $self->open_table(\%settings);
	push (@$html_ref, "$open_table\n");	

	foreach my $row (@$text_ref) {
		my @row = split("\t", $row);
		my $row = $self->create_row(\@row);
		push (@$html_ref, "$row\n");	
	}

	my $close_table = $self->close_table(\%settings);
	push (@$html_ref, "$close_table\n");	
}

#***************************************************************************
# Subroutine:  convert text table to html
# Description: 
#***************************************************************************
sub convert_text_table_to_html2 { 

	my ($self, $text_ref, $html_ref) = @_;

	my $indent = "\t\t";
	my %settings;
	$settings{class}  = 'data';
	$settings{indent} = $indent;
	my $open_table = $self->open_table(\%settings);
	push (@$html_ref, "$open_table\n");	

	my $header     = @$text_ref[0];
	my @header_row = split("\t", $header);

	my %columns_with_values;
	my $i = 0;
	foreach my $row (@$text_ref) {
		my @row = split("\t", $row);
		# Check which fields have values
		if ($i >= 1) { 
			my $j =0;
			foreach my $item (@row) {
				if ($item ne "") {
					my $field = $header_row[$j];
					print "<BR> Got value '$item' for '$field'";
				}	
				$j++;
			}
		}

		my $row = $self->create_row(\@row);
		push (@$html_ref, "$row\n");	
		$i++;
	}

	my $close_table = $self->close_table(\%settings);
	push (@$html_ref, "$close_table\n");	
}

#***************************************************************************
# Subroutine:  open table 
# Description: 
#***************************************************************************
sub open_table { 

	my ($self, $table_set_ref) = @_;

	my $indent = $table_set_ref->{indent};	
	my $class  = $table_set_ref->{class};
	my $width  = $table_set_ref->{width};	

	my $table_open = "\n$indent <table";
	if ($class) {
		$table_open .= "  class=$class";
	}
	if ($width) {
		$table_open .= "  width=$width";
	}
	$table_open .= ">\n";
	return $table_open;
}

#***************************************************************************
# Subroutine:  close table
# Description: 
#***************************************************************************
sub close_table { 

	my ($self, $table_set_ref) = @_;
	my $indent = $table_set_ref->{indent};	
	my $table_close = "\n$indent </table><br>\n";	
	return $table_close;
}

#***************************************************************************
# Subroutine:  create row
# Description: 
#***************************************************************************
sub create_row { 

	my ($self, $item_ref, $href) = @_;

	my @row;
	my $done_first_item = undef;
	push (@row, "<tr>");
	
	my $i = 0;
	foreach my $item (@$item_ref) {
		
		$i++;
		
		# insert non breaking space if value empty
		unless ($item) { $item = '&nbsp;'; }
		
		if  ($i eq 1)         {	push (@row, '<td>');	}
		elsif ( $item =~ ',') {	push (@row, '<td>');	}
		else {	push (@row, '<td class=center>');       }
		
		unless ($done_first_item) {
			if ($href) { 
				my $link  = "<a href='$href'> ";
				push (@row, $link); 
			}
			$done_first_item = 'true';
		}
		push (@row, $item);
		push (@row, '</td>');
	}
	push (@row, "</tr>");

	my $row = join('', @row);
	return $row;
}

#***************************************************************************
# Subroutine:  create_row_centred
# Description: 
#***************************************************************************
sub create_row_centred { 

	my ($self, $item_ref, $href) = @_;

	my @row;
	my $done_first_item = undef;
	push (@row, "\n<tr>");
	
	my $i = 0;
	foreach my $item (@$item_ref) {
		
		$i++;
		
		# insert non breaking space if value only whitespace"
		my $item_copy = $item;
		$item_copy =~ s/\s//g;
		if ($item_copy eq '') {
			$item = '&nbsp;';
		}
		
		if  ($i eq 1)         {	push (@row, '<td>');	}
		else {	push (@row, '<td class=center>');       }
		
		unless ($done_first_item) {
			if ($href) { 
				my $link  = "<a href='$href'> ";
				push (@row, $link); 
			}
			$done_first_item = 'true';
		}
		push (@row, $item);
		push (@row, '</td>');
	}
	push (@row, "</tr>");

	my $row = join('', @row);
	return $row;
}

#***************************************************************************
# Subroutine:  create row
# Description: 
#***************************************************************************
sub create_aligned_row { 

	my ($self, $item_ref, $align_prefs, $href) = @_;

	my @row;
	my $done_first_item = undef;
	push (@row, "\n<tr>");
	
	my $i = 0;
	foreach my $item (@$item_ref) {
		
		$i++;
		
		# insert non breaking space if value only whitespace"
		my $item_copy = $item;
		$item_copy =~ s/\s//g;
		if ($item_copy eq '') {
			$item = '&nbsp;';
		}
		
		my $align = $align_prefs->{$i};
		push (@row, "<td class=$align>");
		
		unless ($done_first_item) {
			if ($href) { 
				my $link  = "<a href='$href'> ";
				push (@row, $link); 
			}
			$done_first_item = 'true';
		}
		push (@row, $item);
		push (@row, '</td>');
	}
	push (@row, "</tr>");

	my $row = join('', @row);
	return $row;
}


#***************************************************************************
# Subroutine:  create header row
# Description: 
# Arguments:   
#***************************************************************************
sub create_header_row { 

	my ($self, $item_ref) = @_;

	my @row;
	push (@row, "<tr>");
	foreach my $item (@$item_ref) {
		
		# insert non breaking space if value only whitespace"
		my $item_copy = $item;
		$item_copy =~ s/\s//g;
		if ($item_copy eq '') {
			$item = '&nbsp;';
		}
		push (@row, '<th>');
		push (@row, $item);
		push (@row, '</th>');
	}
	push (@row, "</tr>");

	my $row = join('', @row);
	return $row;
}

############################################################################
# Load files
############################################################################

#***************************************************************************
# Subroutine:  load_file 
# Description: 
#***************************************************************************
sub load_file {

	my ($self, $filename, $file) = @_;
	
	# GET REFERENCE SEQUENCE DATA AND FORMAT DB
	my ($bytesread, $buffer);
	my $num_bytes = 1024;
	my $totalbytes;
	my $untainted_filename;
	open (OUTFILE, ">", "$file") or die "Couldn't open $file for writing: $!";
	while ($bytesread = read($filename, $buffer, $num_bytes)) {
		$totalbytes += $bytesread;
		print OUTFILE $buffer;
	}
	die "Read failure" unless defined($bytesread);
	#unless (defined($totalbytes)) {
		#print "<p>Error: Could not read file ${untainted_filename}, ";
		#print "or the file was zero length.";
	#} else {
	#	print "<p>Done. File $filename uploaded to $file ($totalbytes bytes)";
	#}
	close OUTFILE or die "Couldn't close $file: $!";
	
	# Convert mac line breaks
	my $command = "perl -pi -e 's/\r/\n/g' $file";
	system $command;
	
}

############################################################################
# Tidy HTML 
############################################################################

#***************************************************************************
# Subroutine:  
# Description: 
#***************************************************************************
sub format_html {

	my ($self, $file) = @_;
	
	# Read in file 

	# Create params
	my $indent = 2;
	my $line_wrap;
	my @indent_markers = qw [ p div ]
}

#***************************************************************************
# Subroutine:  
# Description: 
#***************************************************************************
sub create_page_top {

	my ($self, $refseq_data, $site_ref, $output_ref) = @_;

	# Pages reference create refseq HTML
	my $main_path    = $site_ref->{main_path};
	my $site_path    = $site_ref->{site_path};
	my $header_path  = $site_ref->{header_path};
	my $footer_path  = $site_ref->{footer_path};
	my $sidebar_path = $site_ref->{sidebar_path};
	my $image_path   = $site_ref->{image_path};
	my $headers      = $site_ref->{headers};
	#$devtools->print_web_hash($site_ref);

	# write the first part of the header
	my $html_obj = HTML_Utilities->new();
	$html_obj->incorporate_ssi($header_path, $output_ref);

	if ($headers) {
		#$devtools->print_hash($headers);
		my $divide = "\n\n<!-- HEADER IMAGE-->\n\n";
		push (@$output_ref, $divide);
		my $image_name = 'fossil'; 
		my $string = "\n" . '<div id="headimage"><img src="';
		$string .= "$image_path/" . $image_name;
		$string .=  '.png">' . "\n</div>\n";
		push (@$output_ref, $string);
		$divide = "\n\n<!-- HEADER END-->\n\n";
		#print "\n\t string $string\n\n";;
		push (@$output_ref, $divide);
	}
	push (@$output_ref, "\n\n<!-- MAIN-->\n\n");

	# write the sidebar
	if ($sidebar_path) {
		$html_obj->incorporate_ssi($sidebar_path, $output_ref);
	}
	my $html  = "\n<!-- LEFT CONTENT -->\n";
	$html    .= '<div class="left">';
	push (@$output_ref, $html);
}


#***************************************************************************
# Subroutine:   show_output_message 
# Description:  
#***************************************************************************
sub show_output_message {

	my ($self, $message, $output_type) = @_;

	my $f_message;
	unless ($output_type) { 
		$output_type = 'command_line';
	}
	if ($output_type eq 'html') {
		$f_message = "\n\t<BR>$message";
	}
	else {
		$f_message = "\n\t$message";
	}
	print $f_message;
}

#***************************************************************************
# Subroutine:  close_output_section
# Description:  
#***************************************************************************
sub close_output_section {

	my ($self, $output_type) = @_;

	if ($output_type eq 'html') { 
		print "</ul></font>"; 

	}
	else { 
		print "\n\n"; 
	}

}

############################################################################
# EOF
############################################################################

############################################################################
# EOF
############################################################################
