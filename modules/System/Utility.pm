#!/usr/bin/perl -w
############################################################################
# Module:      Utility.pm
# Description: The Utility module for Pipeline databases
#              Contains routines for creating, and interacting with 
#              Utilitys in the Pipeline framework
# History:     June 2011: Created by Robert Gifford 
############################################################################
package Utility;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::SeqIO;
use Base::FileIO;
use Base::DevTools;
use Base::Console;

############################################################################
# Globals
############################################################################

# Create base objects
my $seqio     = SeqIO->new();
my $fileio    = FileIO->new();
my $devtools  = DevTools->new();
my $console   = Console->new();
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
		
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# Public member functions
############################################################################

#***************************************************************************
# Subroutine:  run_utility_function 
# Description: top level handler 
#***************************************************************************
sub run_utility_function {

	my ($self, $glue_mode) = @_;
	
	# Set the usage statement
   	my $USAGE .= "\n\t  :\n";
   	$USAGE    .= "\n\t     1 = consolidate library";
   	$USAGE    .= "\n\t     2 = clone site";
   	$USAGE    .= "\n\t     3 = merge libraries"; 
   	$USAGE    .= "\n\t     4 = sort library for screening";
   	$USAGE    .= "\n\t     5 = astrid sort";
	$USAGE    .= "\n\n";

	# Load data structure describing basic website
	# TO DO: Move from here, into fxns, dpoesnt always need 'load site'
	$self->show_title();
	
	if ($glue_mode eq 1) { 
		$self->consolidate_library();
	}
	# Refresh entire website except screening DB pages
	elsif ($glue_mode eq 2) { 
		$self->clone_site();
	}
	# Temp utility function
	elsif ($glue_mode eq 3) { 
		my $path1 = './db/refseq_flat/';
		my $path2 = './db/refseq_flat2/';
		$self->merge_libraries($path1, $path2);
	}
	elsif ($glue_mode eq 4) {
		$self->sort_for_screening();
	}
	elsif ($glue_mode eq 5) {
		my $in_dir  = $self->{input_dir};
		my @infile  = split("\/", $in_dir);
		my $dir = pop @infile; # print $dir; die;
		my $out_dir = "./test/ASTRID/converted/$dir";
		my $command = "mkdir $out_dir";
		system $command;
		$self->run_astrid_process($out_dir);
	}
	else { die $USAGE; }	
}

#***************************************************************************
# Subroutine:  create_report_directory 
# Description: create a unique directory for process and results 
#***************************************************************************
sub create_report_directory {
	
	my ($self, $report_dir, $mode) = @_;

	my $mkdir_cmd = "mkdir $report_dir";
	#print "\n\t$mkdir_cmd\n\n";
	my $result = system $mkdir_cmd;
	if ($result > 0) {
		if ($mode) {
			die "\n\t ### Error: couldn't create output directory - check permissions\n\n";	
		}
		else {
			die "<br>Internal error, contact webmaster<br>";
		}
	}
}

############################################################################
# PARSE CONTROL FILE
############################################################################

#***************************************************************************
# Subroutine:  parse control file
# Description: read an input file to get parameters for screening
#***************************************************************************
sub parse_control_file {

	my ($self, $ctl_file_path, $hash_ref) = @_;
	
	# see above	
	my @file;
	$fileio->read_input_file($ctl_file_path, \@file);
	
	#my %data;
	my $start = 'BEGIN PARAMS';
	my $stop  = 'ENDBLOCK';
	$fileio->read_standard_field_value_block(\@file, $start, $stop, $hash_ref);
	#$devtools->print_hash($pipeline_obj);

	# Get genome paths
	my @target_block;
	$start = 'BEGIN TARGETS';
	$stop  = 'ENDBLOCK';
	$fileio->extract_text_block(\@file, \@target_block, $start, $stop);
	my @targets;
	foreach my $line (@target_block) {
		
		if ($line =~ /^\s*#/)   { next; } # discard comment line 
		chomp $line;
		push (@targets, $line);
	}
	$hash_ref->{targets} = \@targets;	

	# Get refseq list
	my @refseq_block;
	$start = 'BEGIN REFSEQLIST';
	$stop  = 'ENDBLOCK';
	$fileio->extract_text_block(\@file, \@refseq_block, $start, $stop);
	my @refseqs;
	foreach my $line (@refseq_block) {
		
		if ($line =~ /^\s*#/)   { next; } # discard comment line 
		chomp $line;
		push (@refseqs, $line);
	}
	$hash_ref->{refseqlist} = \@refseqs;	
}

############################################################################
# ASTRID PROCESS: Set up the data
############################################################################

#***************************************************************************
# Subroutine:  run_astrid_process
# Description: run asrid batch process 
#***************************************************************************
sub run_astrid_process {

	my ($self, $report_dir) = @_;
	
	my $in_dir = $self->{input_dir};
	my @dir;
	$fileio->read_directory_to_array($in_dir, \@dir);
	#$devtools->print_array(\@dir); die;	
    
	# Write raw data to report directory 
    foreach my $subdir (@dir) {

		my @subdir;
		my $subdir_path = $in_dir . "/$subdir";
		$fileio->read_directory_to_array($subdir_path, \@subdir);
		#$devtools->print_array(\@subdir); die;	

		# Make the directories for sorting
		my $subdirec  = $report_dir . "/$subdir/";
		my $command = " mkdir $subdirec";
		system $command;

		my $tmp_dir = $report_dir . "/$subdir/sort_1";
		my $command1     = " mkdir $tmp_dir";
		system $command1;
		my $batch_dir = $report_dir . "/$subdir/sort_2";
		my $command2     = " mkdir $batch_dir";
		system $command2;

		#$infile =~ s/\s+//g;
		# Run first step 
		my $input_subdirec = $in_dir . "/$subdir";
		$self->do_first_step($input_subdirec, $tmp_dir);
		# Run second step 
		$self->do_second_step($tmp_dir, $batch_dir);
		
		# Clean up first directory		
		my $command3     = "rm $tmp_dir/*";
		system $command3;
		
		my $command4     = "mv $batch_dir/* $report_dir/$subdir/";
		system $command4;
	
		my $command5     = "rmdir $batch_dir $tmp_dir";
		system $command5;
	
	}
}

#***************************************************************************
# Subroutine:  do_first_step
# Description: do first step of preprocessing 
#***************************************************************************
sub do_first_step {

	my ($self, $directory, $output_path) = @_;
	
	# process the directory
    print "\n\t ## Preparing data in directory '$directory'\n";
	my @file_bits = split (/\//, $directory);
	my $toplevel = pop @file_bits;
	my @dir_split = split(/_/, $toplevel);
	my $gene = shift @dir_split;
	#print "\n\t GENE$gene"; die;

	my %directory;
	$fileio->read_directory_to_hash($directory, \%directory);
	my @keys = keys %directory;
	# $devtools->print_array(\@keys); die;

	foreach my $file (@keys) {
		
		# Create data structures for output
		my @data_out;
		my @seqs_out;
		my @file_bits = split (/\./, $file);
		my $stem = shift @file_bits;
		my @split = split('-', $stem);
		my $mod = pop @split;	
		my @mod = split('', $mod);
		my $char = pop @mod;
		print "\n\t Processing FILE $file -i GENE $gene KEY $stem:  ($char)";

		# Read in the sequences
		my @seqfile;
		my $seqpath = $directory . '/' . $file;
		$seqio->read_fasta($seqpath, \@seqfile);
		
		foreach my $seq_ref (@seqfile) {
			my $header = $seq_ref->{header};
			my $sequence_id  = $seq_ref->{sequence_id};
			my @header  = split(/\s+/, $header); 
			my $subtype = $header[0];
			my $patnum  = $header[2];
			my $country = $header[1];
			my $year    = $header[3];
			my $hla_A1  = $header[4];
			my $hla_A2  = $header[5];
			my $hla_B1  = $header[6];
			my $hla_B2  = $header[7];
			unless ($hla_A1) { $hla_A1 = 'NULL'; }
			unless ($hla_A2) { $hla_A2 = 'NULL'; }
			unless ($hla_B1) { $hla_B1 = 'NULL'; }
			unless ($hla_B2) { $hla_B2 = 'NULL'; }
			$header =~ s/^\s+//;  #remove leading spaces
			$header =~ s/\s+$//;  #remove trailing spaces			
			$header =~ s/\s+//g; # fill internal whitespace
			#$devtools->print_array(\@header);
			my $sequence  = $seq_ref->{sequence};
			$sequence =~ s/-//g;
			my @data_line;
			#push (@data_line, $header);
			push (@data_line, $sequence_id);
			#push (@data_line, $subtype);
			#push (@data_line, $patnum);
			#push (@data_line, $country);
			#push (@data_line, $year);
			#push (@data_line, $hla_A1);
			#push (@data_line, $hla_A2);
			#push (@data_line, $hla_B1);
			#push (@data_line, $hla_B2);
			#push (@data_line, $header);
			push (@data_line, $mod);
			#$devtools->print_array(\@data_line); die;
			my $data_line = join("\t", @data_line);
			$data_line =~ s/#//g;
			my $seq_line  = ">$sequence_id\n$sequence\n\n";
			#my $seq_line  = ">$header\n$sequence\n\n";
			push (@seqs_out, $seq_line); 
			push (@data_out, "$data_line\n"); 
		}
		# Write files out
		$file =~ s/-//g;
		my $seqs_out = $output_path . '/' . $gene . '_' . $file . '.fas';
		$fileio->write_output_file($seqs_out, \@seqs_out);
		$file =~ s/\.fasta//g;
		my $data_out = $output_path . '/' . $gene . '_' . $file . '.data.txt';
		$fileio->write_output_file($data_out, \@data_out);
	}
}

#***************************************************************************
# Subroutine:  do_second_step 
# Description: do second step of preprocessing 
#***************************************************************************
sub do_second_step {

	my ($self, $in_directory, $out_directory) = @_;

	# Read in the directory
	my %directory;
	$fileio->read_directory_to_hash($in_directory, \%directory);
	#my @keys = keys %directory;
	#$devtools->print_hash(\%directory); die;	
	my @keys = keys %directory;

	my %merged;
	foreach my $file (@keys) {
	
		#print "\n\t  $file";
		# Create data structures for output
		my @data_out;
		my @seqs_out;
		my @file_bits = split (/\./, $file);
		my $stem = shift @file_bits;
		my $stem2 = $stem;
		my @split = split('-', $stem);
		my $mod = pop @split;	
		my @mod = split('', $mod);
		my $char = pop @mod;
		my @split2 = split('',$stem2);
		pop @split2;	
		my $short_stem = join ('',@split2);
		my $suffix = pop @file_bits;
		my $key = $short_stem . '_' . $suffix;
		#print "\n\t STEM $stem2 SHORT STEM $short_stem KEY $key:  ($char)";

		if ($merged{$key}) {
			
			my @file1;
			my $file1 = $in_directory . '/' . $file;
			$fileio->read_input_file($file1, \@file1);
			
			my $file2 = $merged{$key};
			my $path  = $in_directory . '/' . $file2;
			my @file2;
			$fileio->read_input_file($path, \@file2);
		
			my @combined = (@file1, @file2);
        	if ($suffix eq 'txt') {
				unshift (@combined, "sequence_id\tHLA\n"); 
			}
			
			my $file_name = $short_stem . 'ab.' . $suffix;
			my $set_dir = $out_directory . '/' . $short_stem . '_ab';
			unless ($fileio->check_directory_exists($set_dir)) {	
				system "mkdir $set_dir";
			}
			my $combined_path = $set_dir . '/' . $short_stem . 'ab.' . $suffix;
			#my $combined_path = $out_directory . '/' . $short_stem . 'ab.' . $suffix;
		 	#print "\n\t # key '$key': Writing $file and $file2 to $combined_path";
			$fileio->write_output_file($combined_path, \@combined);
		
		}
		else {
			$merged{$key} = $file;	
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
	my $description = 'GLUE utilities';
	my $author      = 'Robert J. Gifford';
	my $contact		= '<rgifford@adarc.org>';
	$console->show_about_box($title, $version, $description, $author, $contact);
}

############################################################################
# EOF
############################################################################
