#!/usr/bin/perl -w
############################################################################
# Module:      FileIO.pm 
# Description: Functions for working with ASCII text files 
# History:     Rob Gifford, Novemeber 2006: Creation
############################################################################
package FileIO;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

############################################################################
# Globals
############################################################################
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create a new FileIO object
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Member variables
	my $self = {
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# Public Member Functions
############################################################################

#***************************************************************************
# Subroutine:  get_infile_type
# Description: returns the extension of an infile as a 'file_type'
#***************************************************************************
sub get_infile_type {

	my ($self, $file) = @_;
	my @file_bits = split (/\./, $file);
	my $file_type = pop @file_bits;
	return $file_type;
}

#***************************************************************************
# Subroutine:  check_directory_exists
# Description: routine to check if a directory path is valid
#***************************************************************************
sub check_directory_exists {

	my ($self, $directory) = @_;
	unless (opendir(DIR, $directory)) {
		#die "\n\t Cannot open directory \"$directory\"\n\n";
		return undef;
	}
	return 1;
}

#***************************************************************************
# Subroutine:  read_file
# Description: read an input file to an array
# Arguments:   $file: the name of the file to read
#              $array_ref: array to copy to
#***************************************************************************
sub read_file {

	my ($self, $file, $array_ref) = @_;

	unless (-f $file) {
		if (-d $file) {
			print "\n\t Cannot open file \"$file\" - it is a directory\n\n";
			return 0;
		}
		else {
			print "\n\t Cannot open file \"$file\"\n\n";
			return 0;
		}

 	}
	unless (open(INFILE, "$file")) {
		print "\n\t Cannot open file \"$file\"\n\n";
		return 0;
	}
	@$array_ref = <INFILE>;
	close INFILE;

	return 1;
}

#***************************************************************************
# Subroutine:  write_file
# Description: write an array to an ouput file
# Arguments:   $file: the name of the file to write to 
#              $array_ref: array to copy
#***************************************************************************
sub write_file {

	my ($self, $file, $array_ref) = @_;
	unless (open(OUTFILE, ">$file")) {
		print "\n\t Couldn't open file \"$file\" for writing\n\n";
		return 0;
	}
	print OUTFILE @$array_ref;
	close OUTFILE;
}

#***************************************************************************
# Subroutine:  write_text_to_file
# Description: write a formatted text string to an ouput file
# Arguments:   $file: the name of the file to write to 
#              $text; the string to write out
#***************************************************************************
sub write_text_to_file {

	my ($self, $file, $text) = @_;

	unless (open(OUTFILE, ">$file")) {
		print "\n\t Couldn't open file \"$file\" for writing\n\n";
		return 0;
	}
	print OUTFILE $text;
	close OUTFILE;
}

#***************************************************************************
# Subroutine:  append_text_to_file
# Description: append a text string to an ouput file
# Arguments:   $file: the name of the file to write to 
#              $text; the string to append
#***************************************************************************
sub append_text_to_file {

	my ($self, $file, $text) = @_;

	unless (open(OUTFILE, ">>$file")) {
		print "\n\t Couldn't open file \"$file\" for writing\n\n";
		return;
	}
	print OUTFILE $text;
	close OUTFILE;
	#print "\n\t File \"$file\" created!\n\n";
}

#***************************************************************************
# Subroutine:  read_standard_field_value_block
# Description: read a NEXUS style block containing [field]=[value] lines 
# Returns:     1 if block was found and lines were read
#***************************************************************************
sub read_standard_field_value_block {

	my ($self, $file_data_ref, $start_mark, $end_mark, $extract_ref) = @_;

	# Extract the block	
	my @block_data;
	$self->extract_text_block($file_data_ref, \@block_data, $start_mark, $end_mark);
	my $block_size = scalar @block_data; 
	unless ($block_size) { # Nothing read
		return 0;
	} 

	# Get the field-value pairs from the block
	my $captured = 0;
	foreach my $line (@block_data) {
		
		chomp $line;
		if    ($line =~ /^\s*$/)   { next; } # discard blank line
		elsif ($line =~ /^\s*#/)   { next; } # discard comment line 
		
		my @bits = split('=', $line);
		my $field = $bits[0];
		$field =~ s/\s+//;
		my $value = $bits[1];
		$value =~ s/\s+//;
		$value =~ s/;//;
		if ($field and $value) {
			$captured++;
			$extract_ref->{$field} = $value;
		}
	}
	return $captured;
}

############################################################################
# Reading from directories
############################################################################

#***************************************************************************
# Subroutine:  read_directory_to_array
# Description: read a directory and copy the names of all the files in it 
#              to an array
# Arguments:   $directory: the path to the directory we're reading
#              $array_ref: array to copy to
#***************************************************************************
sub read_directory_to_array {

	my ($self, $directory, $array_ref) = @_;
		
	unless (opendir(DIR, $directory)) {
		print "\n\t Cannot open directory \"$directory\"\n\n";
		return;
	}
	my $file;
	while( defined ($file = readdir(DIR))) {
		# don't copy anything that starts with a '.'
		unless ($file =~ /^\./) {
			push(@$array_ref, $file);
		}
	}
}

#***************************************************************************
# Subroutine:  read_directory_tree_leaves_simple
# Description: read a directory tree, and store each 'leaf' 
# Arguments:   $path: path to directory
#              $leaves: array to store 'leaves' (files) as hashes
#***************************************************************************
sub read_directory_tree_leaves_simple {

	my ($self, $path, $leaves_ref) = @_;
	
	unless ($path and $leaves_ref) { die; }
	
	# Read in the top level directory
	my @directory;
	$self->read_directory_to_array($path, \@directory);
	foreach my $file (@directory) {
		
		# Create the file path	
		my $file_path = $path . '/' . $file;
		if (opendir(DIR, $file_path)) {
			$self->recursive_read2($file_path, $leaves_ref); 
		}
		else { # If its not a directory store it
			my %file;
			$file{path} = $file_path;
			$file{file} = $file;
			push (@$leaves_ref, \%file);
		}
	}
}

#***************************************************************************
# Subroutine:  read_directory_tree_leaves
# Description: read a directory tree, and store each 'leaf' 
# Arguments:   $path: path to directory
#              $leaves: array to store 'leaves' (files data) as hashes
# [optional]   $level_codes: hash with correspondence between directory levels and values/labels/classifiers etc
#***************************************************************************
sub read_directory_tree_leaves {

	my ($self, $path, $leaves, $level_codes, $file_code) = @_;
	
	unless ($path and $leaves) { die; }
	
	# Set thE current state
	my $level = 1;
	my $value;
	if ($level_codes) {
		$value = $level_codes->{$level};
	}
	
	# Read in the top level directory
	my @directory;
	$self->read_directory_to_array($path, \@directory);
	foreach my $file (@directory) {
		my %branch_data;
		$branch_data{$value} = $file;
		my $file_path = $path . '/' . $file;
		if (opendir(DIR, $file_path)) {
			$self->recursive_read($file_path, $leaves, \%branch_data, $level, $level_codes, $file_code); 
		}
	}
}

############################################################################
# FASTA IO
############################################################################

#***************************************************************************
# Subroutine:  read_fasta
# Description: read a fasta file into an array of hashes. 
#***************************************************************************
sub read_fasta {

	my ($self, $file, $array_ref, $length_only) = @_;
	
	# Read in the file or else return
	unless (open(INFILE, $file)) {
		print "\n\t Cannot open file \"$file\"\n\n";
		return undef;
	}

	# Iterate through lines in the file
	my @raw_fasta = <INFILE>;
	close INFILE;
	my $header;
    my $sequence;
	foreach my $line (@raw_fasta) {
		
		chomp $line;
		if    ($line =~ /^\s*$/)   { next; } # discard blank line
		elsif ($line =~ /^\s*#/)   { next; } # discard comment line 
		elsif ($line =~ /^>/) {
			
			$line =~ s/^>//g;
					
			# new header, store any sequence held in the buffer
			if ($header and $sequence) {
				$sequence = uc $sequence;
				my $length = length $sequence;
				my %seq;
				$seq{header} = $header;
				unless ($length_only) {
					$seq{sequence} = $sequence;			
				}
				$seq{seq_length} = $length;				
				push(@$array_ref, \%seq);
			}
		
			# reset the variables 
			$line =~ s/^>//;
			$header = $line;
			$sequence = undef;
		}
		else {
			# keep line, add to sequence string
            $sequence .= $line;
     	}
    }
	
	# Before exit, store any sequence held in the buffer
	if ($header and $sequence) {
		$sequence =~ s/\s+//g; # Remove whitespace
		$sequence = uc $sequence;
		my $length = length $sequence;
		my %seq;
		$seq{header} = $header;
		unless ($length_only) {
			$seq{sequence} = $sequence;			
		}
		$seq{seq_length} = $length;				
		push(@$array_ref, \%seq);
	
	}
}

#***************************************************************************
# Subroutine:  recursive_read2
# Description: read everything under a given path
#***************************************************************************
sub recursive_read2 {

	my ($self, $path, $leaves_ref) = @_;
	
	# Read in the top level directory
	my @directory;
	$self->read_directory_to_array($path, \@directory);
	
	# Iterate through
	foreach my $file (@directory) {
	
		my $file_path = $path . '/' . $file;

		# Recurse if's a directory
		if (opendir(DIR, $file_path)) {
			$self->recursive_read2($file_path, $leaves_ref); 
		}
		else {
			my %file;
			$file{path} = $file_path;
			$file{file} = $file;
			push (@$leaves_ref, \%file);
		}
	}
}

#***************************************************************************
# Subroutine:  recursive_read
# Description: read everything under a given directory path 
#***************************************************************************
sub recursive_read {

	my ($self, $path, $leaves, $branch_data, $level, $level_codes, $file_code) = @_;
	
	# Increment level 
	$level++;

	# Read in the top level directory
	my @directory;
	$self->read_directory_to_array($path, \@directory);
	
	# Get the total number of levels
	my @levels = keys %$level_codes;
	my $levels = scalar @levels;

	# Iterate through
	foreach my $file (@directory) {
		my $file_path = $path . '/' . $file;
		# Recurse if's a directory
		if (opendir(DIR, $file_path)) {
			my $value = $level_codes->{$level};
			$branch_data->{$value} = $file;
			$self->recursive_read($file_path, $leaves, $branch_data, $level, $level_codes, $file_code); 
			# Reset branch data
			delete $branch_data->{$value};

		}
		elsif ($level eq ($levels + 1)) {
			unless ($file_code) { $file_code = 'file'; }
			$branch_data->{$file_code} = $file;
			my %data = %$branch_data;
			$data{path} = $file_path;
			push (@$leaves, \%data);
		}
		else { 
			print "\n\t File '$file' is located in internal node";
		}
	}
}

#***************************************************************************
# Subroutine:  extract_text_block
# Description: extract a block of text denoted by start and stop tokens 
#***************************************************************************
sub extract_text_block { 

	my ($self, $input_ref, $block_ref, $start_token, $stop_token) = @_;

	# Make start and stop tokens upper case so tokens are not case sensitive
	$start_token = uc $start_token;
	$stop_token  = uc $stop_token;
	#print "\n\t ###  Stop token:  $stop_token";

	my $extract = undef;
	foreach my $line (@$input_ref) {
		#print "\n\t ###  Start token: $start_token";
		chomp $line;
		my $uc_line = uc $line;
		#print "\n\t ###  LINE $uc_line;";
		if ($stop_token) {
			if ($uc_line =~ $stop_token) {
				$extract = undef;          
			}
		}	
		if  ($extract) { 
			push (@$block_ref, $line); 
		}
		elsif ($uc_line =~ $start_token) { 
			$extract = 'true';         
		}
	}
}

#***************************************************************************
# Subroutine:  create_unique_directory 
# Description: create a unique directory 
#***************************************************************************
sub create_directory {
	
	my ($self, $unique_dir) = @_;

	my $mkdir_cmd = "mkdir $unique_dir";
	my $result = system $mkdir_cmd;
	if ($result > 0) {
		print "\n\t ### Error: couldn't create output directory using command 'mkdir $unique_dir'\n\n";
		die;
	}
}


############################################################################
# EOF
############################################################################