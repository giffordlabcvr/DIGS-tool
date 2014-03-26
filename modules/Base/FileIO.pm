#!/usr/bin/perl -w
############################################################################
# Module:       FileIO.pm 
# Description:  Functions for working with ASCII text files 
# History:      Rob Gifford, Novemeber 2006: Creation
############################################################################
package FileIO;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::DevTools;

############################################################################
# Globals
############################################################################
	
my $devtools = DevTools->new();
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

############################################################################
# 1. Managing files and directories 
############################################################################

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
# Subroutine:  check_file_exists
# Description: check a file exists by attempting to read it
# Arguments:   $file: the name of the file to read
# Returns:     0 if file does exist, $error string otehrwise
#***************************************************************************
sub check_file_exists {

	my ($self, $file) = @_;

	unless (-f $file) {
		if (-d $file) {
			my $error = "file '$file' is directory";
			return $error;
		}
		else {
			my $error = "'$file' cannot be opened";
			return $error;
		}
 	}
	return 0;
}

#***************************************************************************
# Subroutine:  get infile type
# Description: returns the extension of an infile as a 'file_type'
#***************************************************************************
sub get_infile_type {

	my ($self, $file) = @_;
	my @file_bits = split (/\./, $file);
	my $file_type = pop @file_bits;
	return $file_type;
}

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
# Subroutine:  read_directory_tree_leaves
# Description: read a directory tree, and store each 'leaf' 
# Arguments:   $path: path to directory
#              $leaves: array to store 'leaves' (files data) as hashes
# [optional]   $level_codes: hash with correspondence between directory 
#                            levels and values/labels/classifiers etc
#***************************************************************************
sub read_directory_tree_leaves {

	my ($self, $path, $leaves, $level_codes, $file_code) = @_;
	
	unless ($path and $leaves) { die; }
	
	# Set the current state
	my $level = 1;
	my $value;
	if ($level_codes) {
		$value = $level_codes->{$level};
	}
	
	# Read in the top level directory
	my @directory;
	$self->read_directory_to_array($path, \@directory);
	foreach my $file (@directory) {
		
		#print "\n\t test $file";

		my %branch_data;
		$branch_data{$value} = $file;
		
		my $file_path = $path . '/' . $file;
		if (opendir(DIR, $file_path)) {
			$self->recursive_read($file_path, $leaves, \%branch_data, $level, $level_codes, $file_code); 
			#$devtools->print_array($leaves);
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
	
		#print "\n$file: \t #### $level";
		my $file_path = $path . '/' . $file;
		# Recurse if's a directory
		if (opendir(DIR, $file_path)) {
			my $value = $level_codes->{$level};
			$branch_data->{$value} = $file;
			#print "\n\t ### going into \t $file_path";
			#$devtools->print_hash($branch_data);
			$self->recursive_read($file_path, $leaves, $branch_data, $level, $level_codes, $file_code); 
			
			# Reset branch data
			delete $branch_data->{$value};

		}
		elsif ($level eq ($levels + 1)) {
			unless ($file_code) { $file_code = 'file'; }
			# If setting has been given, recode 'file' to some other variable
			# specified in the hash reference $branch_data
			$branch_data->{$file_code} = $file;
			my %data = %$branch_data;
			$data{path} = $file_path;
			push (@$leaves, \%data);
		}
		else { 
			print "\n\t File '$file' is located in internal node";
			print "\n\t Directory depth incorrect - '$level' ne '$levels'";
		}
	}
}

#***************************************************************************
# Subroutine:  read_directory_tree_leaves_simple
# Description: read a directory tree, and store each 'leaf' 
# Arguments:   $path: path to directory
#              $leaves: array to store 'leaves' (files data) as hashes
# [optional]   $level_codes: hash with correspondence between directory 
#                            levels and values/labels/classifiers etc
#***************************************************************************
sub read_directory_tree_leaves_simple {

	my ($self, $path, $leaves_ref) = @_;
	
	unless ($path and $leaves_ref) { die; }
	
	# Read in the top level directory
	my @directory;
	$self->read_directory_to_array($path, \@directory);
	foreach my $file (@directory) {
		
		#print "\n\t test $file";
		my $file_path = $path . '/' . $file;
		if (opendir(DIR, $file_path)) {
			$self->recursive_read2($file_path, $leaves_ref); 
			#$devtools->print_array($leaves);
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
	
		#print "\n$file: \t #### $level";
		my $file_path = $path . '/' . $file;
		# Recurse if's a directory
		if (opendir(DIR, $file_path)) {
			#print "\n\t ### going into \t $file_path";
			#$devtools->print_hash($branch_data);
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


############################################################################
# 2. Basic IO
############################################################################

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
	#print "\n\t File \"$file\" created!\n\n";
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
	#print "\n\t File \"$file\" created!\n\n";
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

############################################################################
# 3. Extracting sections of text using search tokens
############################################################################

#***************************************************************************
# Subroutine:  read_sql_block
# Description: read a NEXUS style block with elements of an SQL statement
#***************************************************************************
sub read_sql_block {

	my ($self, $file_data_ref, $start_mark, $end_mark, $extract_ref) = @_;

	# Extract the block	
	my @block_data;
	$self->extract_text_block($file_data_ref, \@block_data, $start_mark, $end_mark);
	my $block_size = scalar @block_data; 
	unless ($block_size) { # Nothing read
		return 0;
	} 

	# Get the SQL elements from the block
	my $captured = 0;
	foreach my $line (@block_data) {
		chomp $line;
		if    ($line =~ /^\s*$/)   { next; } # discard blank line
		elsif ($line =~ /^\s*#/)   { next; } # discard comment line 
		my @bits = split(':', $line);

		my $field = $bits[0];
		$field =~ s/\s+//;
		my $value = $bits[1];
		$value =~ s/;//;
		if ($field and $value) {
			$captured++;
			$extract_ref->{$field} = $value;
			#print "\n\t $field='$value'"; # DEBUG
		}
	}	
	return $captured;
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
			#print "\n\t $field='$value'"; # DEBUG
		}
	}
	return $captured;
}

#***************************************************************************
# Subroutine:  extract_text_blocks
# Description: extract block(s) of text denoted by start and stop tokens 
#***************************************************************************
sub extract_text_blocks { 

	my ($self, $input_ref, $blocks_ref, $start_token, $stop_token) = @_;

	my $extract = undef;
	my $i;
	my @block;
	foreach my $line (@$input_ref) {
		if ($line =~ $stop_token and $extract) { 
			#$devtools->print_array(\@block);die;
			$extract = undef;
			my @copy = @block;
			push (@$blocks_ref, \@copy);
			undef @block;
		}
		if  ($extract) { 
			push (@block, $line); 
		}
		elsif ($line =~ $start_token) { 
			$extract = 'true';         
		}
	}
	return $extract;
}

#***************************************************************************
# Subroutine:  extract_text_block
# Description: extract a block of text denoted by start and stop tokens 
#***************************************************************************
sub extract_text_block { 

	my ($self, $input_ref, $block_ref, $start_token, $stop_token) = @_;

	my $extract = undef;
	foreach my $line (@$input_ref) {
		if ($stop_token) {
			if ($line =~ $stop_token) { 
				$extract = undef;          
			}
		}	
		if  ($extract) { 
			push (@$block_ref, $line); 
		}
		elsif ($line =~ $start_token) { 
			$extract = 'true';         
		}
	}

}

#***************************************************************************
# Subroutine:  create_unique_directory 
# Description: create a unique directory 
#***************************************************************************
sub create_unique_directory {
	
	my ($self, $unique_dir, $mode) = @_;

	my $mkdir_cmd = "mkdir $unique_dir";
	#print "\n\t$mkdir_cmd\n\n";
	my $result = system $mkdir_cmd;
	if ($result > 0) {
		if ($mode) {
			die "<br>Internal error, contact webmaster<br>";
		}
		else {
			die "\n\t ### Error: couldn't create output directory - check permissions\n\n";	
		}
	}
}

############################################################################
# EOF
############################################################################
