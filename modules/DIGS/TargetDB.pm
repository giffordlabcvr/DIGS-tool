#!/usr/bin/perl -w
############################################################################
# Module:      TargetDB.pm
# Description: Genome/organism sequence data management functions
# History:     December 2009: Created by Robert Gifford 
############################################################################
package TargetDB;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::FileIO;
use Base::Console;

############################################################################
# Globals
############################################################################

my $s_length      = 70;              # Standard line length (for FASTA)
my $line_limit    = 10000000000;     # Maximum number lines in file 
my $blast_program  = 'makeblastdb';
#my $blast_program = 'formatdb';

# Create base objects
my $fileio    = FileIO->new();
my $console   = Console->new();
my $devtools   = DevTools->new();

1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Create a new TargetDB.pm 'object'
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	my %levels;
	$levels{1} = 'grouping';
	$levels{2} = 'organism';
	$levels{3} = 'source_type';
	$levels{4} = 'version';
	
	# Set member variables
	my $self = {
		
		# Paths
		directory_levels     => \%levels,
		line_limit           => $line_limit,
		genome_use_path      => $parameter_ref->{genome_use_path},
		skipindexing_paths   => $parameter_ref->{skipindexing_paths},

	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# Functions
############################################################################

#***************************************************************************
# Subroutine:  format_targets_for_blast
# Description: index unformatted FASTA files for BLAST searching
#***************************************************************************
sub format_targets_for_blast {

	my ($self) = @_;
	
	# Get variables from self
	my $genome_path = $self->{genome_use_path};
	unless ($genome_path) { die "\n\t Path to target files is not set\n\n\n"; }

	# Show warning
	print "\n\t ### WARNING: This function requires standard file extensions for target FASTA files.\n";
	print "\n\t ### i.e. '.fas', '.fa', '.fasta', '.faa', '.fna'\n";
	print "\n\t ### Note: large files may generate split index files (01.nin, 02.nin ... etc)";
	print "\n\t ### When this happens this utility will prompt for formatting, even ";
	print "if index files already exist.\n";

	# Ask user how to handle the process
	my $question = "\n\t Prompt before formatting each target file? ";
	my $prompt = $console->ask_yes_no_question($question);
	if ($prompt eq 'n') { $prompt = undef; }

	# Index genomes by key ( organism | type | version )
	my %server_data;	
	print "\n\n\t Refreshing target data under path '$genome_path'\n";
	$self->read_target_directory(\%server_data);

	my $skip_ref = $self->{skipindexing_paths};
	#$devtools->print_hash($skip_ref); #die;	
	
	# Iterate through and check formatting in each genome
	print "\n\t #~#~# Loading target data\n";
	my $skipped = '0';
	my @keys = sort keys %server_data;
	foreach my $key (@keys) {
		
		# Get genome data
		my $genome_ref  = $server_data{$key};
		my $organism    = $genome_ref->{organism};
		my $type        = $genome_ref->{source_type};
		my $version     = $genome_ref->{version};
		my $path        = $genome_ref->{version_path};
		my $group       = $genome_ref->{grouping};
	
		my @target = ( $group, $organism , $type, $version );
		my $target_id = join ('/', @target);
		if ($skip_ref->{$target_id}) {  
			print "\n\t       Skipping '$organism': '$type' '$version'";
			$skipped++;
			next;
		}
		
		print "\n\t #~#~# Checking: $group;\t'$organism'\t'$type'\t'$version'";
		$self->check_genome_formatting($genome_ref);
		my $unformatted_ref = $genome_ref->{unformatted};
		my $num_unformatted = scalar @$unformatted_ref;
		
		if ($num_unformatted) {  # Do formatting
			print "\n\n\t #~#~# Format target file: $organism, $type, $version";
			foreach my $file (@$unformatted_ref) {
				print "\n\t file '$file'";
			}
			if ($prompt) {
				my $question = "\n\n\t Do you want to format the above files?";
				my $answer   = $console->ask_yes_no_question($question);
				if ($answer eq 'y') {
					$self->format_target_for_blast($genome_ref);
				}	
			}
			else { $self->format_target_for_blast($genome_ref); }
		}
	}
	print "\n\n\t #~#~# Skipped '$skipped' files";
}

#***************************************************************************
# Subroutine:  summarise_targets_long
# Description: summarise data in the 'target genomes' directory
#***************************************************************************
sub summarise_targets_long {

	my ($self) = @_;
	
	# Get variables from self
	my $genome_path = $self->{genome_use_path};
	unless ($genome_path) { die "\n\t Path to genomes is not set\n\n\n"; }
	
	# Index genomes by key ( organism | type | version )
	my %server_data;
	$self->read_target_directory(\%server_data);

	# add header row
	my @summary;
	my @header = ('File', 'Organism', 'Group', 'Data-type', 'Version', 'Scaffolds',
                  '# bases', '# lines');
	my $header = join ("\t", @header);
	$header .= "\n";
	push (@summary, $header);

	# Iterate through, summarising target files
	print "\n\n\t ### Generating detailed summary of target files under '\$DIGS_GENOMES'";
	my @keys = keys %server_data;
	foreach my $key (@keys) {
		
		# Get genome data
		my $genome_ref  = $server_data{$key};
		my $group       = $genome_ref->{grouping};
		my $organism    = $genome_ref->{organism};
		my $type        = $genome_ref->{source_type};
		my $version     = $genome_ref->{version};
		my $path        = $genome_ref->{version_path};
		#$devtools->print_hash($genome_ref); die;
		print "\n\t ### Summarising target files for '$organism'";

		# Iterate through indexing by file stems
		my $files_ref   = $genome_ref->{files};
		my @files;
		foreach my $file (@$files_ref) {
			my @file_bits = split (/\./, $file);
			my $type = pop @file_bits;
			if ($type eq 'fa') {
				push (@files, $file);
			}
		}
		my @sorted = sort @files;
		
		# Get statistics for each file
		foreach my $file (@sorted) {
		
			# Get path to file
			my $chunk_path = $genome_path . "/$path/$file";
			my %data;
			$self->get_target_file_statistics(\%data, $chunk_path);
			my $total_bases    = $data{total_bases};
			my $line_count     = $data{total_lines};
			my $num_scaffolds  = $data{number_scaffolds};
			#print "\n\t FILE $file: $chunk_path";
			
			my @line;
			push (@line, $file);
			push (@line, $organism);
			push (@line, $group);
			push (@line, $type);
			push (@line, $version);
			push (@line, $num_scaffolds);
			push (@line, $total_bases);
			push (@line, $line_count);
			my $line = join("\t", @line);
			push (@summary, "$line\n");
		}
	}

	# Write results to file
	my $summary = "digs-target-dbs-detailed-summary.txt";
	print "\n\n\t ### Writing detailed summary to '$summary'\n";
	$fileio->write_file($summary, \@summary);
}

#***************************************************************************
# Subroutine:  summarise_targets_short
# Description: summarise data in the 'targets' directory
#***************************************************************************
sub summarise_targets_short {

	my ($self) = @_;
	
	#  Get member data structure that describes expected directory structure
	my $genome_path = $self->{genome_use_path};
	unless ($genome_path) { die  "\n\t Path not set\n\n\n"; }
	my $levels_ref  = $self->{directory_levels};
	unless ($levels_ref)  { die  "\n\t Levels not set\n\n\n"; }
	my @levels = sort by_number keys %$levels_ref;

	# Index current, locally-held genome data 
	my @genome_files;
	$fileio->read_directory_tree_leaves($genome_path, \@genome_files, $levels_ref);

	# Iterate through, summarising target files
	print "\n\n\t ### Generating brief summary of target files under '\$DIGS_GENOMES'";
	my %target_keys;
	foreach my $file_ref (@genome_files) {
		
		my $grouping = $file_ref->{grouping};
		my $organism = $file_ref->{organism};
		my $type     = $file_ref->{source_type};
		my $version  = $file_ref->{version};
		my $key = $grouping . '|' .$organism . '|' . $type . '|' . $version;
		$target_keys{$key} = 1;
	}
	
	# Add header row
	my @summary;
	my @header = ('Grouping', 'Organism', 'Data type', 'Version');
	my $header = join("\t", @header);
	push (@summary, "$header\n");

	# Add the summary information rows 
	my @keys = sort keys %target_keys;
	foreach my $key (@keys) {
		$key =~ s/\|/\t/g;
		push (@summary, "$key\n");
	}
	
	# Write results to file
	my $summary = "digs-target-dbs-brief-summary.txt";
	print "\n\n\t ### Writing brief summary to '$summary'\n";
	$fileio->write_file($summary, \@summary);

}

#***************************************************************************
# Subroutine:  get_target_file_statistics
# Description: get metrics on a target genome sequence file 
#***************************************************************************
sub get_target_file_statistics {

	my ($self, $data_ref, $chunk_path) = @_;

	# open the file for reading, or die on failure
	open FILE, "<$chunk_path" 
	or die "\n\tCan't open '$chunk_path'\n";
	
	# Open filehandle
	my $num_lines;
	my $num_scaffolds = 0;
	my $line_count    = 0;
	my $total_bases   = 0;
	
	# Count the number of lines in the file
	while ( <FILE> ) { 
		
		my $line = $_;
		if ($line =~ /^\s*$/)      { next; } # discard blank line
		$line_count++;
		chomp $line;
		if  ($line =~ /^>/)  { $num_scaffolds++; }
		else {
			$line =~ s/\s+//g; # remove whitespace
			my $line_length = length $line;
			$total_bases = $total_bases + $line_length;	
		}
	}
	close FILE;
	
	#print "\n\t # Paths:       $chunk_path: ";
	#print "\n\t # Line count:  $line_count";
	#print "\n\t # Base count:  $total_bases";
	#print "\n\t # Scaffolds:   $num_scaffolds";
	$data_ref->{total_bases}      = $total_bases;
	$data_ref->{total_lines}      = $line_count;
	$data_ref->{number_scaffolds} = $num_scaffolds;
}

############################################################################
# Internals
############################################################################

#***************************************************************************
# Subroutine:  read_target_directory
# Description: read the contents of a 'target DB' directory containing 
#***************************************************************************
sub read_target_directory { 

	my ($self, $data_ref) = @_;
	
	#  Get member data structure that describes expected directory structure
	my $genome_path = $self->{genome_use_path};
	unless ($genome_path) { die  "\n\t Path not set\n\n\n"; }
	my $levels_ref  = $self->{directory_levels};
	unless ($levels_ref)  { die  "\n\t Levels not set\n\n\n"; }
	my @levels = sort by_number keys %$levels_ref;

	# Index current, locally-held target data 
	my @genome_files;
	$fileio->read_directory_tree_leaves($genome_path, \@genome_files, $levels_ref);

	# Iterate through the files
	foreach my $file_ref (@genome_files) {
		
		my $file     = $file_ref->{file};
		my $type     = $file_ref->{source_type};
		my $organism = $file_ref->{organism};
		my $version  = $file_ref->{version};
		
		# Create the path to file
		my $path;
		foreach my $level (@levels) {
			my $field = $levels_ref->{$level};
			my $value = $file_ref->{$field};
			#print "\n\t LEVEL $level\tFIELD $field\tVALUE $value";;	
			$path .= "$value/";
		}
		$file_ref->{version_path}  = $path;
		
		# Store in a hash
		my $key = $organism . '|' . $type . '|' . $version;
		if ($data_ref->{$key}) {
			my $genome_ref = $data_ref->{$key};
			my $files_ref  = $genome_ref->{files};
			push (@$files_ref, $file);
		}
		else {
			my @files;
			push (@files, $file);
			$file_ref->{files} = \@files;
			$data_ref->{$key} = $file_ref;
		}
	}
}

#***************************************************************************
# Subroutine:  check_genome_formatting
# Description: look for pre-existing BLAST indexes in a genome directory 
#***************************************************************************
sub check_genome_formatting {

	my ($self, $genome_ref) = @_;

	# Iterate through indexing by file stems
	my $files_ref   = $genome_ref->{files};
	my %stems;
	my @formatted;
	my @unformatted;
	foreach my $file (@$files_ref) {
		
		my @file_bits = split (/\./, $file);
		my $type = pop @file_bits;
		my $stem = join ('.', @file_bits);
		#print "\n\t NAME $file:\t\t $stem";
		unless ($stem and $type) { 
			die "\n\n\t  filename structure error for '$file'\n\n";
		}
		
		# Parsing the file extension of the target sequence file
		if ($type eq 'fasta' or $type eq 'fas' or $type eq 'fa'
		 or $type eq 'fna' or $type eq 'faa'
		) { 
			$stem .= ".$type";
			$type = 'fa'; # standardise type for FASTA files
		}
	
		# Check if the file looks like a BLAST sequence database file	
		my $result = $self->check_filetype($type);
		unless ($result) {
			print "\n\t   File '$file' is not a recognized BLAST-associated filetype";
		}
		
		# Add this filestem to the list of filestems we have seen the stem
		if ($stems{$stem}) { 
			my $files_ref = $stems{$stem};
			if ($files_ref->{$type}) { 
                print "\n\t Showing hash of file extensions for this file stem\n";
                $devtools->print_hash($files_ref);
                print "\n\t NAME $file, TYPE '$type': $stem";
                print "\n\t This file appears to be a duplicate, exiting....\n\n";
                exit;
			}
			$files_ref->{$type} = 1;
			if ($type eq 'fa') {
				$files_ref->{stem_target} = $file;
			}
		}
		else {
			my %stem_files;
			$stem_files{$type} = 1;
			$stem_files{stemfile_name} = $file;
			if ($type eq 'fa') {
				$stem_files{stem_target} = $file;
			}
			$stems{$stem} = \%stem_files;
		}

	}

	my @stems = keys %stems;
	foreach my $stem (@stems) {
		
		my $types_ref = $stems{$stem};
		my $nin = $types_ref->{nin};
		my $nsq = $types_ref->{nsq};
		my $nhr = $types_ref->{nhr};
		#my $nsi = $chunk_data->{nsi};
		#my $nsd = $chunk_data->{nsd};
		my $target = $types_ref->{stem_target};
		unless ($nhr and $nin and $nsq) {
			
			my $stemfile_data = $stems{$stem};
			my $stemfile_name = $stemfile_data->{stemfile_name};
			if ($target) {
				print "\n\t  adding unformatted target  $target for '$stemfile_name'";
				push (@unformatted, $target); 
			}
			else {
				print "\n\t   No target found for:";
				print "\n\t   filestem '$stem'";
				print "\n\t   from file '$stemfile_name' (ignoring)";
			}
		}
		else {
			push (@formatted, $target); 
		}
	}

	$genome_ref->{formatted}   = \@formatted;
	$genome_ref->{unformatted} = \@unformatted;
}

#***************************************************************************
# Subroutine:  check_filetype
# Description: Check if string $type is a known file extension of BLAST database
#***************************************************************************
sub check_filetype {

	my ($self, $type) = @_;

	my $result = undef;

	if( $type eq 'fa'  or $type eq 'ol' or $type eq 'gz' or

		$type eq 'nnd' or $type eq 'nni' or 
		$type eq 'nsd' or $type eq 'nsi' or 
		$type eq 'nsq' or $type eq 'ntm' or
		$type eq 'nin' or $type eq 'nhr' or
		$type eq 'nin' or $type eq 'nhr' or
		$type eq 'nhd' or $type eq 'nhi' or
		$type eq 'nog' or $type eq 'nal' or
		$type eq 'nog' or $type eq 'nos' or
		$type eq 'nog' or $type eq 'ntf' or
		$type eq 'nog' or $type eq 'nto' or
		$type eq 'nog' or $type eq 'not' or
		$type eq 'fai' or $type eq 'ndb' or
		$type eq 'txt') {

		$result = 1;
	
	}

	return $result;

}

#***************************************************************************
# Subroutine:  format_target_for_blast
# Description: format sequence files in a directory for BLAST 
#***************************************************************************
sub format_target_for_blast {

	my ($self, $genome_ref) = @_;

	# Get data 
	my $genome_path   = $self->{genome_use_path};
	my $line_limit    = $self->{line_limit};
	my $organism      = $genome_ref->{organism};
	my $type          = $genome_ref->{source_type};
	my $version       = $genome_ref->{version};
	my $version_path  = $genome_ref->{version_path};
	my $path          = $genome_path . $version_path;
	my $bin_path = $blast_program;

	# Iterate through the files
	my $formatted_ref   = $genome_ref->{formatted};
	my $unformatted_ref = $genome_ref->{unformatted};
	foreach my $file (@$unformatted_ref) { 
		
		# Get path to file
		my $chunk_path = $path . "/$file";
		my %data;
		#print "\n\t FILE $file: $chunk_path";
		$self->get_target_file_statistics(\%data, $chunk_path);
		my $total_bases    = $data{total_bases};
		my $line_count     = $data{total_lines};
		my $num_scaffolds  = $data{number_scaffolds};
		print "\n\t Target file line count: $line_count";

		if ($num_scaffolds > 50 and $line_count > $line_limit) {
			
			# Work out file splitting params
			my $split_num = int ($line_count / $line_limit);
			if ($split_num < 1) { $split_num++; }
			my @split_chunks;
			my $done_split = undef;
			print "\n\t Splitting $file to approximately $split_num files";
			unless ($split_num eq 1) {
				$done_split = 1;
				print "\n\t Splitting $file to approximately $split_num files";
				$self->split_genome_chunk($path, $file, \@split_chunks);
			}
			else {
				my $file_path = $path . "/$file";
       			my $command = "makeblastdb -in $file_path -dbtype nucl -parse_seqids > /dev/null";
				#my $command = "formatdb -i $file_path -p F -o T > /dev/null";
				print "\n\t$command\n";
				system $command;
			}

			if ($done_split) {
				foreach my $new_chunk_path (@split_chunks) {
       				my $command = "makeblastdb -in $new_chunk_path -dbtype nucl -parse_seqids > /dev/null";
					#my $command = "formatdb -i $new_chunk_path -p F -o T > /dev/null";
					print "\n\t$command\n";
					system $command;
				}
				my $file_path = $path . "/$file";
				my $command = "mv $file_path ./";
				#print "\n\tMoving original file '$file_path' to ./\n";
				#system $command;
			}
		}
		else {
       		my $command = "makeblastdb -in $chunk_path -dbtype nucl -parse_seqids > /dev/null";
			#my $command = "formatdb -i $chunk_path -p F -o T > /dev/null";
			print "\n\t$command\n";
			system $command;
			push (@$formatted_ref, $file);
		}
	}

}

############################################################################
# Splitting large files
############################################################################

#***************************************************************************
# Subroutine:  split_genome_chunk
# Description: split a large sequence file into several files 
#***************************************************************************
sub split_genome_chunk {

	my ($self, $path, $file, $split_chunks_ref) = @_;

	# Get data
	my $line_limit = $self->{line_limit};
	
	# Split the file
	my @scaffold_chunk;
	my $line_count = 0;
	my $chunk_count = 0;
	my $scaffold_count = 0;
	my $approx_divider;
	my $chunk_path = $path . $file;
	open GENOME, "<$chunk_path" or die "\n\tCan't open $chunk_path\n";
	while ( <GENOME> ) {

		my $line = $_;
		if  ($line =~ /^>/)  { 
			$scaffold_count++;
			
			if ($line_count > $line_limit) {
				$chunk_count++;
				my $new_path = $path . $file . '_' . $chunk_count . '.fa';
				$fileio->write_file($new_path, \@scaffold_chunk);			
				push (@$split_chunks_ref, $new_path);
				
				$scaffold_count = 1; # reset scaffold count
				$line_count = 0; # reset line_count
				@scaffold_chunk = (); # reset chunk
				
				# Store the header line in the new chunk
				push (@scaffold_chunk, $line);
			}
			else {
				push (@scaffold_chunk, $line);
			}
		}
		else {
			$line_count++;
			push (@scaffold_chunk, $line);
		}
	}

	# Store last chunk 
	$chunk_count++;
	my $new_path = $path . $file . '_' . $chunk_count . '.fa';
	$fileio->write_file($new_path, \@scaffold_chunk);			
	push (@$split_chunks_ref, $new_path);
	
	# Remove the original file
	my $command = "rm $chunk_path";
	print "\n\t REMOVED file:\n\t $command\n\n";
	system $command;
	
}

#***************************************************************************
# Subroutine:  split_longline_contig
# Description: split single line of sequence into several lines
#***************************************************************************
sub split_longline_contig {

	my ($self, $chunk_path) = @_; 

	# Create path to chunk, and open filehandle
	print "\n\t $chunk_path";
	open TARGET, "<$chunk_path" or die "\n\tCan't open $chunk_path\n";
	my @file_bits = split(/\./, $chunk_path);
	pop @file_bits;
	my $chunk_file_stem = join ('', @file_bits);	
	my $reformat_path = $chunk_file_stem . '_rf.fa';
	print "\n\t $reformat_path";
	open REFORMATTED_CHUNK, ">$reformat_path" or die "\n\tCan't open $reformat_path\n";

	# Run through chunk, extracting matches as we encounter them
	my $i;
	while ( <TARGET> ) {
		
		$i++;
		print "\n\t line $i";
		my $line = $_;
		if    ($line =~ /^\s*$/)   { next; } # discard blank line
		chomp $line;

		if   ($line =~ /^>/)  { 
			# write line
			print REFORMATTED_CHUNK "$line\n"; 
		}
		else {
		
			# Work out coordinates in this line
			$line =~ s/\s+//g; # remove whitespace
			my $line_length = length $line;
			if ($line_length > $s_length) {
				
				my @line = split ('', $line);
				my $nt_counter = 0;
				my $short_line = '';
				foreach my $char (@line) {
					$nt_counter++;
					if ($nt_counter eq $s_length) {
						$nt_counter = 0;
						# write line
						print REFORMATTED_CHUNK "$short_line\n"; 
						$short_line = '';
					}
					else {
						$short_line .= $char;
					}
				}
				if ($short_line) {
					print REFORMATTED_CHUNK "$short_line\n"; 
				}
			}
			else {
				# write line
				print REFORMATTED_CHUNK "$line\n"; 
			}
		}
	}
	close TARGET;
	close REFORMATTED_CHUNK;
}

############################################################################
# Basic
############################################################################

#***************************************************************************
# Subroutine:  by number
# Description: sort an array of integers by ascending numerical order 
#***************************************************************************
sub by_number { $a <=> $b }	

############################################################################
# EOF
############################################################################
