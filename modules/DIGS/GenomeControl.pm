#!/usr/bin/perl -w
############################################################################
# Module:      GenomeControl.pm
# Description: Genome/organism sequence data management functions
# History:     December 2009: Created by Robert Gifford 
############################################################################
package GenomeControl;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::FileIO;
use Base::SeqIO;
use Base::DevTools;
use Base::Console;
use Base::Sequence;

############################################################################
# Globals
############################################################################

my $s_length    = 70;           # Standard line length
my $line_limit  = 10000000;     # Maximum number lines in file 

# Create base objects
my $fileio    = FileIO->new();
my $seqio     = SeqIO->new();
my $devtools  = DevTools->new();
my $console   = Console->new();
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
		blast_bin_path       => $parameter_ref->{blast_bin_path},
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# Functions
############################################################################

#***************************************************************************
# Subroutine:  summarise_genomes
# Description: summarise data in the 'target genomes' directory
#***************************************************************************
sub summarise_genomes {

	my ($self) = @_;
	
	# Get variables from self
	my $genome_path = $self->{genome_use_path};
	unless ($genome_path) { die "\n\t Path to genomes is not set\n\n\n"; }
	
	# Index genomes by key ( organism | type | version )
	my %server_data;
	$self->read_genome_directory(\%server_data);

	# add header row
	my @summary;
	my @header = ('File', 'Organism', 'Data type', 'Version', 'Scaffolds',
                  '# nucleotides', '# lines');
	my $header = join ("\t", @header);
	$header .= "\n";
	push (@summary, $header);

	# Iterate through and check formatting in each genome
	print "\n\n\t ### Summarising genomes\n";
	my @keys = keys %server_data;
	foreach my $key (@keys) {
		
		# Get genome data
		my $genome_ref  = $server_data{$key};
		my $organism    = $genome_ref->{organism};
		my $type        = $genome_ref->{source_type};
		my $version     = $genome_ref->{version};
		my $path        = $genome_ref->{version_path};
		print "\n\t ### Summarising target files for '$organism'";
		#$devtools->print_hash($genome_ref);

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
		#$devtools->print_array(\@sorted);
		
		# Get statistics for each file
		foreach my $file (@sorted) {
		
			# Get path to file
			my $chunk_path = $genome_path . "/$path/$file";
			my %data;
			$self->get_genome_chunk_stats(\%data, $chunk_path);
			my $total_bases    = $data{total_bases};
			my $line_count     = $data{total_lines};
			my $num_scaffolds  = $data{number_scaffolds};
			#print "\n\t FILE $file: $chunk_path";
			#$devtools->print_hash(\%data); die;
			
			my @line;
			push (@line, $file);
			push (@line, $organism);
			push (@line, $type);
			push (@line, $version);
			push (@line, $num_scaffolds);
			push (@line, $total_bases);
			push (@line, $line_count);
			my $line = join("\t", @line);
			push (@summary, "$line\n");
		}
	}

	# write results to file
	my $summary = "target_genomes_summary.txt";
	print "\n\n\t ### Writing summary to '$summary'\n";
	$fileio->write_file($summary, \@summary);
}

#***************************************************************************
# Subroutine:  refresh_genomes
# Description: refresh the 'target genomes' directory
#***************************************************************************
sub refresh_genomes {

	my ($self) = @_;
	
	# Get variables from self
	my $genome_path = $self->{genome_use_path};
	unless ($genome_path) { die "\n\t Path to genomes is not set\n\n\n"; }
	
	# Index genomes by key ( organism | type | version )
	my %server_data;
	$self->read_genome_directory(\%server_data);
	#$devtools->print_hash(\%server_data); die;

	# Iterate through and check formatting in each genome
	my @keys = keys %server_data;
	foreach my $key (@keys) {
		
		# Get genome data
		my $genome_ref  = $server_data{$key};
		my $organism    = $genome_ref->{organism};
		my $type        = $genome_ref->{source_type};
		my $version     = $genome_ref->{version};
		my $path        = $genome_ref->{version_path};
		#$devtools->print_hash($genome_ref);
		
		print "\n\t #~#~# Checking '$organism' genome: '$type' version '$version'";
		$self->check_genome_formatting($genome_ref);
		my $unformatted_ref = $genome_ref->{unformatted};
		my $num_unformatted = scalar @$unformatted_ref;
		
		if ($num_unformatted) {  # Do formatting
			print "\n\n\t #~#~# Format genome: $organism, $type, $version";
			my $question = "\n\t Do you want to format unformatted chunks in this genome?";
			my $answer   = $console->ask_yes_no_question($question);
			if ($answer eq 'y') {
				$self->format_genome($genome_ref);
			}	
		}
	}
}

#***************************************************************************
# Subroutine:  read_genome_directory
# Description: read in the contents of a structured directory containing
#              genome data files
#***************************************************************************
sub read_genome_directory { 

	my ($self, $data_ref) = @_;
	
	#  Get member data structure that describes expected directory structure
	my $genome_path = $self->{genome_use_path};
	unless ($genome_path) { die  "\n\t Path not set\n\n\n"; }
	my $levels_ref  = $self->{directory_levels};
	unless ($levels_ref)  { die  "\n\t Levels not set\n\n\n"; }
	my @levels = sort by_number keys %$levels_ref;

	# Index current, locally-held genome data 
	my @genome_files;
	$fileio->read_directory_tree_leaves($genome_path, \@genome_files, $levels_ref);

	# Iterate through the files
	foreach my $file_ref (@genome_files) {
		
		#$devtools->print_hash($file_ref);
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
# Description: 
#***************************************************************************
sub check_genome_formatting {

	my ($self, $genome_ref) = @_;

	# Iterate through indexing by file stems
	my $files_ref   = $genome_ref->{files};
	my %stems;
	my @formatted;
	my @unformatted;
	foreach my $file (@$files_ref) {
		
		#$devtools->print_hash($file_ref);
		my @file_bits = split (/\./, $file);
		my $type = pop @file_bits;
		my $stem = join ('.', @file_bits);
		
		# SORTING THINGS OUT
		if ($type eq 'fasta' or $type eq 'fa') { 
			$stem .= ".$type";
		}
		if ($type eq 'fasta') { 
			$type = 'fa'; # standardise type for FASTA files
		} 
		
		#print "\n\t NAME $file:\t\t $stem";
		unless ($stem and $type) { 
			die "\n\n\t  filename structure error for '$file'\n\n";
		}
		unless  (
			$type eq 'fa'  or $type eq 'ol' or
			$type eq 'nnd' or $type eq 'nni' or 
			$type eq 'nsd' or $type eq 'nsi' or 
			$type eq 'nsq' or $type eq 'ntm' or
			$type eq 'nin' or $type eq 'nhr') {  
			print "\n\n\t  FILE '$file' IS OF UNKNOWN TYPE\n\n";
		}
		
		if ($stems{$stem}) { 
			my $files_ref = $stems{$stem};
			if ($files_ref->{$type}) { 
				print "\n\t NAME $file: $stem";
				#$devtools->print_hash($file_ref);
				#$devtools->print_hash($files_ref);
				die;
			}
			$files_ref->{$type} = 1;
			if ($type eq 'fa') {
				$files_ref->{stem_target} = $file;
			}
		}
		
		else {
			my %stem_files;
			$stem_files{$type} = 1;
			if ($type eq 'fa') {
				$stem_files{stem_target} = $file;
			}
			$stems{$stem} = \%stem_files;
		}
	}
	#$devtools->print_hash(\%stems);

	my @stems = keys %stems;
	foreach my $stem (@stems) {
		
		my $types_ref = $stems{$stem};
		#$devtools->print_hash($types_ref);
		my $nin = $types_ref->{nin};
		my $nsq = $types_ref->{nsq};
		my $nhr = $types_ref->{nhr};
		#my $nsi = $chunk_data->{nsi};
		#my $nsd = $chunk_data->{nsd};
		my $target = $types_ref->{stem_target};
		unless ($nhr and $nin and $nsq) {
			
			if ($target) {
				push (@unformatted, $target); 
			}
			else {
				print "\n\t No target for $stem";
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
# Subroutine:  format_genome
# Description: 
#***************************************************************************
sub format_genome {

	my ($self, $genome_ref) = @_;

	# Get data 
	my $genome_path   = $self->{genome_use_path};
	my $line_limit    = $self->{line_limit};
	my $blast_bin_dir = $self->{blast_bin_path};
	my $organism      = $genome_ref->{organism};
	my $type          = $genome_ref->{source_type};
	my $version       = $genome_ref->{version};
	my $version_path  = $genome_ref->{version_path};
	my $path          = $genome_path . $version_path;

	my $blast_program = 'makeblastdb';
	my $bin_path;
	if ($blast_bin_dir) {
		$bin_path = $self->{blast_bin_path} . $blast_program;
	}
	else {
		$bin_path = $blast_program;
	}

	# Iterate through the files
	my $formatted_ref   = $genome_ref->{formatted};
	my $unformatted_ref = $genome_ref->{unformatted};
	foreach my $file (@$unformatted_ref) { 
		
		# Get path to file
		my $chunk_path = $path . "/$file";
		my %data;
		#print "\n\t FILE $file: $chunk_path";
		$self->get_genome_chunk_stats(\%data, $chunk_path);
		my $total_bases    = $data{total_bases};
		my $line_count     = $data{total_lines};
		my $num_scaffolds  = $data{number_scaffolds};
		print "\n\t line count $line_count";
		#$devtools->print_hash(\%data); die;

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
       			my $command = "$bin_path -in $file_path -dbtype nucl> /dev/null";
				print "\n\t$command\n";
				system $command;
			}

			#$devtools->print_array(\@split_chunks); die;
			if ($done_split) {
				foreach my $new_chunk_path (@split_chunks) {
       				my $command = "$bin_path -in $new_chunk_path -dbtype nucl> /dev/null";
					print "\n\t$command\n";
					system $command;
				}
				my $file_path = $path . "/$file";
				my $command = "mv $file_path ./";
				print "\n\tMoving original file '$file_path' to ./\n";
				system $command;
			}
		}
		else {
       		my $command = "$bin_path -in $chunk_path -dbtype nucl> /dev/null";
			print "\n\t$command\n";
			system $command;
			push (@$formatted_ref, $file);
		}
	}

}

############################################################################
# Write output
############################################################################

#***************************************************************************
# Subroutine:   get_genome_chunk_stats
# Description:  
#***************************************************************************
sub get_genome_chunk_stats {

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
	
	#open REFORMATTED_CHUNK, ">$reformat_path" or die "\n\tCan't open $reformat_path\n";
	#print "\n\t # Paths:       $chunk_path: ";
	#print "\n\t # Line count:  $line_count";
	#print "\n\t # Base count:  $total_bases";
	#print "\n\t # Scaffolds:   $num_scaffolds";
	
	$data_ref->{total_bases}      = $total_bases;
	$data_ref->{total_lines}      = $line_count;
	$data_ref->{number_scaffolds} = $num_scaffolds;
}

############################################################################
# Splitting large files
############################################################################

#***************************************************************************
# Subroutine:  split_genome_chunk
# Description: 
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
	print "\n\t $command";
	system $command;
	
}

#***************************************************************************
# Subroutine:  split_longline_contig
# Description: reformat any contigs where the sequence is formatted as
#              one long line of text, so that they are split over multiple
#              lines (much easier to work with in a command line context)
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

#***************************************************************************
# Subroutine:  by number
# Description: by number - for use with perl 'sort' 
#***************************************************************************
sub by_number { $a <=> $b }	

############################################################################
# EOF
############################################################################
