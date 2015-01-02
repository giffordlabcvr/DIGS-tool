############################################################################
# Module:       SeqIO.pm 
# Description:  Provides fxns for reading and writing biological 
#               sequence data 
# History:      Rob Gifford, May 2009: Creation
############################################################################
package SeqIO;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::FileIO;
use Base::DevTools;
use Base::Sequence;

############################################################################
# Globals
############################################################################
	
my $fileio   = FileIO->new();
my $devtools = DevTools->new();
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Parameters
#***************************************************************************
sub new {

	my ($invocant) = @_;
	my $class = ref($invocant) || $invocant;

	# Member variables
	my $self = {
	
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# Reading Fxns 
############################################################################

#***************************************************************************
# Subroutine:  read_fasta
# Description: read a fasta file into an array of hashes. 
# Arguments:   $file: the name of the file to read
#              $array_ref: reference to the hash array to copy to
#***************************************************************************
sub read_fasta {

	my ($self, $file, $array_ref, $identifier) = @_;
	
	unless ($identifier) { $identifier = 'SEQ_'; }

	# Read in the file or else return
	unless (open(INFILE, $file)) {
		print "\n\t Cannot open file \"$file\"\n\n";
		return undef;
	}

	# Use process ID and time to create unique ID stem
	my $pid  = $$;
	my $time = time;
	my $alias_stem = $identifier;
	
	# Iterate through lines in the file
	my @raw_fasta = <INFILE>;
	close INFILE;
	my $header;
    my $sequence;
   	my $i = 0;
	foreach my $line (@raw_fasta) {
		
		#print "\n## $i";
		chomp $line;
		if    ($line =~ /^\s*$/)   { next; } # discard blank line
		elsif ($line =~ /^\s*#/)   { next; } # discard comment line 
		elsif ($line =~ /^BEGIN/)  { last; } # stop if we reach a data block
		elsif ($line =~ /^>/) {
			
			$line =~ s/^>//g;
					
			# new header, store any sequence held in the buffer
			if ($header and $sequence) {
				$i++;
				my $alias_id = $alias_stem . "_$i";
				$sequence = uc $sequence;
				my $seq_obj = Sequence->new($sequence, $header, $alias_id);
				#$devtools->print_hash($seq_obj); die;
				push(@$array_ref, $seq_obj);
			
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
		$i++;
		my $alias_id = $alias_stem . "_$i";
		$sequence =~ s/\s+//g; # Remove whitespace
		$sequence = uc $sequence;
		my $seq_obj = Sequence->new($sequence, $header, $alias_id);
		push(@$array_ref, $seq_obj);
	}
}

#***************************************************************************
# Subroutine:  convert_fasta
# Description: 
# Arguments: $fasta_ref: reference to array containing FASTA formatted data
#            $array_ref: reference to the hash array to copy to
#***************************************************************************
sub convert_fasta {

	my ($self, $fasta_ref, $array_ref) = @_;
	
	# Use process ID and time to create unique ID stem
	my $pid  = $$;
	my $time = time;
	my $alias_stem = $pid . '_' . $time . '_';
	
	# Iterate through lines in the file
	close INFILE;
	my $header;
    my $sequence;
   	my $i = 0;
	foreach my $line (@$fasta_ref) {
		
		chomp $line;
		#$line =~ s/\s+//g; # Remove whitespace
		if    ($line =~ /^\s*$/)   { next; } # discard blank line
		elsif ($line =~ /^\s*#/)   { next; } # discard comment line 
		elsif ($line =~ /^>/) {
			
			$line =~ s/^>//g;
					
			# new header, store any sequence held in the buffer
			if ($header and $sequence) {
				$i++;
				my $alias_id = $alias_stem . $i;
				$sequence = uc $sequence;
				#print "\n\t HERE '$header'";
				my $seq_obj = Sequence->new($sequence, $header, $alias_id);
				push(@$array_ref, $seq_obj);
			}
		
			# reset the variables 
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
		$i++;
		my $alias_id = $alias_stem . $i;
		$sequence = uc $sequence;
		#print "\n\t HERE '$header'";
		my $seq_obj = Sequence->new($sequence, $header, $alias_id);
		push(@$array_ref, $seq_obj);
	}
}

#***************************************************************************
# Subroutine:  write fasta 
# Description: write FASTA
#***************************************************************************
sub write_fasta {

	my ($self, $file, $sequences_ref) = @_;

	my @fasta;
	foreach my $sequence_ref (@$sequences_ref) {
	
		my $id       = $sequence_ref->{header};
		my $sequence = $sequence_ref->{sequence};
		my $fasta = ">$id\n$sequence\n\n";
		#print "\n$fasta\n";
		push (@fasta, $fasta);
	}
	$fileio->write_file($file, \@fasta);
}

#***************************************************************************
# Subroutine:  write delimited
# Description: write sequence data to tab-delimited
#***************************************************************************
sub write_delimited {

	my ($self, $file, $data_ref) = @_;

	# Create column headings & check consistency
	#$devtools->print_array($data_ref); die;
	my %fields;
	my $i = 0;
	foreach my $hash_ref (@$data_ref) {
		my @fields = keys %$hash_ref;
		foreach my $field (@fields) {
			if ($i > 1) {
				# Skip fields that are not defined in first row
				unless ($fields{$field}) { next; } 
			}
			else {
				$fields{$field} = 1;
			}
		}
	}

	# Set up the fields for writing the file
	my @fields = keys %fields;
	my @ordered_fields;
	push (@ordered_fields, 'sequence_id');
	foreach my $field (@fields) {
		if ($field eq 'sequence_id') { next; }
		elsif ($field eq 'sequence') { next; }
		else {
			push (@ordered_fields, $field);
		}
	}
	push (@ordered_fields, 'sequence');
	
	# Add column headings
	my $header_line = join("\t", @ordered_fields);
	my @data;
	push (@data, "$header_line\n");

	# Create the file
	foreach my $hash_ref (@$data_ref) {
		my @line;
		
		foreach my $field (@ordered_fields) {
			my $value = $hash_ref->{$field};
			unless ($value) {
				$value = 'NULL';
				if ($field eq 'sequence') {
					$value = 'NULL1';
					print "\n\t ## Warning - empty sequence field"; 
				}
				elsif ($field eq 'sequence_id') {
					$value = 'NULL2';
					print "\n\t ## Warning - empty seq ID field"; 
				}
			}
			push (@line, $value);
		}
		my $line = join("\t", @line);
		push (@data, "$line\n");
	}
	$fileio->write_file($file, \@data);
}

############################################################################
# END OF FILE 
############################################################################
