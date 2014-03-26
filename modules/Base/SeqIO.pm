############################################################################
# Script:       SeqIO.pm 
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

############################################################################
# END OF FILE 
############################################################################
