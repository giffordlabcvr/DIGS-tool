############################################################################
# Script:       SeqIO.pm 
# Description:  Provides fxns for reading and writing biological 
#               sequence data 
# History:      Rob Gifford, May 2009: Creation
# Details:      Sequences will ALWAYS be converted to uppercase as they are
#               read.
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
# Subroutine:  read_GLUE_MSA
# Description: 
#***************************************************************************
sub read_GLUE_MSA {

	my ($self, $file, $array_ref) = @_;
	
	# Read in the file or else return
	unless (open(INFILE, $file)) {
		print "\n\t Cannot open file \"$file\"\n\n";
		return undef;
	}

	# Use process ID and time to create unique ID stem
	my $pid  = $$;
	my $time = time;
	my $alias_stem = $pid . '_' . $time . '_';
	
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
		elsif ($line =~ /^>/) {
			
			$line =~ s/^>//g;
					
			# new header, store any sequence held in the buffer
			if ($header and $sequence) {
				$i++;
				my $alias_id = $alias_stem . $i;
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
		my $alias_id = $alias_stem . $i;
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
# Subroutine:  read_glue_msa 
# Description: read a GLUE MSA file into an array of hashes. 
# Arguments:   $file: the name of the file to read
#              $array_ref: reference to the hash array to copy to
#***************************************************************************
sub read_glue_msa {

	my ($self, $msa_path, $hash_ref) = @_;

	my @file = split(/\//,$msa_path);
    my $aln_name = pop(@file);
	open(IND, "$msa_path") || die "Cant open out vglue file\n";
	my $indb=0;
	my @ind;
	my $indel;
	my @headers;
	my $header;
	my $indel_st;
	my $i=0;

	# Get the indels 
	my %insertions; # 
	while(<IND>){

		if ($_ =~ /^\s*$/){next;}
		if ($_ =~ /^>(.+)(\s*)/){
			$headers[$i] = $1; #get the sequence headers into an array
			$i++;
			next;
		}
		if ($indb == 1){
			if($_ =~ /^ENDBLOCK/){last;}
			$indel = $_;
			chomp($indel);
			@ind = split(/\t/,$indel);
			$header = $ind[0];
			$indel_st = $ind[1];
			$insertions{$indel_st}->{$header} = $ind[2]; # get every inserion into a hash of hashes by position and header
			next;
		}
		if ($_ =~ /^BEGIN INSERTIONS/){
			$indb=1;
			next;
		}	
	}
	if ($indb == 0){
		die "\tNo INDELS found (Please make sure that $msa_path has at the end the INDELS)\n";
	}
	close(IND);
	# FINISH READING INSERTIONS

	# Get the sequences
	my @sequences; 
	$self->read_fasta($msa_path, \@sequences);

	$hash_ref->{sequences}  = \@sequences;
	$hash_ref->{insertions} = \%insertions;
}

############################################################################
# FORMAT CONVERSION
############################################################################

#***************************************************************************
# Subroutine:  fasta_to_phylip 
# Description: fasta to phylip format
#***************************************************************************
sub fasta_to_phylip {

	my ($self, $fasta_path, $phylip_ref, $table_ref) = @_;

	my @fasta;
	$self->read_fasta($fasta_path, \@fasta);
	#$devtools->print_array(\@fasta); die;
	
	my $set_id_len = 20;
	my $phycount;
	my $seq_len;
	foreach my $sequence_ref (@fasta) {
		
		$phycount++;
		my $header      = $sequence_ref->{header};
		my $sequence_id = $sequence_ref->{sequence_id};
		my $sequence    = $sequence_ref->{sequence};
		$seq_len     = length $sequence;
		unless ($seq_len) {
			$seq_len = length $sequence;
		}
		else {
			unless ($seq_len eq length $sequence) {
				die "\n\t sequence $sequence_id is a different length\n
				       \t check FASTA sequences are aligned\n\n";
			}
		}
	
		# Create phylip id
		my @sequence_id = split ('', $sequence_id);
		my $id_len;
		my $phy_id = $phycount . '_'; 
		foreach my $char (@sequence_id) {
			$phy_id .= $char;	
			$id_len = length($phy_id);
			if ($id_len eq 10) { last; }
		}
		my $spacer_len = ($set_id_len - $id_len);
		my $spacer = ' ' x $spacer_len;
		my $phy_seq = $phy_id . $spacer . $sequence . "\n";
		push(@$phylip_ref, $phy_seq);
		
		# store id relationship in translation table 
		$table_ref->{$phy_id} = $header;
		#$table_ref->{$sequence_id} = $phy_id;
	
	}
	
	# Create PHYLIP taxa and characters header
	my $num_taxa  = $phycount;
	my $num_chars = $seq_len;
	my $header_line = $num_taxa . '   ' . $num_chars . "\n";
	unshift(@$phylip_ref, $header_line);

}

############################################################################
# Alternative FASTA reaidng functions
############################################################################

#***************************************************************************
# Subroutine:  read_fasta_file_to_hash
# Description: read a fasta file into a hash, so that the headers 
#              become keys and the sequences the values
# Arguments:   $file: the name of the file to read
#              $hash_ref: reference to the hash array to copy to
#***************************************************************************
sub read_fasta_file_to_hash {

	my ($self, $file, $hash_ref) = @_;
	unless (open(INFILE, $file)) {
		print "\n\t Cannot open file \"$file\"\n\n";
		return undef;
	}

	my @raw_fasta = <INFILE>;
	close INFILE;

	my $header;
    my $sequence;
    foreach my $line (@raw_fasta) {

		if    ($line =~ /^\s*$/)   { next; } # discard blank line
		elsif ($line =~ /^\s*#/)   { next; } # discard comment line 
		elsif ($line =~ /^>/) {
			
			$line =~ s/^>//g;
					
			# new header, store any sequence held in the buffer
			if ($header and $sequence) {
				$sequence =~ s/\s+//g;
				$sequence =~ s/\n//g;
				$header   =~ s/\s+$//;
				$hash_ref->{$header} = uc $sequence;
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
		
		# Remove whitespace
		$sequence =~ s/\s+//g;
		$sequence =~ s/\n//g;
		$header   =~ s/\s+$//;
		chomp $header;	
		$hash_ref->{$header} = uc $sequence;
	}

	# to do: incorporate error checking in 'success'
	# e.g. do num seqs match num chevrons?
	my $success = 1;
	return $success;
}

#***************************************************************************
# Subroutine:  read_fasta_file_to_tab_delimited_array
# Description: read a fasta file into an array, as tab delimited data
# Arguments:   $file: the name of the file to read
#              $array_ref: reference to the hash array to copy to
#***************************************************************************
sub read_fasta_file_to_tab_delimited_array {

	my ($file, $array_ref) = @_;
	unless (open(INFILE, $file)) {
		print "\n\t Cannot open file \"$file\"\n\n";
		return;
	}

	my @raw_fasta = <INFILE>;
	close INFILE;

	my $header;
    my $sequence;
    foreach my $line (@raw_fasta) {

		chomp $line;
		$line =~ s/\s+//g; # Remove whitespace
		if    ($line =~ /^\s*$/)   { next; } # discard blank line
		elsif ($line =~ /^\s*#/)   { next; } # discard comment line 
		elsif ($line =~ /^>/) {
			
			$line =~ s/^>//g;
					
			# new header, store any sequence held in the buffer
			if ($header and $sequence) {
				my $uc_sequence = uc $sequence;
				my $data = "$header\t$uc_sequence\n";
				push(@$array_ref, $data);
			}
		
			# reset the variables 
			$line =~ s/^>//;
			$header = $line;
			$sequence = undef;
		}
		else {
			# keep line, add to sequence string
            $sequence .= $line;
     		$sequence =~ s/~/-/g; # Format gaps correctly
	 }
    }
	
	# Before exit, store any sequence held in the buffer
	if ($header and $sequence) {
		my $uc_sequence = uc $sequence;
		my $data = "$header\t$uc_sequence\n";
		push(@$array_ref, $data);
	}
}

#***************************************************************************
# Subroutine:  index_header_fields
# Description: 
#***************************************************************************
sub index_header_fields {

    my ($index, $header, $data_ref) = @_; 

    # Create an array of the fields in the order they occur in the data file
    my @fields;
    chomp $header;
    @fields = split ("\t", $header);
           
    # set sequence ID to a namespace defined heading ('sequence_id')
    shift @fields; # Remove the 1st column (should be the sequence id)
    unshift (@fields, 'sequence_id'); # Use 'sequence_id' as column heading
    my $field_index = 0;
    foreach my $field (@fields) {
        $data_ref->{$field_index} = $field;
        $field_index++;
    }   
}

############################################################################
# Writing Fxns 
############################################################################

#***************************************************************************
# Subroutine:  write fasta 
# Description: write sequences as FASTA
#***************************************************************************
sub write_fasta {

	my ($self, $file, $seqs_ref) = @_;

	my @fasta;
	my $id;
	my $sequence;
	foreach my $seq_ref (@$seqs_ref) {
		my $id       = $seq_ref->{sequence_id};
		my $header   = $seq_ref->{header};
		my $sequence = $seq_ref->{sequence};
		my $fasta     = ">$header\n$sequence\n\n";
		push (@fasta, $fasta);
	}
	
	$fileio->write_output_file($file, \@fasta);
}

############################################################################
# END OF FILE 
############################################################################
