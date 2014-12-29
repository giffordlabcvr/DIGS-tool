############################################################################
# Script:       BioIO.pm 
# Description:  Provides commonly used reading and writing services for 
#               bioinformatics, in particular, reading and writing sequence 
#               data in differengt formats 
# History:      Rob Gifford, May 2007: Creation
############################################################################
package BioIO;

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

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Member variables
	my $self = {
	
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# Newick tree parsing 
############################################################################

#***************************************************************************
# Subroutine:  back_translate_newick
# Description: use a translation table to back translate aliased newick tree
#              to the original taxa names
#***************************************************************************
sub back_translate_newick {

	my ($self, $newick, $trans_table_ref) = @_;

	my @newick = split('', $newick);
	my $taxon = '';
	my $translate_newick;
	my $in_branch_len = undef;
	my $in_taxon      = undef;
	#print $newick; die;

	foreach my $char (@newick) {
		
		chomp $char;
		if ($char eq '(') {
			$translate_newick .= $char;
		}
		elsif ($char eq ',' or $char eq ')' or $char eq ';') {
			$translate_newick .= $char;
			$in_branch_len = undef;
		}
		elsif ($char eq ':') {
			$in_branch_len = 'true';
			 if ($in_taxon) {
				my $translate = $trans_table_ref->{$taxon};
				#print "\ntaxon '$taxon': '$translate'";
				unless ($translate) {
					$devtools->print_hash($trans_table_ref);
					die "\n\t Nothing found for '$taxon'\n\n";
				}
				$translate_newick .= $translate;
				$in_taxon = undef;
				$taxon = '';
			}
			$translate_newick .= $char;
		}
		elsif ($in_branch_len) {
			$translate_newick .= $char;
		}
		else {
			$in_taxon = 'true';
			$taxon .= $char;
		}
	}
		
	return $translate_newick;

}

############################################################################
# NEXUS parsing 
############################################################################

#***************************************************************************
# Subroutine:  extract newick
# Description: 
#***************************************************************************
sub extract_newick {

	my ($self, $file_ref) = @_;

	# Extract the newick string
	my @newick;
	my $raw_newick = pop @$file_ref;
	my @raw_newick = split('', $raw_newick);
	my $newick;
	my $copy = undef;
	foreach my $char (@raw_newick) {
	
		if ($char eq '(' and !$copy) {
			$copy = 'true';
		}
		if ($char eq ';') {
			$newick .= $char;
			$copy = undef;
		}
		if ($copy) {
			$newick .= $char;
		}
	}
	return $newick;
}

#***************************************************************************
# Subroutine:  extract_features
# Description: 
#***************************************************************************
sub extract_features {

	my ($self, $file_ref, $features_ref) = @_;

	# Extract the taxa block
	$fileio->extract_text_block($file_ref, $features_ref, 'Features', 'Endblock;');
	#$devtools->print_array($features_ref);




}

############################################################################
# Nexus data parsing functions
############################################################################

#**************************************************************************
# Subroutine:  extract_nexus_type_blocks
# Description: does what it says 
#***************************************************************************
sub extract_nexus_type_blocks { 

	my ($self, $input_ref, $blocks_ref) = @_;

	my $extract = undef;
	my $block_name = undef;
	my @block;
	foreach my $line (@$input_ref) {
		
		if    ($line =~ /^\s*$/)   { next; } # discard blank line
		my $test_line = uc $line;
		if ($test_line =~ 'ENDBLOCK') { 
			$extract = undef;          
			my @blockcopy = @block;
			$blocks_ref->{$block_name} = \@blockcopy; 
			#print "\n\t blockname $block_name END";
			#print "\n\t block @block";
			@block = (); # empty the array for capturing blocks
		}
		elsif ($test_line =~ 'BEGIN') { 
			
			chomp $test_line;
			$test_line =~ s/BEGIN//g;
			$test_line =~ s/\s+//g;
			$block_name = $test_line;
			#print "\n\t blockname $block_name";
			$extract = 'true';         
		}
		elsif  ($extract) { 
			push (@block, $line); 
		}
	}
}


############################################################################
# Genbank file parsing
############################################################################

#***************************************************************************
# Subroutine:  parse_gb_to_refseq
# Description: 
#***************************************************************************
sub parse_gb_to_refseq {

	my ($self, $file, $data_ref) = @_;

    # read the file into an array
	my @file;
    $fileio->read_file($file, \@file);
	
	# Create the metadata section
	my %refseq_data;
	$self->create_metadata(\@file, \%refseq_data);
	
	# Create the features section
	$self->create_features(\@file, \%refseq_data);
	
	# Create sequence block
	my $sequence = $self->get_sequence(\@file);
	
	# Create a refseq file from the parsed data	
	my $metadata_block_ref = $refseq_data{'metadata_block'};
	my $feature_block_ref  = $refseq_data{'features_block'};
	my $sequence_block     = "\nORIGIN\n$sequence\n//\n\n"; 
	
	push (@$data_ref, @$metadata_block_ref);
	push (@$data_ref, @$feature_block_ref);
	push (@$data_ref, $sequence_block);
}

#***************************************************************************
# Subroutine:  
# Description:
#***************************************************************************
sub create_metadata {

	my ($self, $gb_file_ref, $data_ref) = @_;
	
	my %seqdata;
	$self->get_headsect($gb_file_ref, \%seqdata);
	my @header;
	push (@header, "Begin Metadata;\n");
	#$devtools->print_hash(\%seqdata);

	my $name_string = $seqdata{organism};
	my @string = split (/\s/, $name_string);
	#$devtools->print_array(\@string);
	my $name = shift @string;
	push (@header, "name = $name_string;\n");

	my $version_string = $seqdata{version};
	my @version = split (/\s/, $version_string);
	#$devtools->print_array(\@version);
	my $accession = shift @version;
	push (@header, "accession = $accession;\n");
	push (@header, "Endblock;\n\n");

	$data_ref->{metadata_block} = \@header;
}

#***************************************************************************
# Subroutine:  get_headsect  
# Description: Get all from the 'header section' of the Genbank file
#              The 'header section comprises everything before 'FEATURES'
#***************************************************************************
sub get_headsect {

	my ($self, $gb_file_ref, $data_ref) = @_;
	
	# Extract the relevant block
	my $start_token = 'LOCUS';	
	my $stop_token  = 'FEATURES';	
	my @block;
	$fileio->extract_text_block($gb_file_ref, \@block, $start_token, $stop_token);

	my $sequence    = '';
	my $in_sequence = undef;
	my $i = 0;
	foreach my $line (@block) {
		
		$i++;
		my $compress_line = $line;
		$compress_line =~ s/\s//g;
		#$data_ref->{$i} = $compress_line;
		chomp $line;
		
		if($compress_line =~ /^LOCUS/) {
			$line =~ s/^LOCUS\s*//;
			$data_ref->{locus} = $line;
		}
		elsif($compress_line =~ /^DEFINITION/) {
			$line =~ s/^\s*DEFINITION\s*//;
			$data_ref->{accession} = $line;
		}
		elsif($compress_line =~ /^ACCESSION/) {
			$line =~ s/^\s*ACCESSION\s*//;
			my @line = split (/\s+/, $line);
			my $acc = shift @line;
			$data_ref->{accession} = $acc;
		}
		elsif($compress_line =~ /^VERSION/) {
			$line =~ s/^\s*VERSION\s*//;
			$data_ref->{version} = $line;
		}
		elsif($compress_line =~ /^ORGANISM/) {
			$line =~ s/^\s*ORGANISM\s*//;
			$data_ref->{organism} = $line;
		}
	}

	#my @keys = keys %$data_ref;
	#my @sorted = sort by_number @keys;
	#foreach my $key (@sorted) {
	#	my $line = $data_ref->{$key};
	#	#print "$line\n";
	#}
}

#***************************************************************************
# Subroutine:  
# Description:
#***************************************************************************
sub create_features {

	my ($self, $file_ref, $data_ref) = @_;

	my %parsed_cds;
	$self->get_features($file_ref, \%parsed_cds);
	#$devtools->print_hash(\%parsed_cds);

	# Make features block
	my @features_block;
	push (@features_block, "Begin Features;\n");
	my @cd_keys = keys %parsed_cds;
	foreach my $cd_key (@cd_keys) {
		
		my $cd_ref = $parsed_cds{$cd_key};
		#$devtools->print_hash($cd_ref);
		
		my $gene      = $cd_ref->{'locus_tag'};
		my $gene_name = $cd_ref->{'product'};
		unless ($gene_name) {
			#my $note = $cd_ref->{'note'};
			#if ($note) { print "\n\t NOTE '$note'\n"; }
		}
		unless ($gene) {
			$gene = $cd_ref->{'gene'};
			unless ($gene) {
				#$devtools->print_hash($cd_ref);
				$gene = $gene_name;
			}
		}
		if ($gene and $gene_name) {
			my $line = "ORF\t$gene_name\t$gene\t";
			#print $line;
			my $cd_coordinates_ref = $cd_ref->{'cd_coordinates'};
			my $coordinates = '';
			if ($cd_coordinates_ref) {
				$coordinates = join ("\t", @$cd_coordinates_ref);
			}
			$line .= $coordinates;
			$line =~ s/\.\./\t/g;
			$line =~ s/\(/\t/g;
			$line =~ s/\)/\t/g;
			$line =~ s/complement/\t/g;
			push (@features_block, "$line\n");	
		}
	}	
	push (@features_block, "Endblock;\n\n");
	#$devtools->print_array(\@features_block);
	$data_ref->{'features_block'} = \@features_block;
}

#***************************************************************************
# Subroutine:  get_features  
# Description: Get data from the features section of a GenBank file
#***************************************************************************
sub get_features {

	my ($self, $file_ref, $parsed_cds_ref) = @_;
	
	# Extract the Features block
	my $start_token = 'FEATURES';	
	my $stop_token  = "^ORIGIN";	
	my @features;
	$fileio->extract_text_block($file_ref, \@features, $start_token, $stop_token);
	#$devtools->print_array(\@features);
	#exit;

	# Split out coding domain chunks
	my %split_cds;
	$self->split_out_coding_domains(\@features, \%split_cds);
	#$devtools->print_hash(\%split_cds);
	#exit;

	# Split out coding domain chunks
	my @cd_keys = keys %split_cds;
	foreach my $cd_key (@cd_keys) {
		my %parsed_cd;
		my $cd_string = $split_cds{$cd_key};
		$self->parse_coding_domain($cd_string, \%parsed_cd);
		$parsed_cds_ref->{$cd_key} = \%parsed_cd;
	}
	
	#$devtools->print_hash($parsed_cds_ref);
	#exit;
}	
	
#***************************************************************************
# Subroutine:  split_out_cds 
# Description: Get each codon domain section separately in an indexed hash 
#***************************************************************************
sub split_out_coding_domains {

	my ($self, $features_ref, $split_cds_ref) = @_;

	my $i = 0;
	my $cd_string = '';
	my $in_cd = undef;
	foreach my $line (@$features_ref) {
		
		#chomp $line;
		my @line = split(/\s+/, $line);
		if ($line[1] eq 'CDS') {
			$in_cd = 'true';
			$cd_string .= $line;
		}
		elsif ($line =~ /translation/) { # CD always ends with translation
			$i++;
			#print "\n\n$cd_string\n\n\n";
			$split_cds_ref->{$i} = $cd_string;
			$cd_string = '';
			$in_cd = undef;
		}
		elsif ($in_cd) {
			$cd_string .= $line;
		}
	}
}

#***************************************************************************
# Subroutine:  parse_coding_domain
# Description: 
#***************************************************************************
sub parse_coding_domain {

	my ($self, $cd_string, $parsed_cd_ref) = @_;

	#print "\n\n$cd_string\n\n\n";
	my @cd_block = split("\n", $cd_string);
	foreach my $line (@cd_block) {
		
		chomp $line;
		my @line = split(/\s+/, $line);
		if ($line[1] eq 'CDS') {
			my $range = $line[2]; 
			$range =~ s/\s//g;
			my @cd_coordinates;
			if ($range =~ 'join') {
				#print "\n\t join '$range'";
				$range =~ s/join//g;
				$range =~ s/\)//g;
				$range =~ s/\(//g;
				my @join_range = split(/,/, $range);
				foreach my $join_bit (@join_range) {
					#print "\n\t join '$join_bit'";
					push(@cd_coordinates, $join_bit);
				}
			}
			else {
				my @range = split(/\.\./, $range);
				my $start = $range[0]; 
				my $stop  = $range[1];
				push (@cd_coordinates, $start);
				push (@cd_coordinates, $stop);
			}
			$parsed_cd_ref->{'cd_coordinates'} = \@cd_coordinates;
		}
		elsif ($line =~ /gene=/) {
			$line =~ s/\s//g;
			$line =~ s/"//g;
			my @key_value_pair = split(/=/, $line);
			my $gene = $key_value_pair[1];
			$parsed_cd_ref->{'gene'} = $gene;
		}
		elsif ($line =~ /locus_tag=/) {
			$line =~ s/\s//g;
			$line =~ s/"//g;
			my @key_value_pair = split(/=/, $line);
			my $locus_tag = $key_value_pair[1];
			$parsed_cd_ref->{'locus_tag'} = $locus_tag;
		}
		elsif ($line =~ /codon_start=/) {
			$line =~ s/\s//g;
			$line =~ s/"//g;
			my @key_value_pair = split(/=/, $line);
			my $codon_start = $key_value_pair[1];
			$parsed_cd_ref->{'codon_start'} = $codon_start;
		}
		elsif ($line =~ /note=/) {
			my @key_value_pair = split(/=/, $line);
			my $codon_start = $key_value_pair[1];
			$parsed_cd_ref->{'note'} = $codon_start;
		}
		elsif ($line =~ /product=/) {
			$line =~ s/"//g;
			my @key_value_pair = split(/=/, $line);
			my $product = $key_value_pair[1];
			$parsed_cd_ref->{'product'} = $product;
		}
		elsif ($line =~ /protein_id=/) {
			$line =~ s/\s//g;
			$line =~ s/"//g;
			my @key_value_pair = split(/=/, $line);
		}
	}
}

#***************************************************************************
# Subroutine:  get_sequence
# Description: 
#***************************************************************************
sub get_sequence {

	my ($self, $file_ref, $f_width) = @_;

	# Set default width for sequence if width not set
	unless ($f_width) { $f_width = 80; }

	# Extract seq matrix
	my @sequence;
	$fileio->extract_text_block($file_ref, \@sequence, "^ORIGIN", '//');
	#$devtools->print_array(\@sequence);
	my $sequence = join("\n", @sequence);
	$sequence =~ s/\d//g; # remove numbers
    $sequence =~ s/[\s\/]//g;
	$sequence = uc $sequence; # We like sequences to be upper case

	# Format width
	my @f_sequence = split('', $sequence); 
	my $i = 0;
	my $f_sequence;
	foreach my $char (@f_sequence) {
		$i++;
		$f_sequence .= $char;
		my $modulo = $i % $f_width;
		if ($modulo eq 0) { # Insert a line break
			$f_sequence .= "\n";
		} 
			
	}
	return $f_sequence;
}


#***************************************************************************
# Subroutine:  by number
# Description: by number - for use with perl 'sort'  (cryptic but works) 
#***************************************************************************
sub by_number { $a <=> $b }	

############################################################################
# DEPRECATED 
############################################################################


############################################################################
# END OF FILE
############################################################################


############################################################################
# END OF FILE 
############################################################################
