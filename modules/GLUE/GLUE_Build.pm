#!/usr/bin/perl -w
############################################################################
# Script:       GLUE_Build 
# Description:  Build a GLUE analysis 
# History:      Rob Gifford, May 2014: Creation
############################################################################
package GLUE_Build;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::IO;
use Base::BioIO;
use Base::SeqIO;
use Base::FileIO;
use Base::DevTools;
use Base::Sequence;

# Component Classes
use GLUE::RefSeqParser;
use GLUE::RefSeq;
use GLUE::RefSeqAlignment;

############################################################################
# Globals
############################################################################

# Create Base objects
my $io       = IO->new();
my $fileio   = FileIO->new();
my $bioio    = BioIO->new();
my $seqio    = SeqIO->new();
my $devtools = DevTools->new();
my $console  = Console->new();
my $writer   = HTML_Utilities->new();

# For displaying alignment progress
my $increment = 50;

1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create new GLUE_Build program class 
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Store sequence IDs in the order they were submitted
	my @seq_ids;
	my $sequences = $parameter_ref->{sequences};
	foreach my $seq_ref (@$sequences) {
		my $seq_id = $seq_ref->{header};
		push (@seq_ids, $seq_id);
	}

	# Set parameters 
	my $self = {
		refseq_use_path => $parameter_ref->{refseq_use_path},
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# SECTION: COMMAND LINE BUILD FXNS
############################################################################

#***************************************************************************
# Subroutine:  setup_glue_process
# Description: set up for GLUE processing 
#***************************************************************************
sub setup_glue_process {

	my ($self, $glue_obj) = @_;

	# Get parameters and objects from self
	my $process_id  = $glue_obj->{process_id};
	my $output_type = $glue_obj->{output_type};
	my $output_path = $glue_obj->{output_path};
	unless ($process_id)  { die; } 
	unless ($output_path) { die; } 
	unless ($output_type) { die; } 

	# Create a unique directory to store the output of this process
	my $report_dir = $output_path . $process_id . '/';
	$fileio->create_unique_directory($report_dir, $output_type);
	$glue_obj->{report_dir} = $report_dir;

	# Create the reference sequence
	my $refseq_path  = $glue_obj->{refseq_path};
	my $parser_obj = RefSeqParser->new();
	my %params;
	$parser_obj->parse_refseq_flatfile($refseq_path, \%params);
	my $refseq = RefSeq->new(\%params);
	my $refseq_name = $refseq->{name};
	$glue_obj->{refseq} = $refseq; # Store the refseq
	# DEBUG #$refseq->describe();  die;  
	
	# Get the input sequences and associated data
	$self->parse_input_seqdata($glue_obj);
	
}

#***************************************************************************
# Subroutine:  parse_input_seqdata
# Description: 
#***************************************************************************
sub parse_input_seqdata {

	my ($self, $glue_obj) = @_;

	my @sequences;
	
	# Get sequence data
	if ($glue_obj->{query_fasta}) {
		my $raw_path = $glue_obj->{query_fasta};
		$seqio->read_fasta($raw_path, \@sequences);
		$glue_obj->{num_raw_seqs} = scalar @sequences;
		if ($glue_obj->{data_path}) {
			my $datatool = DataTool->new();
			my %results;
			my @data;
			$fileio->read_file($glue_obj->{data_path}, \@data);
			$datatool->link_sequences_and_data(\@sequences, \@data, \%results);
			$glue_obj->{raw_tab_data} = \@data;
		}
		$glue_obj->{sequences} = \@sequences;
	}
	# Get sequence data
	if ($glue_obj->{SequenceData}) {
		my $dataset_name = $self->{dataset_name};
		unless ($dataset_name) { $dataset_name = 'submitted_seqs.fas'; }
		my $raw_data = $glue_obj->{SequenceData};
		my @dataset = split("\n", $raw_data); 
		my @f_dataset;
		foreach my $line (@dataset) { 
			chomp $line;
			$line =~ s/\r//g;
			push (@f_dataset, "$line\n"); 
		}
		my $report_dir = $glue_obj->{report_dir};
		my $query_fasta_path = $report_dir . "/$dataset_name" . ".fas";
		$fileio->write_file($query_fasta_path, \@f_dataset);
		
		$seqio->read_fasta($query_fasta_path, \@sequences);
		$glue_obj->{num_raw_seqs} = scalar @sequences;
		if ($glue_obj->{data_path}) {
			die;
			my $datatool = DataTool->new();
			my %results;
			my @data;
			$fileio->read_file($glue_obj->{data_path}, \@data);
			$datatool->link_sequences_and_data(\@sequences, \@data, \%results);
			$glue_obj->{raw_tab_data} = \@data;
		}
		$glue_obj->{sequences} = \@sequences;
		#$glue_obj->{query_fasta} = $query_fasta_path;
		#$seqio->read_fasta($query_fasta_path, \@sequences);
		#$glue_obj->{num_raw_seqs} = scalar @sequences;
	}

	# Get input GLUE MSA
	if ($glue_obj->{glue_msa_path}) {
		my $datatool = DataTool->new();
		my $msa_path = $glue_obj->{glue_msa_path};
		unless ($msa_path) { die;}
		my @sequences;
		my %data;
		$datatool->read_glue_msa($msa_path, \@sequences, \%data);
		$glue_obj->{num_raw_seqs} = scalar @sequences;
		my $refseq = $glue_obj->{refseq};
		unless ($refseq) { die; }
		$refseq->{start} = $data{lowest_start};
		$refseq->{stop}  = $data{highest_stop};
		if ($glue_obj->{data_path}) {
			my $datatool = DataTool->new();
			my %results;
			my @data;
			$fileio->read_file($glue_obj->{data_path}, \@data);
			$datatool->link_sequences_and_data(\@sequences, \@data, \%results);
			$glue_obj->{raw_tab_data} = \@data;
		}
		# Create GLUE alignment object 
		#$devtools->print_array(\@sequences); die;
		$glue_obj->create_msa_object(\@sequences);
		unless ($glue_obj->{glue_msa_obj}) { die; }
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

	my ($self, $ctl_file_path, $glue_obj) = @_;
	
	# see above	
	my @file;
	$fileio->read_file($ctl_file_path, \@file);

	# Parse the 'SCREENDB' block
	$self->parse_screendb_block(\@file, $glue_obj);

	#my %data;
	my $start = 'BEGIN PARAMS';
	my $stop  = 'ENDBLOCK';
	$fileio->read_standard_field_value_block(\@file, $start, $stop, $glue_obj);
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
	$glue_obj->{targets} = \@targets;	

	# Get refseq list
	my @refseq_block;
	$start = 'BEGIN REFSEQLIST';
	$stop  = 'ENDBLOCK';
	$fileio->extract_text_block(\@file, \@refseq_block, $start, $stop);
	my @refseqs;
	my $refseqlist = undef;
	foreach my $line (@refseq_block) {
		
		if ($line =~ /^\s*#/)   { next; } # discard comment line 
		chomp $line;
		push (@refseqs, $line);	
		$refseqlist ='true';
	}
	if ($refseqlist) {
		$glue_obj->{refseqlist} = \@refseqs;	
	}
}

#***************************************************************************
# Subroutine:  parse_screendb_block
# Description: read an input file to get parameters for screening
#***************************************************************************
sub parse_screendb_block {

	my ($self, $file_ref, $glue_obj) = @_;
	
	# Parse the 'SCREENDB' block
	my $start = 'BEGIN SCREENDB';
	my $stop  = 'ENDBLOCK';
	my $db_block = $fileio->read_standard_field_value_block($file_ref, $start, $stop, $self);
	unless ($db_block)  {
		return;
	}
	
	# Get the 'SCREENDB' block values and validate
	my $db_name  = $self->{db_name};
	my $server   = $self->{mysql_server};
	my $user     = $self->{mysql_username};
	my $password = $self->{mysql_password};
	$glue_obj->{db_name}  = $db_name;
	$glue_obj->{server}   = $server;
	$glue_obj->{username} = $user;
	$glue_obj->{password} = $password;

	unless ($db_name)  {
		die "\n\t Control file error: 'db_name' undefined in 'SCREENDB' block\n\n\n";
	}
	unless ($server)  {
		die "\n\t Control file error: 'mysql_server' undefined in 'SCREENDB' block\n\n\n";
	}
	unless ($user)  {
		die "\n\t Control file error: 'mysql_username' undefined in 'SCREENDB' block\n\n\n";
	}
	unless ($password)  {
		die "\n\t Control file error: 'mysql_password' undefined in 'SCREENDB' block\n\n\n";
	}
}

############################################################################
# SECTION: WEB BUILD FXNS
############################################################################

#***************************************************************************
# Subroutine:  get_glue_form_params
# Description: initialise processs with parameters received via CGI 
#***************************************************************************
sub get_glue_form_params {

	my ($self, $query, $glue_obj, $glue_form) = @_;
	
	# Set defaults
	my $refseq_name;
	my $refseq_use_path = $self->{refseq_use_path};
	unless ($refseq_use_path) { die; }

	#$devtools->print_hash($query); die;

	$glue_obj->{output_type}  = 'html';
	my $dataset_name = $query->param('DatasetName');
	$glue_obj->{dataset_name} = $dataset_name;

	# print the formatted input page if no sequnce data has been received
	unless ( $query->param('SequenceData') and $query->param('ReferenceSequence')
		or   $query->param('TabDataFile')  and $query->param('ReferenceSequence')) {
		my @html;
		$fileio->read_file($glue_form, \@html);
		my $html  = join('', @html);
		print $html;
		exit;
	}

	# if we've got data collect it
	else {
		
		# Get pasted sequence data if its there
		my $raw_data;
		if ( $query->param('SequenceData') ) {
			$raw_data = $query->param('SequenceData');
 			$glue_obj->{SequenceData} = $raw_data;
    	}
		if ($query->param('glue_mode') eq 'create') {
			$self->{glue_mode} = 'create';
			if ($query->param('use_first_as_refseq')) {
 				$glue_obj->{use_first_as_refseq} = 'true';
            }
            else {
				# get reference sequence
				$refseq_name = $query->param('ReferenceSequence');
				$glue_obj->{reference_vglue} = $refseq_name;
				$glue_obj->{refseq_path}     = $refseq_use_path ."/$refseq_name";
			}
		}
	}

	my $glue_mode = $query->param('glue_mode');
	#print "REFSEQ: $refseq_name "; print "GLUE MODE: $glue_mode ";
	unless ($glue_mode and $refseq_name) { die; }; 

	if ( $query->param('TabDataFile') ) {
		my @raw_tab_data;
		my $file_name = $query->upload('TabDataFile');
		while (my $line = <$file_name>) { push (@raw_tab_data, $line); }
		$glue_obj->{raw_tab_data} = \@raw_tab_data;
	}
	if ($query->{phylogeny}) {
		$glue_obj->{phylogeny} = 'true';
	}
	if ($query->{glue_refset}) {
		my $path = $self->get_default_refset_path($refseq_name);
		unless ($path) { die "\n\t No reference set found for '$refseq_name'"; }
		$glue_obj->{glue_refset} = $path;
		#$glue_obj->{glue_refset} = 'true';
	}
	if ($query->{mut_profile}) {
		$glue_obj->{mut_profile} = 'true';
	}
	if ($query->{mut_frequencies}) {
		$glue_obj->{mut_frequencies} = 'true';
		my $array_ref = $query->{field_num};
		#$devtools->print_array($array_ref); die; 
		foreach my $line (@$array_ref) {
			chomp $line;
			$glue_obj->{field_num} = $line;
		}
	}
	if ($query->{show_mutlist_seqs}) {
		$glue_obj->{show_mutlist_seqs} = 'true';
	}
	if ($query->{show_synonymous_changes}) {
		$glue_obj->{show_synonymous_changes} = 'true';
	}

}

#***************************************************************************
# Subroutine:  get_default_refset_path
# Description: 
#***************************************************************************
sub get_default_refset_path { 

	my ($self, $refseq_name) = @_;

	# Load path to default reference set setting
	my %translations;
	#$translations{"HIV-1-NL43"} = "./projects/HIV-1/HIV-genomes.glu"; 
	$translations{"HIV-1"} = "./projects/HIV-1/HIV-1-HXB2-genomes.glu"; 
	$translations{"HCV"}   = "./projects/HCV/HCV-genomes.glu";
	$translations{"ZEBOV"} = "./projects/ZEBOV/ZEBOV-genomes.glu";
	$translations{"BTV"}   = "./projects/BTV/BTV-genomes.glu";
	$translations{"FMDV"}  = "./projects/FMDV/FMDV-genomes.glu";
	$translations{"CHIKV"} = "./projects/CHIKV/CHIKV-genomes-simpleheader.glu";
	$translations{"HRV"}   = "./projects/HRV/HRV-genomes.glu";
	$translations{"ASFV"}  = "./projects/ASFV/ASFV-genomes.glu";
	$translations{"HepE"}  = "./projects/HepE/HepE-genomes.glu";
	
	my $refseq_path =  $translations{$refseq_name};

	return $refseq_path;
	die;
}	

#***************************************************************************
# Subroutine:  by number
# Description: by number - for use with perl 'sort'  (cryptic but works) 
#***************************************************************************
sub by_number { $a <=> $b }	

############################################################################
# EOF
############################################################################
