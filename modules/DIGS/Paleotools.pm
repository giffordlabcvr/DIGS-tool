#!/usr/bin/perl -w
############################################################################
# Module:      Paleotools.pm
# Description: Tools for excavating viral fossils using screening results 
# History:     January 2012: Created by Robert Gifford 
############################################################################
package Paleotools;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;
use DBI;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::FileIO;
use Base::SeqIO;
use Base::DevTools;
use Base::Console;
use Base::Text_Utilities;

# Third party program interface modules
use Component::Interface::BLAST;   # Interface to BLAST
use Component::GLUE::Sequence;     # Basic sequence manipulations

############################################################################
# Globals
############################################################################

# Create base objects
my $fileio    = FileIO->new();
my $seqio     = SeqIO->new();
my $devtools  = DevTools->new();
my $console   = Console->new();
my $textutils = Text_Utilities->new();
my $html_utils = HTML_Utilities->new();

# Alignments directory
my $alignments_dir = "db/alignments/";
my $genome_path = '/genomes/macrolineage/'; 
my $env_threshold = 50; # for env assign
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

	# Set member variables
	my $self = {
		
		# Member variables
		mode                 => $parameter_ref->{mode},
		process_id           => $parameter_ref->{process_id},
		
		# Paths and constants
		process_id            => $parameter_ref->{process_id},
		output_type           => $parameter_ref->{output_type},
		tmp_path              => $parameter_ref->{tmp_path},
		blast_db_path         => $parameter_ref->{blast_db_path},
		blast_orf_lib_path    => $parameter_ref->{blast_orf_lib_path},
		blast_genome_lib_path => $parameter_ref->{blast_genome_lib_path},
		refseq_lib_path       => $parameter_ref->{refseq_lib_path},
		refseq_use_path       => $parameter_ref->{refseq_use_path},
		rss_path              => $parameter_ref->{rss_path},
		output_path           => $parameter_ref->{output_path},
		
		# Flags
		output_type          => $parameter_ref->{output_type},
		
		# Paths
		output_path          => $parameter_ref->{output_path},
		header_path          => $parameter_ref->{header_path},
		tmp_path             => $parameter_ref->{tmp_path}, 
		refseq_use_path      => $parameter_ref->{refseq_use_path}, 
		alignment_path       => $parameter_ref->{alignment_path}, 
		
		# Component objects
		vglue_obj            => $parameter_ref->{vglue_obj},
		refseq_obj           => $parameter_ref->{refseq_obj},
		blast_obj            => $parameter_ref->{blast_obj},
		muscle_obj           => $parameter_ref->{muscle_obj},
		seqgen_obj           => $parameter_ref->{seqgen_obj},
		water_obj            => $parameter_ref->{water_obj},
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# EXCAVATE OPTIONS 
############################################################################

#***************************************************************************
# Subroutine:  run_excavate_function 
# Description: just a handler, handing off to various BLAST routines
#              extract loci from loci table and save them into fasta file
#***************************************************************************
sub run_excavate_function1 {

	my ($self, $ctl_file) = @_;
	
	unless ($ctl_file) { die; }
	# Show title
	$console->refresh();
	my $title       = 'VGLUE Excavation tools';
	my $version     = '2.0';
	my $description = 'Sequence Database Mining using the GLUE framework';
	my $author      = 'Robert J. Gifford';
	my $contact		= '<robert.gifford@glasgow.ac.uk>';
	$console->show_about_box($title, $version, $description, $author, $contact);
	
	# Read the control file
	my $loader_obj = Utility->new($self);
	$loader_obj->parse_control_file($ctl_file, $self);
	my $path_stem = $self->{genome_path_stem};

	# Get the mode after parsing control file
	my $excavate_mode = $self->{excavate_mode};
	unless ($excavate_mode) { 
		die "\n\t Error parsing control file\n\n";
	}
	# Hand off to functions as specified
	if ($excavate_mode eq 2) { # Try to find LTRs de novo
		print "\n\t # Search for LTRs de novo";
		unless ($path_stem) { die "\n\t No genome path stem specified in ctl file\n\n"; }
		$self->search_for_ltrs($path_stem);
	}
	elsif ($excavate_mode eq 3) { # FIND PBS
		print "\n\t # Search for PBS";
		unless ($path_stem) { die "\n\t No genome path stem specified in ctl file\n\n"; }
		$self->pbs_search($path_stem);
	}
	#elsif ($excavate_mode eq 1) { # Reassess a reference sequence based on GLUE align
	#	print "\n\t # Reassess reference sequence";
	#	$self->reassess_reference_seqs();
	#}
	#elsif ($excavate_mode eq 7) {  # Search for ortholgs 
	#	print "\n\t # Search for orthologs";
	#	$self->search_for_orthologs();
	#	die;
	#}
	else {
		die "\n\t Error parsing control file\n\n";
	}	
}	

#***************************************************************************
# Subroutine:  run_excavate_function 
# Description: just a handler, handing off to various BLAST routines
#              extract loci from loci table and save them into fasta file
#***************************************************************************
sub run_excavate_function2 {

	my ($self, $excavate_mode, $db_name) = @_;
	
	# Show title
	$console->refresh();
	my $title       = 'VGLUE Excavation tools';
	my $version     = '2.0';
	my $description = 'Sequence Database Mining using the GLUE framework';
	my $author      = 'Robert J. Gifford';
	my $contact		= '<rgifford@adarc.org>';
	$console->show_about_box($title, $version, $description, $author, $contact);

  	# USAGE statement	
	my $USAGE  .= "\n\t ### a  = Paleotools - general processes\n";
  	$USAGE  .= "\n\t -a=1 -d=[db]  = classify envelopes"; 
  	$USAGE  .= "\n\t -a=2 -d=[db]  = classify envelopes"; 
 	$USAGE  .= "\n\n";
	unless ($excavate_mode and $db_name) { die "\n\t $USAGE\n\n"; }

	# Get the sequences from the table
	my $db_obj = DB->new();
	$db_obj->load_screening_db($db_name);
	$self->{db_obj} = $db_obj;
	
	# Hand off to functions as specified
	if ($excavate_mode eq 1) { # Reassess a reference sequence based on GLUE align
		print "\n\t # Searching for/classifying envelopes";
		$self->classify_envelopes();
	}
	elsif ($excavate_mode eq 2) { # Try to find ORFs de novo
		print "\n\t # Search for ORFs";
		#$self->search_for_orfs($path_stem);
	}
	elsif ($excavate_mode eq 3) {  # Search for ortholgs 
		print "\n\t # Do GLUE alignments";
		$self->align_loci();
	}
	else {
		die "\n\t $USAGE \n\n";
	}	
}	

############################################################################
# Control file guided Paleotools functions 
############################################################################

#***************************************************************************
# Subroutine:  align_loci
# Description: 
#***************************************************************************
sub align_loci {

	my ($self) = @_;

	# Get settings from the ctl file
	my $refseq  = $self->{refseq};
	my $db_name = $self->{db_name};
	unless ($refseq and $db_name) { die "\n\t Control file error\n\n"; }

	# Get the sequences from the table
	my $db_obj = DB->new();
	$db_obj->load_screening_db($db_name);
	$self->{db_obj} = $db_obj;

	# GET THE LOCUS COORDINATES FROM THE DB
	my $where  = "WHERE Genome_structure = 'X-LTR-X' AND assigned_to = 'Ornithorhynchus-ERV-R'";
	$self->select_from_loci_table($where);

	# GET THE SEQUENCES TO ALIGN
}

############################################################################
# SEARCH FOR LTRs
############################################################################

#***************************************************************************
# Subroutine:  search_for_ltrs
# Description: 
#***************************************************************************
sub search_for_ltrs {

	my ($self, $path_stem) = @_;

	# Get the DB name from self
	my $db_name = $self->{db_name};
	
	# Now get the loci 
	my $db_obj = DB->new();
	$db_obj->load_screening_db($db_name);
	my $ltr_table = $db_obj->{ltr_candidate_table};
	$ltr_table->flush();
	$ltr_table->reset_primary_keys();

	# Get the extracted table reference names
	my @refseqs;
	my $extracted_table = $db_obj->{extracted_table};
	my @fields = qw [ assigned_to ];
	$extracted_table->select_distinct(\@fields, \@refseqs);
	
	# Get the loci table reference names
	my @loci;
	my $loci_table = $db_obj->{loci_table};
	$loci_table->select_distinct(\@fields, \@loci);
	my %loci; # Index the loci in a hash by the reference name
	foreach my $locus_refseq_ref (@loci) {
		my $locus_refseq_name = $locus_refseq_ref->{assigned_to};;
		$loci{$locus_refseq_name} = 1;
	}
	
	# Get chunk paths
	my %chunk_paths;
	$db_obj->get_chunk_paths(\%chunk_paths, $path_stem);
	#$devtools->print_hash(\%chunk_paths);

	# Iterate through the loci
	my %ltr_results;
	foreach my $refseq_ref (@refseqs) {
		my $refseq_name = $refseq_ref->{assigned_to};;
		#my $orientation = $refseq_ref->{orientation};;
		if ($loci{$refseq_name}) { next; } # Skip refseqs that are in the Loci table
		print "\n\t Looking for LTRs in $refseq_name";
		$self->locate_LTRs($db_obj, $refseq_ref, \%chunk_paths);
	}
}

#***************************************************************************
# Subroutine:  locate_LTRs
# Description: locate paired LTRs
#***************************************************************************
sub locate_LTRs {

	my ($self, $db_obj, $data_ref, $chunk_paths) = @_; 

	# Get the name of the retroviral reference sequence
	my $assigned_to = $data_ref->{assigned_to};
	if ($assigned_to =~ /ERV-L/) {
		return;
	}	
	
	# Create objects
	my $blast_obj   = $self->{blast_obj};
	my $process_id  = $self->{process_id};
	my $tmp_path    = $self->{tmp_path};
	my $blast_bin_path = $blast_obj->{blast_bin_path};

	# GET/SET  PARAMS
	my $buffer = 5000; # how many flanking bases to extract

	# Get the extracted
	my @data;
	my @fields = qw [ assigned_to record_id extract_start extract_end 
	                  organism orientation chunk_name scaffold  ];
	my $where =  " WHERE Assigned_to = '$assigned_to' ";
	   $where .= " AND   Assigned_to_gene = 'Pol' ";
	my $extracted_table = $db_obj->{extracted_table};
	$extracted_table->select_rows(\@fields, \@data, $where);
	foreach my $row_ref (@data) {

		# Look for LTRs in this refseq
		my $record_id     = $row_ref->{record_id};
		my $organism      = $row_ref->{organism};
		my $extract_start = $row_ref->{extract_start};
		my $extract_end   = $row_ref->{extract_end};
		my $orientation   = $row_ref->{orientation};
		my $chunk_name    = $row_ref->{chunk_name};
		my $scaffold      = $row_ref->{scaffold};
		
		# Adjustments for scaffold name
		my $truescaf='';
		my @gi = split(/\|/, $scaffold);
		if (scalar(@gi) > 1) { $truescaf = $gi[1]; }
		else                 { $truescaf = $scaffold;  }
		
		# Extract the sequence
		print "\n\t Looking for LTRs in assigned hit: $record_id";

		# Extract a sequence that may contain flanks
		# Extract command example: 
		# /bin/blast/blastdbcmd -db hs_alt_HuRef_chrX.fa -entry 157734237 
		# -range 10-60 -strand minus
		my $start = $extract_start - $buffer;
		my $end   = $extract_end   + $buffer;
		if ($start < 1) { 
			print "\n\t # Setting start '$start' to '1'";
			#next;
			$start = 1;
		}

		if (($end - $start) < ($buffer * 2)) { 
			print "\n\t # Too short, skipping";
			next; 
		}
		if ($start > $end) { die; }
		my $target_path = $chunk_paths->{$chunk_name};
		my $command = $blast_bin_path . "blastdbcmd -db $target_path";
		$command .= " -entry $truescaf ";
		$command .= " -range $start-$end";
		if ($orientation eq '-ve') { $command .= ' -strand minus '; }
		my @sequence = `$command`; # Execute the command
		
		shift @sequence;  # Remove header
		my $sequence = join ('', @sequence);
		$sequence =~ s/\n//g;
		
		#$devtools->print_hash(\%loci_by_scaffold);
		my $seq_len = length $sequence;
		unless ($seq_len) {
			print "\n\t\t no hit for 'range $start-$end' skipping....";
			next;
		}
		my $f_sequence = ">$record_id\n$sequence\n";
		#print "$f_sequence"; die;
		my $file1_path  = $tmp_path . $process_id . '_ltr_search_in1.fas';
		$fileio->write_text_to_file($file1_path, $f_sequence); 
		my $file2_path  = $tmp_path . $process_id . '_ltr_search_in2.fas';
		$fileio->write_text_to_file($file2_path, $f_sequence); 
		my $result_path = $tmp_path . $process_id . '_ltr_result_out.fas';
		
		# Create and execute the BLAST command
		my $method_path  = $blast_bin_path . 'blastn';
		$command  = "$method_path -query $file1_path -subject $file2_path ";
		   $command .= " -out $result_path -outfmt 5 -task blastn"; # XML output
		system $command;		
		#print $command;

		# Parse out the alignment
		my @data;
		my $valid = $blast_obj->parse_xml_format_step_one($result_path, \@data);
		#$devtools->print_array(\@data); die;
		unless ($valid) { 
			print "\n\t\t # INVALID BLAST RESULTS, skipping....";
			next;
		}
		
		# Do LTR parsing step
		$self->parse_ltr_results($db_obj, \@data, $row_ref);
		
		# Clean up
		$command = "rm $result_path";
		system $command;
		$command = "rm $file1_path";
		system $command;
		$command = "rm $file2_path";
		system $command;
	}
}

#***************************************************************************
# Subroutine:  parse_ltr_results 
# Description: 
#***************************************************************************
sub parse_ltr_results {

	my ($self, $db_obj, $data_ref, $extracted_row_ref) = @_; 

	my @ltr_candidates;
	shift @$data_ref; # remove first hit
	my $seq_obj = Sequence->new();
	foreach my $hsp_ref (@$data_ref) {

		my $alnseq = $hsp_ref->{hsp_qseq};
		my @alnseq = split ('', $alnseq);
		my $query_start  = $hsp_ref->{query_start};
		my $query_stop   = $hsp_ref->{query_stop};
		my $aln_stop     = $hsp_ref->{aln_stop};
		my $aln_start    = $hsp_ref->{aln_start};
		my $aln_seq      = $hsp_ref->{hsp_hseq};
		my $query_seq    = $hsp_ref->{hsp_qseq};
		unless ($aln_stop and $query_stop and $aln_start and $query_start) { die; }	
		
		# What orientation are the matches in, relative to one another?
		# What position are matches in, relative to one another?
		my $highest_start;
		my $highest_stop;
		my $highest_seq;
		my $lowest_start;
		my $lowest_stop;
		my $lowest_seq;
		my $highest;
		
		my $orientation  = $extracted_row_ref->{orientation};
		if ($aln_start > $query_start and $orientation eq '+ve'
	       or  $aln_start < $query_start and $orientation eq '-ve') {
			$highest_start = $query_start;
			$highest_stop  = $query_stop;
			$highest_seq   = $query_seq;
			$lowest_start  = $aln_start;
			$lowest_stop   = $aln_stop;
			$lowest_seq    = $aln_seq;
			$highest       = 'query';
		}
		elsif ($aln_start > $query_start and $orientation eq '+ve'
	       or  $aln_start < $query_start and $orientation eq '-ve') {
			$highest_start = $aln_start;
			$highest_stop  = $aln_stop;
			$highest_seq   = $aln_seq;
			$lowest_start  = $query_start;
			$lowest_stop   = $query_stop;
			$lowest_seq    = $query_seq;
			$highest       = 'aln';
		}
		else { next; }	
		
		# Calculate stats for this match
		my $length = $highest_stop - $highest_start; 
		my $intervening = $highest_start - $lowest_stop;
		if ($orientation ne '-ve') {
			$intervening = $intervening * -1;
			#$lowest_seq  = $seq_obj->reverse_and_complement($lowest_seq);
			#$highest_seq = $seq_obj->reverse_and_complement($highest_seq);
		}
		# Skip if conditions not met
		if ($length < 300) { 
			print "\n\t\t # match too short ($length), skipping....";
			next; 
		}
		if ($intervening < 4000) { 
			print "\n\t\t # Match in $orientation orientation intervening too short ($intervening), skipping....";
			next; 
		}

		
		# Store the info
		my %data;
		$data{extract_id}           = $extracted_row_ref->{record_id};	
		$data{organism}             = $extracted_row_ref->{organism};	
		$data{assigned_to}          = $extracted_row_ref->{assigned_to};	
		$data{chunk_name}           = $extracted_row_ref->{chunk_name};	
		$data{scaffold}             = $extracted_row_ref->{scaffold};	
		$data{orientation}          = $extracted_row_ref->{orientation};	
		$data{five_prime_start}     = $lowest_start;
		$data{five_prime_end}       = $lowest_stop;
		$data{five_prime_sequence}  = $lowest_seq;	
		$data{three_prime_start}    = $highest_start;
		$data{three_prime_end}      = $highest_stop;
		$data{three_prime_sequence} = $highest_seq;	
		$data{intervening}          = $intervening;
		$data{aln_len}              = $length;
		#$devtools->print_hash(\%data); exit;
		push (@ltr_candidates, \%data);
	}
	
	# INSERT LTR CANDIDATE DATA
	#$devtools->print_array(\@ltr_candidates); 
	my $ltr_table = $db_obj->{ltr_candidate_table};
	foreach my $candidate_ref (@ltr_candidates) {
		$ltr_table->insert_row($candidate_ref);
	}
	#die;	
}

############################################################################
# PBS search 
############################################################################

#***************************************************************************
# Subroutine:  pbs_search
# Description: 
#***************************************************************************
sub pbs_search {
	
	my ($self, $path_stem) = @_;

	# Get data from self 
	$self->show_title();

	# Load the screening database
	my $db_name = $self->{db_name};
	unless ($db_name) {die; }
	my $db_obj = DB->new();
	$db_obj->load_screening_db($db_name);
	$self->{screening_db} = $db_obj;

	# Read in all the tRNA sequences
	my @tRNAs;
	#my $seqpath = './db/extra/tRNA_human.fas';
	#my $seqpath = './db/extra/Loxodonta_tRNA.fa';
	my $seqpath = $self->{tRNA_database};
	unless ($seqpath) {die; }
	$seqio->read_fasta($seqpath, \@tRNAs);

	# Flush PBS table
	my $pbs_table = $db_obj->{pbs_candidate_table};
	$pbs_table->flush();
	$pbs_table->reset_primary_keys();

	# Get the refseqs
	my $refseq_use_path = $self->{refseq_use_path};
	my $loci_table = $db_obj->{loci_table};
	my @refseqs;
	$loci_table->select_distinct_single_field('assigned_to', \@refseqs);

	# Get paths
	my %chunk_paths;
	$db_obj->get_chunk_paths(\%chunk_paths, $path_stem);
	
	# Iterate through refseqs	
	foreach my $refseq_name (@refseqs) {
		
		# Parse the input file and create the reference sequence
		#unless ($refseq_name =~ /ERV-9/) { next; }
		#unless ($refseq_name =~ /HERV-B7H6/) { next; }
		print "\n\t ### Looking for PBS in '$refseq_name'";
		my $parser_obj = RefSeqParser->new();
		my %params;
		my $path = $refseq_use_path . "$refseq_name";
		$parser_obj->parse_refseq_flatfile($path, \%params);
		my $refseq = RefSeq->new(\%params);
	
		# Get the leader
		my %utrs;
		$refseq->get_utrs(\%utrs);
		my $leader = $utrs{LEA};
		my $ltr    = $utrs{LTR};
		unless ($ltr and $leader) { 
			print "\n\t No leader, skipping";
			next;
		}
		
		# Get the leader end
		my $leader_end;
		my $features_ref = $refseq->{features};
		foreach my $feature_ref (@$features_ref) {
			my $gene_name  = $feature_ref->{name};
			if ($gene_name eq 'LEA') {
				$leader_end  = $feature_ref->{stop};
			}
		}	
		
		my $leader_seq = $ltr . $leader;
		my @erv_leaders;
		my %ref_leader;
		$ref_leader{locus_id}    = $refseq_name . '_REFSEQ';
		$ref_leader{sequence}    = $leader_seq;
		$ref_leader{assigned_to} = $refseq_name;
		push(@erv_leaders, \%ref_leader);
		$self->get_leader_loci($refseq, \%chunk_paths, \@erv_leaders);

		# Get all the loci with leaders for this sequence
		foreach my $erv_seq_ref (@erv_leaders) {

			# Write to file
			my $name = $erv_seq_ref->{locus_id};
			my $seq  = $erv_seq_ref->{sequence};
			my $lib_path = './site/tmp/tmp_leader.fas';
			my $leader_fasta = ">$name\n$seq";
			$fileio->write_text_to_file($lib_path, $leader_fasta);
			$erv_seq_ref->{fasta_path} = $lib_path;

			# Do the WATER PBS search
			$self->water_pbs_search($erv_seq_ref, $leader_end, \@tRNAs);
		}
	}
}

#***************************************************************************
# Subroutine:  water_pbs_search
# Description: 
#***************************************************************************
sub water_pbs_search {
	
	my ($self, $seq_ref, $leader_end, $tRNA_ref) = @_;

	# Get data from self 
	my $tmp_path    = $self->{tmp_path};
	my $water_obj   = $self->{water_obj};
	my $seq_obj     = Sequence->new();
	my $refseq_name = $seq_ref->{assigned_to};
	my $leader_path = $seq_ref->{fasta_path};
	
	# PBS results
	my @pbs_results;
	foreach my $tRNA_seqref (@$tRNA_ref) {
			
		# Get tRNA seq data
		my $header    = $tRNA_seqref->{header};
		my $query_seq = $tRNA_seqref->{sequence};
		unless ($query_seq) {die;} 
		
		# Get the 'acceptor stem' region from the tRNA sequence
		$query_seq = $seq_obj->reverse_and_complement($query_seq);
		my @query = split('', $query_seq);
		my $tail = 'TGG';
		my $i;
		do {
			$tail .= shift(@query);
			$i++;
		} until ($i eq 20);
		$query_seq = $tail;

		# VALIDATION with known acceptor stems
		#if ($query_seq =~ /TGGTGCCGTGACTCGGAT/) {
		#	print "\n# Histidine $header: $tail";
		#}
		#if ($query_seq =~ /TGGTTCCCTGACTGGGAA/) {
		#	print "\n# Glutamic acid $header: $tail";
		#}

		# Set up for the alignment	
		my $fasta1      = ">$header\n$query_seq";
		my $probe_path  = $tmp_path . 'tmp_nt1.txt';
		$fileio->write_text_to_file($probe_path, $fasta1);
		
		# Run and parse smith-waterman align 
		my $result_path  = $tmp_path . 'PBS_result.txt';
		$water_obj->run_water($probe_path, $leader_path, $result_path);
		my %data;
		$water_obj->parse_water($result_path, \%data);
		my $length_of_match = $data{length_of_match};
		my $identity        = $data{identity};

		# Ignore short matches
		unless ($length_of_match)             { next; } 
		if ($length_of_match < 15)            { next; } 
		if ($data{query_start} > $leader_end) { next; }
		
		my $divergence = (100-$identity)/100;
		$divergence = sprintf("%.3f", $divergence);
		
		$data{div} = $divergence;
		$data{PBS} = $header;
		$data{assigned_to} = $seq_ref->{assigned_to};
		$data{locus_id}    = $seq_ref->{locus_id};
		#print "\nMatch to $refseq_name: $header is $div ";
		#print "with match of length $length_of_match";
		#if ($header eq 'Homo_sapiens_chr1.trna43-GlyGCC') { die; }	
		push (@pbs_results, \%data);
	}
	
	# Get the best match
	my %best_match;	
	$water_obj->sort_water_pbs_results(\@pbs_results, \%best_match);
	my $filename = $tmp_path . "/$refseq_name"  . '_water_pbs_results';
	$fileio->write_output_file($filename, \@pbs_results);

	# INSERT NEW DATA CREATE NEW DB OBJ AVOID DISCONNECT 
	my $db_name = $self->{db_name};
	my $db_obj_tmp = DB->new();
	$db_obj_tmp->load_screening_db($db_name);
	my $pbs_table = $db_obj_tmp->{pbs_candidate_table};
	if ($best_match{pbs_assign}) {
		#$devtools->print_hash(\%best_match);
		$pbs_table->insert_row(\%best_match);
	}
}

############################################################################
# SEARCH FOR ORFs
############################################################################

#***************************************************************************
# Subroutine:  search_for_orfs
# Description: 
#***************************************************************************
sub search_for_orfs {

	my ($self, $path_stem) = @_;

	# Get the DB name from self
	my $db_name = $self->{db_name};
	
	# Now get the loci 
	my $db_obj = DB->new();
	$db_obj->load_screening_db($db_name);

	# Populate a hash with paths to the sequence DB files
	my %chunk_paths;
	$db_obj->get_chunk_paths(\%chunk_paths, $path_stem);
	$self->{chunk_paths} = \%chunk_paths;
	
	# Get the extracted table reference names
	my @refseqs;
	my $loci_table = $db_obj->{loci_table};
	$loci_table->select_distinct_single_field('assigned_to', \@refseqs);
	
	# Iterate through the loci
	foreach my $refseq_name (@refseqs) {
		print "\n\t Looking for ORFs in $refseq_name";
		#unless ($refseq_name eq 'HERV-K-HML2') { next; }
		$self->derive_orfs($db_obj, $refseq_name, \%chunk_paths);
	}
}

#***************************************************************************
# Subroutine:  derive_orfs
# Description: 
#***************************************************************************
sub derive_orfs {

	my ($self, $db_obj, $refseq_name, $chunk_paths_ref) = @_; 

	my $loci_table = $db_obj->{loci_table};
	my $orf_table  = $db_obj->{orf_table};
	
	my $blast_obj = $self->{blast_obj};
	my $threshold = $self->{minimum_orf_length};
	unless ($threshold)  { die; }
	my $blast_bin_path = $blast_obj->{blast_bin_path};
	
	# Get the extracted table reference names
	my @loci;
	my @fields = qw [ record_id chunk_name scaffold orientation
	                  assigned_to organism
					  extract_start extract_end genome_structure ];
	my $where = " WHERE assigned_to = '$refseq_name' 
	                AND genome_structure != 'LTR'
	                AND genome_structure != 'LEA'
	                AND genome_structure != 'LTR-LEA' ";
	$loci_table->select_rows(\@fields, \@loci, $where);
	my $i=0;
	foreach my $loci_ref (@loci) {
		
		# Get locus data
		$i++;
		my $locus_id      = $loci_ref->{record_id};
		my $chunk_name    = $loci_ref->{chunk_name};
		my $scaffold      = $loci_ref->{scaffold};
		my $orientation   = $loci_ref->{orientation};
		my $extract_start = $loci_ref->{extract_start};
		my $extract_end   = $loci_ref->{extract_end};
		my $structure     = $loci_ref->{genome_structure};
		my $assigned_to   = $loci_ref->{assigned_to};
		my $organism      = $loci_ref->{organism};
		print "\n\t Looking for ORFs in $refseq_name locus: id $locus_id ($structure)";

		# Get target path
		my $target_path = $chunk_paths_ref->{$chunk_name};
		unless ($target_path) { die "\n\t Path to target file not found\n\n"; }

		# Parsing for blastdbcmd
		my @gi = split(/\|/,$scaffold);	
		if (scalar(@gi) > 1) {
			$scaffold = $gi[1];
		}

		# Create the command
		# Command example: 
		# /bin/blast/blastdbcmd -db hs_alt_HuRef_chrX.fa -entry 157734237 
		# -range 10-60 -strand minus
		my $command = $blast_bin_path . "blastdbcmd -db $target_path";
		$command .= " -entry $scaffold ";
		$command .= " -range $extract_start-$extract_end ";
		if ($orientation eq '-ve') { $command .= ' -strand minus '; }
		
		# Execute the command
		my @sequence = `$command`;
		shift @sequence;  # Remove header
		my $sequence = join ('', @sequence);
		$sequence =~ s/\n//g;
		#print "\n\t SEQUENCE $sequence ";

		my $header = $scaffold . "_$extract_start-" . $extract_end;
		my $seq_obj = Sequence->new($sequence, $header, $locus_id);
		my %data;
		$seq_obj->get_putative_coding_sections(\%data);
		#if ($i eq 2) { die;	}
		#$devtools->print_hash(\%data);
		
		# Insert any long ORFs
		my @frames = 1..3;
		foreach my $frame (@frames) {
			
			my $frame_results = $data{$frame};
			my $orfs_ref      = $frame_results->{orfs};
			my @start_keys = sort by_number keys %$orfs_ref;
			foreach my $start (@start_keys) {
				
				my $orf_seq = $orfs_ref->{$start};
				my $seq_len = length $orf_seq;
				my $end = $start + ($seq_len * 3);
				
				if ($seq_len >= $threshold) {
					my %insert_data;
					$insert_data{locus_id}             = $locus_id;
					$insert_data{assigned_to}          = $assigned_to;
					$insert_data{organism}             = $organism;
					$insert_data{chunk_name}           = $chunk_name;
					$insert_data{scaffold}             = $scaffold;
					$insert_data{orf_extract_start}    = $start;
					$insert_data{orf_extract_end}      = $end;
					$insert_data{long_orf_gene_assign} = 'tmp';
					$insert_data{long_orf}             = $orf_seq;
					$insert_data{orf_length}           = $seq_len;
					$orf_table->insert_row(\%insert_data);
				}
			}
		}
	}
}

############################################################################
# CpG CORRECTION
############################################################################

#***************************************************************************
# Subroutine:  get_cpg_form_params 
# Description: initialise CpG program with parameters received via CGI 
#***************************************************************************
sub get_cpg_form_params {

	my ($self, $query, $data_ref) = @_;

	if ( $query->param('filename')) {
		$data_ref->{fasta}    = $query->upload('filename');
		$data_ref->{cpg_prop} = $query->param('cpg_prop');
	}
	# print the formatted input page if no sequence data has been received
	else {
		my @html;
		$fileio->read_input_file('site/html/reconstruct.html', \@html);
		my $html  = join('', @html);
		print $html;
		exit;
	}
}

#***************************************************************************
# Subroutine:  do_cpg_reconstruct 
# Description:  
#***************************************************************************
sub do_cpg_reconstruct {

	my ($self) = @_;
	# Loader get CGI params
	my $seq_obj      = Sequence->new();
	my $process_id   = $self->{process_id};

	# Initialise CGI output
	my $query = new CGI;
	my $loader_obj   = Loader->new();
	$loader_obj->get_cpg_form_params($query, $self);

	my $views_obj   = Views->new();
	my $site        = $views_obj->{paleo_top};
	my $header_path = $views_obj->{paleo_top}->{header_path};
	my $output_path = $self->{output_path};
	unless ($header_path and $output_path) { die; }
	
	# Create a unique directory to store the output of this process
	my $report_dir = $output_path . $process_id . '/';
	my $mkdir_cmd  = "mkdir $report_dir";
	my $result = system $mkdir_cmd;
	if ($result > 0) {
		#die "\n\t ### Error: couldn't create output directory - check permissions\n\n";	
		die "<br>Internal error, contact webmaster<br>";
	}
	$self->{output_path} = $report_dir;
	
	my $file = $self->{fasta};
	$file=~m/^.*(\\|\/)(.*)/; # strip the remote path and keep the filename 
	my $path = $report_dir . $file . '_local.fasta';

	open(LOCAL, ">$path") or die $!;
	while(<$file>) 
	{ print LOCAL $_; } 
	#print "$file has been successfully uploaded... thank you.\n";
	
	# Convert mac line breaks
	my $command = "perl -pi -e 's/\r/\n/g' $path";
	system $command;
	
	my @sequences;
	$seqio->read_fasta($path, \@sequences);
	#print $path;
	#$devtools->print_array(\@sequences); die;
	
	# Analyse alignment
	my $refseq = shift @sequences;
	#$devtools->print_hash($refseq); die;
	#my $rss = RefSeqAlignment->new($refseq, \@sequences);
	#$rss->cpg_reconstruct($refseq, \@sequences);

	#my $cpg_prop     = $self->{cpg_prop};
	#my $seqgen_obj   = $self->{seqgen_obj};
	#my @paleo_header;
	#$fileio->read_input_file($header_path, \@paleo_header);
	#my $paleo_header  = join('', @paleo_header);
	
}

############################################################################
# Classify envelopes
############################################################################

#***************************************************************************
# Subroutine:  classify_envelopes 
# Description: 
#***************************************************************************
sub classify_envelopes {

	my ($self, $ctl_file) = @_;

	my $db_obj = $self->{db_obj};
	my $loci_table = $db_obj->{loci_table};
	my $envelopes_table = $db_obj->{env_table};
	unless ($envelopes_table) {
		$devtools->print_hash($db_obj); die;
	}
	$envelopes_table->flush();
	$envelopes_table->reset_primary_keys();

	# Set up chunk paths
	my @chunk_keys;
	my @fields = qw [ organism chunk_name ]; 
	$loci_table->select_distinct(\@fields, \@chunk_keys);
	my %chunk_paths;
	$self->get_chunk_paths(\@chunk_keys, \%chunk_paths);
	#$devtools->print_hash(\%chunk_paths);	 die;

	my @loci;
	my @fields2 = qw [ record_id assigned_to extract_start extract_end organism
                      orientation chunk_name scaffold genome_structure ];

	$loci_table->select_rows(\@fields2, \@loci);
	foreach my $row (@loci) {
		
		my $record_id   = $row->{record_id};
		my $organism    = $row->{organism};
		my $assigned_to = $row->{assigned_to};
		my $chunk_name  = $row->{chunk_name};
		my $scaffold    = $row->{scaffold};
		my $structure   = $row->{genome_structure};
		my $start       = $row->{extract_start};
		my $end         = $row->{extract_end};
		my $orientation = $row->{orientation};
		unless ( $record_id and $organism and $assigned_to
                 and $start and $end and $orientation ) {
			$devtools->print_hash($row); die;
		}
		if ($structure =~ /^.-LTR-.$/) {next;}
		if ($structure =~ /-T$/) {
			if($structure !~ /LTR-T$/) {next;}
		}
		# Get path to chunk
		my $key = $organism . '|' . $chunk_name;
		#print "\n\t$key\n";
		
		# Set the coordinates to get any flanking downstream from gag/pol
		if($structure !~ /LTR-.$/){
			if(($structure =~ /pol/) || ($structure =~ /env/)){
				if($orientation =~ /-ve/){
					$start-=2000;
				}else{
					$end+=2000;
				}
			}else{
				if($orientation =~ /-ve/){
                    $start-=7000;
                }else{
                    $end+=7000;
                } 	
			}	
		}
		if($start <= 0){
			$start = 1;
		}
		my @header = ( $record_id, $organism, $assigned_to, $start, $end, $orientation );
		my $header = join ('_', @header);
        print "\n\t Extracting $record_id from $organism";
        print "\n\t # $chunk_name";
        print "\n\t # $scaffold, $start-$end";
        print "\n\t # $structure";
		my $path = $chunk_paths{$key}->{path};
		unless ($path)  { die; }

		my $fasta  = $self->blast_extract($header, $path, $scaffold, $start, $end, $orientation);
		my $temp_fasta = './site/tmp/tmp_loci_4env.fasta';
		my $env_blast_ref;
		$fileio->write_text_to_file($temp_fasta, $fasta);		
		$env_blast_ref = $self->get_env_type($temp_fasta);
		#$devtools->print_hash($env_blast_ref); #  die;
		unless ($env_blast_ref->{scaffold})  { print "\n\tNo ENV\n"; next;}
		unless ($env_blast_ref->{bit_score} > $env_threshold)  { print "\n\tNo ENV\n"; next;}
		
		my %env_data;
		$env_data{locus_id}		= $record_id;
		$env_data{assigned_to}	= $assigned_to;
		$env_data{env_assign}	= $env_blast_ref->{scaffold};
		$env_data{genome_structure}	= $structure;
		$env_data{orientation}	= $orientation;
		$env_data{bit_score}	= $env_blast_ref->{bit_score};
		$env_data{identity}		= $env_blast_ref->{identity};
		$env_data{e_value_num}	= $env_blast_ref->{e_value_num};
		$env_data{e_value_exp}	= $env_blast_ref->{e_value_exp};
		$env_data{subject_start}= $env_blast_ref->{aln_start};
		$env_data{subject_end}	= $env_blast_ref->{aln_stop};
		if($orientation =~ /-ve/){
			$env_data{query_start}	= $end - $env_blast_ref->{query_stop} + 1;
			$env_data{query_end}	= $end - $env_blast_ref->{query_start} + 1;
		}else{
			$env_data{query_start}  = $start + $env_blast_ref->{query_start} - 1;
			$env_data{query_end}    = $start + $env_blast_ref->{query_stop} - 1;
		}
		$env_data{align_len}	= $env_blast_ref->{align_len};
		#print "gap: $env_blast_ref->{gap_openings}\n";
		$env_data{gap_openings} = "'NULL'";
		if($env_blast_ref->{gap_openings} == 0){
			$env_data{gap_openings} = "'None'";
		}
		$env_data{mismatches}	= $env_blast_ref->{mismatches};
		#$devtools->print_hash(\%env_data);   die;
		#my $referencia = \%env_data;
		#print "$referencia\n";

		$envelopes_table->insert_row(\%env_data);
	}
}

#***************************************************************************
# Subroutine:   get_env_type
# Description:  
#***************************************************************************
sub get_env_type {

    my ($self, $fasta_path) = @_;
    
    my $result_file = $fasta_path . '.blast_result';
     
    # BLASTx against the env library
    my $lib_path = './db/blast/refseq_envs_AA.fasta';
    my $blast_alg = 'blastx';
    #print "\n\t ##### ENV BLAST for $name";
    my $blast_obj = $self->{blast_obj};
    unless ($blast_obj) { die; }
    $blast_obj->blast($blast_alg, $lib_path, $fasta_path, $result_file);
    my @results;
    $blast_obj->parse_tab_format_step_one($result_file, \@results);
    # Get the best match from this file
    #my $top_match_ref = shift @results;
	foreach my $match (@results){
		if($match->{orientation} eq '-ve'){
			next;
		}else{
			return $match;
			last;
		}
	}
	#return $top_match_ref;
	#$devtools->print_hash($top_match_ref);   die;
}

############################################################################
# SELECT FXNS
############################################################################

#***************************************************************************
# Subroutine:  get_leader_loci
# Description: 
#***************************************************************************
sub get_leader_loci {
	
	my ($self, $refseq,  $chunk_paths_ref, $data_ref) = @_;

	my $blast_bin_path = $self->{blast_obj}->{blast_bin_path};

	# Get all data from loci table
	my $db_name = $self->{db_name};
	my $db_obj_tmp = DB->new();
	$db_obj_tmp->load_screening_db($db_name);
	my $loci_table = $db_obj_tmp->{loci_table};
	my $refseq_name = $refseq->{name};
	my @loci;
	my @fields = qw [ record_id assigned_to extract_start extract_end 
					  organism orientation chunk_name scaffold genome_structure ];
	my $where  = " WHERE Assigned_to = '$refseq_name' ";
	   $where .= ' AND Genome_structure LIKE \'%LTR-LEA%\'';
	$loci_table->select_rows(\@fields, \@loci, $where);
	foreach my $locus_ref (@loci) {
	
		#$devtools->print_hash($locus_ref);
		my $record_id   = $locus_ref->{record_id};
		my $start       = $locus_ref->{extract_start};
		my $end         = $locus_ref->{extract_end};
		my $orientation = $locus_ref->{orientation};
		my $scaffold    = $locus_ref->{scaffold};
		my $chunk_name  = $locus_ref->{chunk_name};
		my $chunk_path = $chunk_paths_ref->{$chunk_name};

		# Parsing for blastdbcmd
		my @gi = split(/\|/,$scaffold);	
		if (scalar(@gi) > 1) {
			$scaffold = $gi[1];
		}

		# Create the command
		# Command example: 
		# /bin/blast/blastdbcmd -db hs_alt_HuRef_chrX.fa -entry 157734237 
		# -range 10-60 -strand minus
		my $command = $blast_bin_path . "blastdbcmd -db $chunk_path";
		$command .= " -entry $scaffold ";
		$command .= " -range $start-$end ";
		if ($orientation eq '-ve') { $command .= ' -strand minus '; }
		
		# Execute the command
		my @sequence = `$command`;
		shift @sequence;  # Remove header
		my $sequence = join ('', @sequence);
		$sequence =~ s/\n//g;
			
		# Write to file
		my %data;
		$data{locus_id} = $record_id;
		$data{sequence} = $sequence;
		$data{assigned_to} = $locus_ref->{assigned_to};
		push(@$data_ref, \%data);
	}
}

#***************************************************************************
# Subroutine:  select_from_loci_table
# Description: 
#***************************************************************************
sub select_from_loci_table {

	my ($self, $loci_ref, $where) = @_;

	# Get settings from the ctl file
	my $db_obj     = $self->{db_obj};
	my $chunks_ref = $self->{chunk_data};
	my $loci_table = $db_obj->{loci_table};
	my @fields = qw [ record_id assigned_to extract_start extract_end organism
                      orientation chunk_name scaffold genome_structure ];
	$loci_table->select_rows(\@fields, $loci_ref, $where);
	#$devtools->print_array(\@loci); die;
	
	my @fasta;
	foreach my $row (@$loci_ref) {
	
		my $record_id   = $row->{record_id};
		my $organism    = $row->{organism};
		my $assigned_to = $row->{assigned_to};
		my $chunk_name  = $row->{chunk_name};
		my $scaffold    = $row->{scaffold};
		my $structure   = $row->{genome_structure};
		my $start       = $row->{extract_start};
		my $end         = $row->{extract_end};
		my $orientation = $row->{orientation};
		unless ( $record_id and $organism and $assigned_to and $start and $end and $orientation ) {
			$devtools->print_hash($row); die;
		}
		print "\n\t Extracting $record_id from $organism";
		print "\n\t # $chunk_name";
		print "\n\t # $scaffold, $start-$end";
		
		# Get path to chunk
		my $key = $organism . '_' . $chunk_name;
		my $data_ref = $chunks_ref->{$key};
		my $path     = $data_ref->{path};
		
		unless ($path)  { die; }
		my @header = ( $record_id, $organism, $assigned_to, $start, $end, $orientation );
		my $header = join ('_', @header);
		my $fasta  = $self->blast_extract($header, $path, $scaffold, $start, $end, $orientation);
		push (@fasta, $fasta);

	}
	$devtools->print_array(\@fasta); die;
}

############################################################################
# General Paleotools functions 
############################################################################

#***************************************************************************
# Subroutine:  get_chunk_paths 
# Description: 
#***************************************************************************
sub get_chunk_paths {

	my ($self, $chunks_ref, $chunk_paths) = @_;

	my $db_obj = $self->{db_obj};
	$db_obj->load_genome_database();
	my $chunks_table = $db_obj->{genome_chunks_table};
	foreach my $chunk_ref (@$chunks_ref) {

		#$devtools->print_array($chunk_ref);	die;
		#$devtools->print_hash($chunk_ref);	die;
		my $organism   = $chunk_ref->{organism};
		my $chunk_name = $chunk_ref->{chunk_name};
		my @fields = qw [ source_type version grouping ];
		my @data;
		my $where  = " WHERE Organism = '$organism' and chunk_name = '$chunk_name'";
		$chunks_table->select_rows(\@fields, \@data, $where);
		my $num_rows = scalar @data;
		if ($num_rows > 1)  { die "\n\t Duplicate entries in genome DB\n\n"; }
		my $row_ref = shift @data;
		my $source_type = $row_ref->{source_type};
		my $version     = $row_ref->{version};
		my $grouping    = $row_ref->{grouping};
		
		my $path =  $genome_path; 
           $path .= "$grouping/";
		   $path .= "$organism/";
		   $path .= "$source_type/";
		   $path .= "$version/";
		   $path .= "$chunk_name";
		
		# Store the data
		$chunk_ref->{grouping}    = $grouping;
		$chunk_ref->{source_type} = $source_type;
		$chunk_ref->{version}     = $version;
		$chunk_ref->{path}        = $path;
		
		my $key = $organism . '|' . $chunk_name;
		$chunk_paths->{$key} = $chunk_ref;
	}
	#$devtools->print_hash($chunk_paths);	 die;
}

#***************************************************************************
# Subroutine:  blast_extract
# Description: 
#***************************************************************************
sub blast_extract {

	my ($self, $header, $chr_path, $scaf, $start, $end, $orientation) = @_;
	#print "\n\t$chr_path\n"; die;
	unless ($header and $chr_path and $scaf and $start and $end and $orientation) { die; }

	# Do neccessary adjustments on contig name	
	my $truescaf='';
	my @gi = split(/\|/,$scaf);
	if (scalar(@gi) > 1)  { $truescaf = $gi[1]; }
    else                  { $truescaf = $scaf;  }

	# Call the BLAST extract process
	my $seq = '';
	my @subseq;
	
	my $scaffold_len = `./bin/blast/blastdbcmd -db $chr_path -entry $truescaf -outfmt %l`;

	if($end > $scaffold_len){
		$end = $scaffold_len;
	}	


	if ($orientation eq '-ve'){
		@subseq = `./bin/blast/blastdbcmd -db $chr_path -entry $truescaf -range $start-$end -strand minus`;
	}
	else {
		@subseq = `./bin/blast/blastdbcmd -db $chr_path -entry $truescaf -range $start-$end`;
	}
	#my $eader = 'test';

	# Remove the BLAST header and concatenate 
	shift @subseq;	
	foreach my $line (@subseq){
		$seq .= $line;
	}
	my $fasta = ">$header\n$seq\n";
	return $fasta;
}

#***************************************************************************
# Subroutine:  show_title
# Description: does what it says 
#***************************************************************************
sub show_title {

	$console->refresh();
	
	my $title       = 'VGLUE Excavation tools';
	my $version     = '2.0';
	my $description = 'Sequence Database Mining using the GLUE framework';
	my $author      = 'Robert J. Gifford';
	my $contact		= '<rgifford@adarc.org>';
	
	$console->show_about_box($title, $version, $description, $author, $contact);
}

############################################################################
# BASIC FXNS
############################################################################

sub by_number { $a <=> $b }	

############################################################################
# EOF
############################################################################
