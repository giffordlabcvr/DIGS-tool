#!usr/bin/perl -w
############################################################################
# Module:      Utility.pm   
# Description: DIGS tool utility functions
# History:     December  2017: Created by Robert Gifford 
############################################################################
package Utility;

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
use Base::DevTools;

# Program components
use DIGS::ScreenBuilder; # Functions to set up screen

############################################################################
# Globals
############################################################################

# Base objects
my $fileio    = FileIO->new();
my $console   = Console->new();
my $devtools  = DevTools->new();
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create new Utility 'object'
#***************************************************************************
sub new {

	my ($invocant, $digs_obj) = @_;
	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {

		# Global settings
		process_id             => $digs_obj->{process_id},
		program_version        => $digs_obj->{program_version},
		
		# DIGS tool object
		digs_obj => $digs_obj,

	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# UTILITY FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  show_utility_help_page
# Description: show help page information for utility functions
#***************************************************************************
sub show_utility_help_page {

	my ($self) = @_;

	# Create utility help menu
	my $program_version = $self->{program_version};

    my $HELP   = "\n\t ### DIGS version $program_version";
       $HELP .= "\n\t ### usage: $0 m=[option] -i=[control file] -e=[utility help]\n";

       $HELP  .= "\n\t ### Utility functions\n"; 
	   $HELP  .= "\n\t -u=1  Add extra tables to screening DB"; 
	   $HELP  .= "\n\t -u=2  Flush screening DB"; 
	   $HELP  .= "\n\t -u=3  Drop screening DB"; 
	   $HELP  .= "\n\t -u=4  Show BLAST chains"; 
	   $HELP  .= "\n\t -u=5  Show locus chains"; 
	   $HELP  .= "\n\t -u=6  Show nomenclature chains"; 
	   $HELP  .= "\n\t -u=7  Summarise genomes (short, by species)";
	   $HELP  .= "\n\t -u=8  Summarise genomes (long, by target file)";
	   $HELP  .= "\n\t -u=9  Translate DB schema"; 
	   $HELP  .= "\n\t -u=10 Create standard locus IDs\n"; 
	   $HELP  .= "\n\t -u=11 Upload data to digs_results table\n"; 
	   #$HELP  .= "\n\t -u=12  Extract sequences using track";
	   $HELP  .= "\n\n"; 

	print $HELP;
}

#***************************************************************************
# Subroutine:  run_utility_process
# Description: handler for DIGS tool utility functions 
#***************************************************************************
sub run_utility_process {

	my ($self, $option, $infile) = @_;

 	# Show title
	my $digs_obj = $self->{digs_obj};
	$digs_obj->show_title();  

	# If we're doing anything else we need a control file as input
	unless ($infile) { die "\n\t Option '$option' requires an infile\n\n"; }

	# Try opening control file
	my @ctl_file;
	my $valid = $fileio->read_file($infile, \@ctl_file);
	unless ($valid) {  # Exit if we can't open the file
		die "\n\t ### Couldn't open control file '$infile'\n\n\n ";
	}
	
	# If control file looks OK, store the path and parse the file
	$self->{ctl_file} = $infile;
	my $loader_obj = ScreenBuilder->new($digs_obj);
	$loader_obj->parse_control_file($infile, $digs_obj, $option);

	# Store the ScreenBuilder object (used later)
	$self->{loader_obj} = $loader_obj;

	# Create the output directories
	$loader_obj->create_output_directories($self);

	# Load/create the screening database
	my $db_name = $loader_obj->{db_name};
	unless ($db_name) { die "\n\t Error: no DB name defined \n\n\n"; }
	$digs_obj->initialise_screening_db($db_name);

	# Hand off to functions 
	if ($option eq 1) { # Add a table of data to the screening database
		$self->extend_screening_db();
	}
	elsif ($option eq 2) { # Flush screening DB
		my $db = $digs_obj->{db};
		my $db_name = $db->{db_name};
		my $question = "\n\n\t  Are you sure you want to flush data in the $db_name database?";
		my $answer1 = $console->ask_yes_no_question($question); # Ask to make sure
		if ($answer1 eq 'y') { $db->flush_screening_db(); }
	}
	elsif ($option eq 3) { # Drop screening DB 
		my $db = $digs_obj->{db};
		$db->drop_screening_db();    
	}
	elsif ($option eq 4) { # Show the BLAST chains for each extracted locus
		$self->show_blast_chains();
	}
	elsif ($option eq 5) { # Consolidate DIGS results into higher level loci 
		$self->show_locus_chains();
	}
	elsif ($option eq 6) { # Consolidate DIGS results into higher level loci 
		$self->show_nomenclature_chains();
	}
	elsif ($option eq 7)    { # Summarise target genome directory (short)
		my $target_db_obj = TargetDB->new($digs_obj);
		$target_db_obj->summarise_targets_short();
	}
	elsif ($option eq 8) { # Summarise target genome directory (long)
		my $target_db_obj = TargetDB->new($digs_obj);
		$target_db_obj->summarise_targets_long();
	}
	elsif ($option eq 9) { # DB schema translation
		die;
		my $db_obj = $self->{db};
		$db_obj->translate_schema();
	}
	elsif ($option eq 10) { # Standardised locus naming
		die;
		$self->upload_data_to_digs_results_table();
	}
	elsif ($option eq 10) { # Standardised locus naming
		die;
		$self->create_standard_locus_ids();
	}
	elsif ($option eq 12) {
		die;
		unless ($infile) {  die "\n\t Option '$option' requires an infile\n\n"; }
		my $loader_obj = ScreenBuilder->new($self);
		my @extracted;
		$self->extract_track_sequences(\@extracted, $infile);
	}
	else {
		print "\n\t  Unrecognized option '-u=$option'\n";
	}
}

#***************************************************************************
# Subroutine:  extend_screening_db
# Description: console managemant of ancillary tables in the screening database
#***************************************************************************
sub extend_screening_db { 

}

#***************************************************************************
# Subroutine:  show_blast_chains
# Description: Show BLAST chains for all extracted sequences
#***************************************************************************
sub show_blast_chains {
	
	my ($self) = @_;

	# Get relevant variables and objects
	my $digs_obj = $self->{digs_obj};
	my $db = $digs_obj->{db};
	unless ($db) { die; } # Sanity checking
	my $digs_results_table = $db->{digs_results_table}; 
	my $blast_chains_table = $db->{blast_chains_table};
	my $extract_where = " ORDER BY record_id ";
	my @extracted_ids;
	my @fields = qw [ record_id assigned_name assigned_gene ];
	$digs_results_table->select_rows(\@fields, \@extracted_ids, $extract_where);	 

	# Iterate through the digs result rows	
	foreach my $hit_ref (@extracted_ids) {
		my $digs_result_id = $hit_ref->{record_id};
		my @chain;
		my $assigned_name = $hit_ref->{assigned_name};
		my $assigned_gene = $hit_ref->{assigned_gene};
		my @chain_fields = qw [ record_id probe_name probe_gene 
		                        organism target_name 
		                        scaffold subject_start subject_end
		                        bitscore identity align_len ];
		my $blast_where  = " WHERE digs_result_id = $digs_result_id ";
		   $blast_where .= " ORDER BY subject_start";
		$blast_chains_table->select_rows(\@chain_fields, \@chain, $blast_where);	
		print "\n\t ### BLAST hit chain for extracted locus $digs_result_id";
		print " ($assigned_name, $assigned_gene):";
		foreach my $hit_ref (@chain) {
			my $blast_id    = $hit_ref->{record_id};
			my $probe_name  = $hit_ref->{probe_name};
			my $probe_gene  = $hit_ref->{probe_gene};
			my $organism    = $hit_ref->{organism};
			my $scaffold    = $hit_ref->{scaffold};
			my $start       = $hit_ref->{subject_start};
			my $end         = $hit_ref->{subject_end};
			my $bitscore    = $hit_ref->{bitscore};
			my $identity    = $hit_ref->{identity};
			my $f_identity  = sprintf("%.2f", $identity);
			my $align_len   = $hit_ref->{align_len};
			print "\n\t\t $blast_id:\t Score: $bitscore, \%$f_identity identity ";
			print "across $align_len aa ($start-$end) to:\t $probe_name ($probe_gene) ";
		}
	}
}

#***************************************************************************
# Subroutine:  show_locus_chains
# Description: Show composition of consolidated loci
#***************************************************************************
sub show_locus_chains {
	
	my ($self) = @_;

	# Get relevant variables and objects
	my $digs_obj = $self->{digs_obj};
	my $db_ref = $digs_obj->{db};
	unless ($db_ref) { die; } # Sanity checking
	my $dbh = $db_ref->{dbh};

	my $loci_exists = $db_ref->does_table_exist('loci');
	unless ($loci_exists) {
		print "\n\t  The locus tables don't seem to exist (have you ran consolidate?\n\n";
		exit;
	}

	$db_ref->load_loci_table($dbh);
	$db_ref->load_loci_chains_table($dbh);
	#$devtools->print_hash($db); die;

	my $digs_results_table = $db_ref->{digs_results_table}; 
	my $loci_table         = $db_ref->{loci_table};
	my $loci_chains_table  = $db_ref->{loci_chains_table};
	unless ($digs_results_table and $loci_chains_table) {
		print "\n\t # Locus tables not found - run consolidate first\n\n\n";
		exit;
	}

	# Get all loci
	my $loci_where = " ORDER BY record_id ";
	my @loci;
	my @fields = qw [ record_id locus_structure ];
	$loci_table->select_rows(\@fields, \@loci, $loci_where);	 
	
	# Iterate through loci
	foreach my $locus_ref (@loci) {

		my $locus_id = $locus_ref->{record_id};
		print "\n\t ### Chain $locus_id ";	
		my $chain_where = " WHERE locus_id = $locus_id ";
		my @results;
		my @fields = qw [ record_id locus_id digs_result_id ];
		$loci_chains_table->select_rows(\@fields, \@results, $chain_where);	 
		
		foreach my $result_ref (@results) {

			my @digs_results;
			my $digs_result_id = $result_ref->{digs_result_id};
			my @result_fields = qw [ assigned_name assigned_gene 
			                         organism target_name 
			                         scaffold extract_start extract_end
			                         bitscore identity align_len ];
			my $where  = " WHERE record_id = $digs_result_id ";
			$digs_results_table->select_rows(\@result_fields, \@digs_results, $where);			

			foreach my $result_ref (@digs_results) {
		
				my $assigned_name = $result_ref->{assigned_name};
				my $assigned_gene = $result_ref->{assigned_gene};
				my $organism      = $result_ref->{organism};	
				my $scaffold      = $result_ref->{scaffold};
				my $start         = $result_ref->{extract_start};
				my $end           = $result_ref->{extract_end};
				my $bitscore      = $result_ref->{bitscore};
				my $identity      = $result_ref->{identity};
				my $align_len     = $result_ref->{align_len};
				my $f_identity    = sprintf("%.2f", $identity);
				print "\n\t\t $digs_result_id:\t Score: $bitscore, \%$f_identity identity ";
				print "across $align_len aa ($start-$end) to:\t $assigned_name ($assigned_gene) ";
	
			}
		}
	}
}

#***************************************************************************
# Subroutine:  show_nomenclature_chains
# Description: Show nomenclature chains for all consolidated annotations
#***************************************************************************
sub show_nomenclature_chains {
	
	my ($self) = @_;

	# Get relevant variables and objects
	my $digs_obj = $self->{digs_obj};
	my $db = $digs_obj->{db};
	unless ($db) { die; } # Sanity checking

	# Get relevant variables and objects
	my $dbh = $db->{dbh};
	$db->load_nomenclature_table($dbh);
	$db->load_nomenclature_chains_table($dbh);
	my $nomenclature_table = $db->{nomenclature_table}; 
	my $chains_table  = $db->{nomenclature_chains_table};
	unless ($nomenclature_table and $chains_table) { die; } # Sanity checking

	# Get all named loci
	my $nom_where = " ORDER BY record_id ";
	my @loci;
	my @fields = qw [ record_id full_id ];
	$nomenclature_table->select_rows(\@fields, \@loci, $nom_where);	 
	
	# Iterate through loci
	foreach my $locus_ref (@loci) {

		my $locus_id = $locus_ref->{record_id};
		print "\n\t ### Chain $locus_id: ";	
		my $chain_where = " WHERE nomenclature_locus_id = $locus_id ";
		my @results;
		@fields = qw [ record_id track_id nomenclature_locus_id ];
		$chains_table->select_rows(\@fields, \@results, $chain_where);	 
		
		foreach my $result_ref (@results) {

			my @digs_results;
			my $track_entry_id = $result_ref->{track_id};
			my @fields = qw [ assigned_name assigned_gene ];
			my $where  = " WHERE record_id = $track_entry_id ";
			$nomenclature_table->select_rows(\@fields, \@digs_results, $where);
			
			foreach my $result_ref (@digs_results) {
		
				my $assigned_name = $result_ref->{assigned_name};
				my $assigned_gene = $result_ref->{assigned_gene};
				print " $track_entry_id ";
	
			}
		}
	}
}

#***************************************************************************
# Subroutine:  extract_track_sequences
# Description: extract FASTA nucs from a genome assembly using an input track 
#***************************************************************************
sub extract_track_sequences {
	
	my ($self, $extracted_ref, $track_path) = @_;

	# Get paths, objects, data structures and variables from self
	my $blast_obj   = $self->{blast_obj};

	# Try to read the tab-delimited infile
	print "\n\n\t #### WARNING: This function expects a tab-delimited data table with column headers in order!";
	my $question1 = "\n\n\t Please enter the path to the file with the table data and column headings\n\n\t";
	my $query_na_track_genome_path = $console->ask_question($question1);
	
	# Read FASTA probe library
	my @track;
	$fileio->read_file($track_path, \@track);

	# Iterate through the tracks extracting
	my $i = 0;	
	my @probe_fasta;
	foreach my $line (@track) {
		
		$i++;
		
		chomp $line; # remove newline
		my @line = split("\t", $line);

		my $j;
		#foreach my $element (@line) {
		#	$j++;
		#	print "\n\t ELEMENT $j : $element"
		#}
		#die;
		
		my $name          = $line[0];
		my $scaffold      = $line[2];
		
		# TODO - hacky - resolve
		if ($scaffold =~ /Unk/) { next; }
		
		my $subject_start = $line[3];
		my $subject_end   = $line[4];
		my $gene          = $line[5];
		my $id            = $line[0];
		my $orientation;

		if ($subject_start < $subject_end) {
			$orientation = '+';
		}
		elsif ($subject_start > $subject_end) {
			$orientation = '-';
			my $switch     = $subject_start;
			$subject_start = $subject_end;
			$subject_end   = $switch;
		}
		
		# Extract the sequence
		my %data;
		$data{subject_start} = $subject_start;
		$data{subject_end}   = $subject_end;
		$data{orientation}   = $orientation;
		$data{scaffold}      = $scaffold;

		my $target_path = "$query_na_track_genome_path" . "/$scaffold" . '.fa';;	
		my $sequence = $blast_obj->extract_sequence($target_path, \%data);
		if ($sequence) {	
			#print "\n\t got seq $sequence \n\n";
			my %probe;
			$probe{probe_name}      = $name;
			$probe{probe_gene}      = $gene;
			$probe{probe_id}        = $name . "_$gene";
			$probe{sequence}        = $sequence;
			$probe{probe_type}      = 'UTR';
			$probe{blast_alg}       = 'blastn';
			
			push(@$extracted_ref, \%probe);
			
			my $header = "$name" . "_$gene";
			$header =~ s/\(/\./g;
			$header =~ s/\)//g;

			print "\n\t\t Getting probe $i: $header";

			my $digs_fasta = ">$header" . "\n$sequence\n";
			push (@probe_fasta, $digs_fasta);;
		
		}
		else {
			die "\n\t Sequence extraction failed";
		}
	}	
	
	my $outfile = 'extracted.DIGS.fna';
	$fileio->write_file($outfile, \@probe_fasta);
}

############################################################################
# EOF
############################################################################