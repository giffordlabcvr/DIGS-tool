#!/usr/bin/perl -w
############################################################################
# Script:      digs_for_eves_tools.pl
# Creator:     R.J. Gifford
# Description: DIGS for EVEs perl tools
# History:     Version 1.0
############################################################################

# Include the PERL module library for DIGS 
use lib ($ENV{DIGS_HOME} . '/modules/'); 

# Other paths
my $genome_use_path = $ENV{DIGS_GENOMES} . '/'; 
$genome_use_path =~ s/\/\//\//g; # Remove any double backslashes
my $blast_bin_path  = '';  # left empty if BLAST+ programs are in your path 


############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;
use Getopt::Long;
use DBI;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base modules
use Base::Console;
use Base::FileIO;
use Base::DevTools;

# DIGS modules
use DIGS::DIGS;
use DIGS::Initialise;


# Third party program interface modules
use Interface::MySQLtable;   # Interface to BLAST 

############################################################################
# Paths & Globals
############################################################################

# Create a unique process ID for this DIGS screening process
my $pid  = $$;
my $time = time;
my $process_id  = $pid . '_' . $time;
my $version = '1.0';

# DIGS for EVEs screening DB ancillary table names
my $eve_data_table_name = 'eve_data';
my $host_data_table_name = 'host_data';

# Flags
my $force;
my $verbose;
my $exclude_phage = 1;

# TAxonomy database connection globals
my $db_name  = 'NCBI_Taxonomy';
my $server   = 'localhost';

# Settings for reports
my @bitscore = qw [ 60 80 100 120 ];
my @groups =   qw [ vertebrates ];
my @families = qw [ Circoviridae Parvoviridae Hepadnaviridae Bornaviridae Filoviridae Chuviridae ];
	
############################################################################
# Instantiations
############################################################################

# Base utilites
my $fileio     = FileIO->new();
my $console    = Console->new();
my $devtools   = DevTools->new();

my $mysql_username = ($ENV{DIGS_MYSQL_USER}); 
my $mysql_password = ($ENV{DIGS_MYSQL_PASSWORD}); 

# Instantiate main program classes using global settings
my %params;
$params{program_version} = 'NULL';
$params{process_id}      = $process_id;
$params{blast_bin_path}  = $blast_bin_path; 
$params{genome_use_path} = $genome_use_path;
$params{mysql_username}  = $mysql_username ; 
$params{mysql_password}  = $mysql_password; 

my $digs_tool_obj = DIGS->new(\%params);

############################################################################
# Set up USAGE statement
############################################################################

# Initialise usage statement to print if usage is incorrect
my ($USAGE) = "\n\t  usage: $0 m=[option] -h=[help]\n\n";

############################################################################
# Main program
############################################################################

# Run script
main();

# Exit script
print "\n\n\t # Exit\n\n";
exit;

############################################################################
# Subroutines
############################################################################

#***************************************************************************
# Subroutine:  main
# Description: top level handler fxn
#***************************************************************************
sub main {
	
	# Options that require a file path
	my $infile   = undef;
	
	# Options that require a numerical value
	my $mode     = undef;
	my $help     = undef;

	# Read in options using GetOpt::Long
	GetOptions ('mode|m=i'       => \$mode, 
			    'help'           => \$help,
			    'verbose'        => \$verbose,
			    'force'          => \$force,
			    'infile|i=s'     => \$infile,
	) or die $USAGE;

	show_title();
	
	# Set flags based on options received
	if ($verbose) {  $digs_tool_obj->{verbose} = 'true'; }
	if ($force)   {  $digs_tool_obj->{force} = 'true';   }
	if ($mode) { 
		if ($mode eq 1) {
			format_ncbi_fasta_for_digs();
		}
		if ($mode eq 2) {
			format_ncbi_fasta_for_digs2($infile);
		}
		if ($mode eq 3) {
			get_taxonomy_from_name_list($infile);
		}
		if ($mode eq 4) {
			create_eve_report($infile);
		}
		if ($mode eq 5) {
			crunch_loci($infile);
		}
	}
	elsif ($help) { # Show help page
		show_help_page();
		exit;
	}
	else {
		die $USAGE;
	}
	
}

#***************************************************************************
# Subroutine:  create_eve_report
# Description:  
#***************************************************************************
sub create_eve_report {

	my ($ctl_file) = @_;

	unless ($ctl_file) {
		die "\n\t  This option requires a control file as input\n\n"; 		
	}

	# Load screening DB
	my $initialise_obj = Initialise->new();
	$initialise_obj->do_general_setup($digs_tool_obj, 1, $ctl_file);
	my $db_name = $digs_tool_obj->{db_name};
	unless ($db_name) {
		die "\n\t  Error - failed to get DB name from \n\n"; 		
	}

	# Initialise screening DB
	$initialise_obj->initialise_screening_db($digs_tool_obj, $db_name);

	# Get handle for the 'digs_results' table
	my $db_ref = $digs_tool_obj->{db};
	#$devtools->print_hash($db_ref); die;

	# Check consistency
	print   "\n\n\t  Checking that all viral references have taxonomy information";
	sleep 1;
	check_eve_data_table_matching($db_ref);
	
	# Create summary report
	print   "\n\n\t  Initiating report....";
	sleep 1;
	create_eve_summary_report($db_ref);
}

#***************************************************************************
# Subroutine:  create_eve_summary_report
# Description:  
#***************************************************************************
sub create_eve_summary_report {

	my ($db_ref) = @_;

	my %family_reports;
	create_virus_family_reports(\%family_reports, $db_ref);
	#$devtools->print_hash(\%family_reports); die;

	my @labels;
	foreach my $family (@families) {

		my $trim = $family;
		$trim =~ s/viridae//;
		push (@labels, "'$trim'");
	}
	
	foreach my $group (@groups) {

		# Create summary
		print   "\n\n\t  Creating summary report for '$group'";	

		my @combined_series;
		foreach my $bitscore (@bitscore) {

			#$devtools->print_hash($all_gene_data_ref); die;
			#$devtools->print_hash($family_data_ref); die;
			
			my @series;
			foreach my $family (@families) {

				my $family_data_ref = $family_reports{$family};	
				my $all_gene_data_ref = $family_data_ref->{'all'};	

				my $count_ref = $all_gene_data_ref->{$bitscore}; 
				my $count = $count_ref->{count}; 
				push (@series, "'$count'");
			}
 			my $series_nums = join(', ', @series);
 			my $series = "\n\t  [$series_nums],";
			push (@combined_series, $series);
 			 
		 }

 		my $labels = join(', ', @labels);
		my @final;
		push (@final, "\n var data = {");
		push (@final, "\n  labels: [ $labels ],");
		push (@final, "\n  series: [");

		push (@final, @combined_series);
		push (@final, "\n  ]");
		push (@final, "\n };");
		print @final;

	}

#        var data = {
#          labels: ['Borna', 'Filo', 'Nya', 'Hepadna', 'Circo', 'Parvo', 'Gemini', 'Ophio', 'Reo', 'Flavi', 'Virga', 'Betaflex'],
#          series: [
#            [542, 443, 320, 780, 553, 453, 326, 434, 568, 610, 756, 895],
#            [412, 243, 280, 580, 453, 353, 300, 364, 368, 410, 636, 695],
#            [41, 24, 300, 80, 253, 153, 100, 364, 268, 210, 336, 295],
#            [10, 20, 40, 60, 70, 80, 80, 80, 90, 90, 120, 129],
#         ]
#        };

}


#***************************************************************************
# Subroutine:  write_virus_family_reports
# Description:  
#***************************************************************************
sub create_virus_family_reports {

	my ($family_reports_ref, $db_ref) = @_;

	# Iterate through virus families to create family level reports
	foreach my $family (@families) {

		my %family_report;
		get_virus_family_report_data($db_ref, $family, \%family_report);
 		#$devtools->print_hash(\%family_report); die;
		
		# store this family rpeort 
		$family_reports_ref->{$family} = \%family_report;

		# Write the report
		write_virus_family_report(\%family_report, $family);
	}

}

#***************************************************************************
# Subroutine:  get_virus_family_report_data
# Description:  
#***************************************************************************
sub get_virus_family_report_data {

	my ($db_ref, $family, $report_ref) = @_;

	# Get the genes
	my $dbh = $db_ref->{dbh};
	my @genes;
	get_genes_by_virus_group($db_ref, 'virus_family', $family, \@genes);
	push(@genes, 'all');

	# Get reports by gene
	get_data_by_gene($db_ref, $report_ref, $family, \@genes);

}

#***************************************************************************
# Subroutine:  get_data_by_gene
# Description:  
#***************************************************************************
sub get_data_by_gene {

	my ($db_ref, $report_ref, $virus_family, $genes_ref) = @_;

	my $dbh = $db_ref->{dbh};

	# Iterate through genes
	foreach my  $gene (@$genes_ref) {
	
		my %gene_data;		
		print   "\n\n\t    Getting all hits to gene '$gene'";
		foreach my $bitscore (@bitscore) {
			
			my %bitscore_gene_data;
			
			print "\n\t Creating report for '$virus_family:$gene' at bitscore '$bitscore'";

			# Summarise the results for this family
			my $query1 = "SELECT assigned_name, assigned_gene, virus_genus, sequence_length, bitscore";
			  $query1 .= " FROM digs_results, $eve_data_table_name";
			  $query1 .= " WHERE digs_results.assigned_name = $eve_data_table_name.virus_species ";
			  $query1 .= " AND virus_family = '$virus_family' ";
			  
			  unless ($gene eq 'all') {
			   	  $query1 .= " AND digs_results.assigned_gene = '$gene' ";
			  }				
		   	  $query1 .= " AND bitscore >= $bitscore";
		
			  #print "\n\n\t $query1 \n\n"; die;
			my $sth = $dbh->prepare($query1);
			$sth->execute();	
			my @data;
			my $number = 0;
			while (my $row = $sth->fetchrow_arrayref) {
			
				#$devtools->print_array($row); die;
				my %data;
				$data{assigned_name}   = @$row[0];
				$data{assigned_gene}   = @$row[1];
				$data{virus_genus}     = @$row[2];
				$data{sequence_length} = @$row[3];
				$data{bitscore}        = @$row[4];
				push (@data, \%data);		
				$number++;
			}

			#$bitscore_gene_data{results} = \@data;
			$bitscore_gene_data{count}   = $number;
			$gene_data{$bitscore} = \%bitscore_gene_data;
		
		}
		$report_ref->{$gene} = \%gene_data;
	}
}

#***************************************************************************
# Subroutine:  get_genes_by_virus_group
# Description:  
#***************************************************************************
sub get_genes_by_virus_group {

	my ($db_ref, $taxonomic_level, $virus_family, $data_ref) = @_;

	my $dbh = $db_ref->{dbh};

	# Get the genes for this family
	my $query  = "SELECT DISTINCT assigned_gene";
	   $query .= " FROM $eve_data_table_name, digs_results";
       $query .= " WHERE assigned_name = $eve_data_table_name.virus_species ";
       $query .= " AND   $eve_data_table_name.$taxonomic_level = '$virus_family' ";
       $query .= " ORDER BY assigned_gene;";
	#print "\n\n\t $query \n\n";
	
	# Get data
	my $sth = $dbh->prepare($query);
	$sth->execute();	
	my @virus_family_genes;
	while (my $row = $sth->fetchrow_arrayref) {
		
		#$devtools->print_array($row); die;
		my $gene = @$row[0];
		if ($gene and $verbose) {
			print "\n\t\t Family: $virus_family: gene name '$gene'";
		}

		push (@$data_ref, $gene);		
	}
}

############################################################################
# Writing reports (ClkNva)
############################################################################

#***************************************************************************
# Subroutine:  write_virus_family_report
# Description:  
#***************************************************************************
sub write_virus_family_report {

	my ($family_report, $family) = @_;

	print   "\n\n\t  Creating report for family '$family'";	
	my $all_gene_data_ref = $family_report->{'all'};	
	#$devtools->print_hash($all_gene_data_ref); die;

	my @text;
	push(@text, "\t var data = {\n");

	# Bitscore labels
	push(@text, "\t   labels: [");
	my @f_bitscores;
	foreach my $bitscore (@bitscore) {
		push (@f_bitscores, "'$bitscore'");
	}
	my $f_bitscore_line = join(',', @f_bitscores);
	push(@text, $f_bitscore_line);
	push(@text, " ],\n");
	
	# Counts
	push(@text, "\t   series: [ ");
	my @counts;
	foreach my $bitscore (@bitscore) {
		my $count_ref = $all_gene_data_ref->{$bitscore};
		my $count = $count_ref->{count};
		push (@counts, $count);
	}
	my $counts_line = join(',', @counts);
	push(@text, $counts_line);
	push(@text, " ],\n");

	push(@text, "\t };\n");

	my $text = join("", @text);

	print "\n\n$text";

}


############################################################################
# Matching
############################################################################

#***************************************************************************
# Subroutine:  check_eve_data_table_matching
# Description:  
#***************************************************************************
sub check_eve_data_table_matching {

	my ($db_ref) = @_;

	# Execute a select distinct on 'assigned_name' and 'assigned_gene'
	my $results_table = $db_ref->{digs_results_table};
	my @fields = qw [ assigned_name ];
	my @results;
	my $order = " ORDER BY assigned_name, scaffold, extract_start ";
	$results_table->select_distinct(\@fields, \@results, $order);
	
	# Find any sequences assigned to references absent from the EVE taxonomy table
	my @missing;
	my %distinct_names;
	foreach my $row_ref (@results) {

		my $assigned_name = $row_ref->{assigned_name};
		my $data_rows = check_exists_in_data_table($db_ref, $assigned_name);
		unless ($data_rows) {
			print  "\n\t virus called '$assigned_name' not found";
			push (@missing, "$assigned_name\n");
		}
		elsif ($data_rows > 1) {
			print "\n\t multiple viruses called '$assigned_name'";
		}
	}
	
	# Output the list
	my $num_missing = scalar @missing;
	if ($num_missing) {
		my $missing_file = "missing_from_eve_data.txt";
		$fileio->write_file($missing_file, \@missing);
		print "\n\t Identified missing taxa. Exiting...\n\n\n";
		exit;	
	}
	else {
		print "\n\n\t  The '$eve_data_table_name' is matched to data in 'digs_results'";
		sleep 1;
	}

}

#***************************************************************************
# Subroutine:  check_exists_in_data_table
# Description:  
#***************************************************************************
sub check_exists_in_data_table {

	my ($db_ref, $assigned_name) = @_;

	
	my $dbh = $db_ref->{dbh};
	my $query = "SELECT virus_species, virus_class, virus_family, virus_genus";
	$query .= " FROM $eve_data_table_name";
    $query .= " WHERE virus_species = '$assigned_name' ";
	#print $query; die;
	my $sth = $dbh->prepare($query);
	$sth->execute();	
	my @data;
	while (my $row = $sth->fetchrow_arrayref) {
		
		#$devtools->print_array($row); exit;
		push (@data, $row);	
		my $species = @$row[0];
		my $family  = @$row[2];
		if ($species and $verbose) {
			print "\n\t\t Assigned name: '$assigned_name' = family '$family'";
		}
	
	}
	my $number = scalar @data;
	return $number;

}

############################################################################
# Taxonomy
############################################################################

#***************************************************************************
# Subroutine:  get_taxonomy_from_name_list
# Description:  
#***************************************************************************
sub get_taxonomy_from_name_list {

	my ($infile) = @_;

    my %db;
	load_ncbi_taxonomy(\%db);

	my @infile;
	$fileio->read_file($infile, \@infile);
	my @tax_strings;
	my @not_found;
	my $header = "Species\tClass\tSuperorder\tOrder\tFamily\tGenus";
	push (@tax_strings, "$header\n");
	foreach my $species (@infile) {
	
		chomp $species;
		my %taxonomy;
		my $tax_string = get_species_taxonomy_data($species, \%taxonomy, \%db);
		unless ($tax_string) {
			push (@not_found, "$species\n");
		}
		else {

			#$devtools->print_hash(\%taxonomy); die;
			my @output_string;

			my $class      = $taxonomy{class};
			my $superorder = $taxonomy{superorder};
			my $order      = $taxonomy{order};
			my $family     = $taxonomy{family};
			my $genus     = $taxonomy{genus};
		
			if ($species) {
				$species =~ s/_/ /g;
				push(@output_string, $species);
			} else {
				push(@output_string, "Unclassified");
			}
			
			if ($class) {
				push(@output_string, $class);
			} else {
				push(@output_string, "Unclassified");
			
			}
			if ($superorder) {
				push(@output_string, $superorder);
			} else {
				push(@output_string, "Unclassified");
			}
			
			if ($order) {
				push(@output_string, $order);
			} else {
				push(@output_string, "Unclassified");
			}
			
			if ($family) {
				push(@output_string, $family);
				} else {
				push(@output_string, "Unclassified");
			}
			
			if ($genus) {
				push(@output_string, $genus);
			} else {
				push(@output_string, "Unclassified");
			}

			my $output_string = join("\t", @output_string);
			print "\n\t ## Species:\t$output_string";	

			push (@tax_strings, "$output_string\n");
		}
			
	}
	
	my $outfile = $infile . '.names.txt';
	$fileio->write_file($outfile, \@tax_strings);

	my $notfound = $infile . '.notfound.txt';
	$fileio->write_file($notfound, \@not_found);


}

#***************************************************************************
# Subroutine:  format_ncbi_fasta_for_digs
# Description: create a screening datbase 
#***************************************************************************
sub format_ncbi_fasta_for_digs {

	my $path = $ENV{NCBI_VIRUS_REFS};	
	my @files;
	$fileio->read_directory_tree_leaves_simple($path, \@files);

    my %db;
	load_ncbi_taxonomy(\%db);

	# Show what the exclusion settings are for phages
	if ($exclude_phage) {
		print "\n\t The setting to exclude phages is ON";
		sleep 2;
	}
	else {
		print "\n\t The setting to exclude phages is OFF";
		sleep 2;
	}

	# Iterate through the files
	my %sorted;
	my $i = 0;
	my @digs_fasta;
	my @tax_table;
	my @unclassified;
	foreach my $file_ref (@files) {
		
		my $file       = $file_ref->{file};
		my $file_path = $file_ref->{path};
		my @fasta;
		$fileio->read_fasta($file_path, \@fasta);
		foreach my $seq_hash (@fasta) {

			# Extract information about this sequence
			my $header   = $seq_hash->{header};
			my $sequence = $seq_hash->{sequence};
			$header=~ s/\//-/g;
			my @split1 = split(/\[/, $header);
			my $name = pop @split1;
			
			# Exclude sequences based on flags
			if ($exclude_phage and $name =~ /phage/) {
				next;
			}
			#if ($exclude_phage and $name =~ /mimivirus/) {
			#	next;
			#}
			
			$i++;
			#print "\n\t # $i:  $header";
			
			
			# Get taxonomy data from NCBI taxonomy
			$name =~ s/]//g;
			my %taxonomy;
			my $tax_string = get_species_taxonomy_data($name, \%taxonomy, \%db);
			$name =~ s/\(/-/g;
			$name =~ s/\)/-/g;
			$name =~ s{\A\s*}{};
			$name =~ s{\s*\z}{};
			$name =~ s/ /-/g;

			my $header2 = pop @split1;
			my @split2 = split(/\|/, $header2);
			my $gene   = pop @split2;
			$gene =~ s{\A\s*}{};
			$gene =~ s{\s*\z}{};
			$gene=~ s/ /-/g;

			$name =~ s/---/-/g; # Clean up
			$gene =~ s/--/-/g;  # Clean up
			
			
			# Normalise gene names 
			my $normalised_gene = normalise_gene_names($name, $gene);

			# Create header
	
			# Create fasta
			print "\n\t # $i DIGS $name ($normalised_gene)";
			#print "\n\t # $i DIGS $name ($normalised_gene) '$tax_string'";
			my $digs_header = $name . '_' . $normalised_gene;
			my $digs_fasta = ">$digs_header\n$sequence\n";	
			push (@digs_fasta, $digs_fasta);

			# Sort by taxonomic level - e.g. family, genus etc)
			add_sequence_to_taxonomy_group_fasta(\%sorted, $name, $gene, $tax_string, $digs_fasta);
		
			# Add the details of this sequence to taxonomy table
			if ($tax_string) {
				create_taxonomy_table_line($name, $gene, $tax_string, \@tax_table);
			}
			else {
				push(@unclassified, "$name\t$gene\n");
			}
		}
	}

	# Write taxonomy table
	$fileio->write_file('ncbi_virus_taxonomy.txt', \@tax_table);
	# Write out the table of viruses where we got not taxonomic info
	$fileio->write_file('unclassified.txt', \@unclassified);
	# Write the FASTA
	$fileio->write_file('NCBI_viruses.faa', \@digs_fasta);

	# Write subgroup files
	#my @subdivisions = keys %sorted;
	#foreach my $subdivision (@subdivisions) {
	#	my $filename = $subdivision . '.DIGS.faa';
	#	my $array_ref = $sorted{$subdivision};
	#	#$fileio->write_file($filename, $array_ref);
	#}
}

#***************************************************************************
# Subroutine:  format_ncbi_fasta_for_digs2
# Description: simple header conversion routine
#***************************************************************************
sub format_ncbi_fasta_for_digs2 {

	my ($infile) = @_;

	my $species_name;
	my $species_name_question = "\n\t What is the name of the species?";
	$species_name = $console->ask_question($species_name_question);

	my @fasta;
	my @digs_fasta;
	$fileio->read_fasta($infile, \@fasta);
	my $i;
	foreach my $seq_hash (@fasta) {
		
		# Extract information about this sequence
		$i++;
		my $header   = $seq_hash->{header};
		my $sequence = $seq_hash->{sequence};
		$header=~ s/\//-/g;
		my @split1 = split(/\[/, $header);
		#$devtools->print_array(\@split1); #die;
		shift @split1;
		my $gene_element = shift @split1;
		$gene_element =~ s/]//g;
		$gene_element =~ s/\s+$//;
		my @gene_element = split('=', $gene_element);
		my $gene = pop @gene_element;
	   	      
		# Normalise gene name
		#my $normalised_gene = normalise_gene_names($name, $gene);

		# Create fasta
		#print "\n\t # $i DIGS $name ($normalised_gene)";
		print "\n\t # $i DIGS NAME: $species_name ($gene)";
		my $digs_header = $species_name . '_' . $gene;
		my $digs_fasta = ">$digs_header\n$sequence\n";	
		push (@digs_fasta, $digs_fasta);

	}

	my $outfile = $infile . '.out.txt';
	$fileio->write_file($outfile, \@digs_fasta);
	
}

#***************************************************************************
# Subroutine:  create_taxonomy_table_line
# Description: create taxonomy table line
#***************************************************************************
sub create_taxonomy_table_line {

	my ($name, $gene, $tax_string, $tax_table) = @_;
	
	#print "\n\t$tax_string";
	my @tax_string = split(':', $tax_string);
	
	# Iterate through the taxonomy string elements
	my $baltimore = 'Unclassified';
	my $family    = 'Unclassified';
	foreach my $element (@tax_string) {
		if ($element eq 'root') { # Skip root elements
			next; 
		}
		# Baltimore
		if ($element =~ /ssDNA viruses/
		or  $element =~ /dsDNA viruses/
		or  $element =~ /dsRNA viruses/
		or  $element =~ /ssRNA viruses/
		or  $element =~ /ssRNA positive-strand viruses/
		or  $element =~ /ssRNA negative-strand viruses/
		or  $element =~ /Retro-transcribing viruses/
		) {
			$baltimore = $element;
		}
		if ($element =~ /viridae/) {
			$family = $element;
		}
	}

	# Special cases
	if ($name eq 'Citrus-endogenous-pararetrovirus') {
		$family = 'Caulimoviridae'
	}
	if ($name eq 'Chickpea-chlorosis-virus-A') {
		$baltimore = 'ssRNA positive-strand viruses, no DNA stage';
		$family    = 'Luteoviridae';
	}
	if ($name eq 'Dragonfly-associated-mastrevirus') {
		$baltimore = 'ssDNA viruses';
		$family    = 'Geminiviridae';
	}
	if ($name eq 'Arhar-cryptic-virus-I') {
		$baltimore = 'dsDNA viruses, no RNA stage';
		$family    = 'Partitiviridae';
	}

	$baltimore =~ s/unclassified //;

	my $seq_details = "$name\t$baltimore\t$family\t$gene\n";
	push (@$tax_table, $seq_details);

}

#***************************************************************************
# Subroutine:  normalise_gene_names
# Description: normalise names
#***************************************************************************
sub normalise_gene_names {

	my ($name, $gene) = @_;
	
	# Circoviridae
	if ($name =~ /[Cc]ircovirus/ or 
		$name =~ /[Cc]yclovirus/ or 
		$name =~ /Beak-and-feather-disease-virus/
		) {
		if ($gene =~ /[Rr]ep/) {
			$gene = 'Rep';
		}
	}
	# Parvoviridae
	if ($name =~ /[Aa]deno-associated/ or
		$name =~ /[Pp]arvovirus/  or
		$name =~ /[Bb]ufavirus/  or
		$name =~ /Minute-virus-of-mice/) {
		if ($gene =~ /[Cc]ap/ or
			$gene =~ /ORF2/ or 
			$gene =~ /VP/) { 
			$gene = 'VP';
		}	
		if ($gene =~ /[Rr]ep/ or
			$gene =~ /REP/ or 
			$gene =~ /ORF1/ or 
			$gene =~ /nonstructural/ or 
			$gene =~ /non-structural/ or 
			$gene =~ /NS/ ) { 
			$gene = 'NS';
		}	
	}
	# Filoviridae
	if ($name =~ /Reston-ebolavirus/ or
		$name =~ /Bundibugyo-virus/  or
		$name =~ /Sudan-ebolavirus/) {
		if ($gene =~ /RNA-dependent-RNA-polymerase/) {
			$gene = 'L-polymerase';
		}
	}
	if ($name =~ /Reston-ebolavirus/ and $gene =~ /polymerase-complex-protein/) {
		$gene = 'L-polymerase';
	}

	return $gene;
}

#***************************************************************************
# Subroutine:  add_sequence_to_taxonomy_group_fasta
# Description: add a sequence to each taxonomic set (e.g. genus, family)
#              it belongs in
#***************************************************************************
sub add_sequence_to_taxonomy_group_fasta {

	my ($sorted_ref, $name, $gene, $tax_string, $digs_fasta) = @_;

	#print "\n\t$tax_string";
	my @tax_string = split(':', $tax_string);
	
	# Remove species (dont want to write out files at species level 
	pop @tax_string;
	pop @tax_string;
	
	# Iterate through the taxonomy string elements
	foreach my $element (@tax_string) {
		if ($element eq 'root') { # Skip root elements
			next; 
		}
		else {
			$element =~ s/\//-/g;
			$element =~ s/ /_/g;
			if ($sorted_ref->{$element}) {
				my $fasta_array_ref = $sorted_ref->{$element};	
				push (@$fasta_array_ref, $digs_fasta);
			}
			else {
				my @fasta;
				push (@fasta, $digs_fasta);
				$sorted_ref->{$element} = \@fasta;
			}
		}
	}
}

#***************************************************************************
# Subroutine:  get_species_taxonomy_data
# Description: 
#***************************************************************************
sub get_species_taxonomy_data {

    my ($species, $taxonomy_ref, $db_ref) = @_;

	# Get tables
	my $ncbi_names = $db_ref->{ncbi_names}; 
	my $ncbi_nodes = $db_ref->{ncbi_nodes}; 

	# Remove subspecies element if present
	my @elements  = split("_", $species);
	my $element_count = scalar @elements;
	if ($element_count eq 3) {
		pop @elements; # Remove the sub-species element
		$species = join ("_", @elements);
	}
	$species =~ s/_/ /g;
	print "\n\t # Getting taxonomy for '$species'";
	
	# Iterate through the species
	$species =~ s/-/ /g;
	my @fields = qw [ tax_id name_txt unique_name name_class ];
	my $where = "WHERE name_txt = '$species'";
	my %data;
	$ncbi_names->select_row(\@fields, \%data, $where);
	my $tax_id = $data{tax_id};
	unless ($tax_id) { 
		print "\n\t No taxon_id found for '$species'" ; 
		return 0;
	}
	my $limit = 'class';
	my $i = 0;
	retrieve_taxonomy($db_ref, $tax_id, $taxonomy_ref, $limit, $i);
	
	my $taxa_string = $taxonomy_ref->{taxonomy_string};
	#print "\n\t TAXA STRING $taxa_string";
	return $taxa_string;
	
}

############################################################################
# SECTION: INTERACTING WITH THE NCBI TAXONOMY DB
############################################################################

#***************************************************************************
# Subroutine:  load_ncbi_taxonomy 
# Description: 
#***************************************************************************
sub load_ncbi_taxonomy {

    my ($db_obj) = @_;

	my $dbh = DBI->connect("dbi:mysql:$db_name:$server", $mysql_username, $mysql_password);
    
    # Main Screening DB tables
    load_ncbi_names($db_obj, $dbh);  
    load_ncbi_nodes($db_obj, $dbh);  

}

#***************************************************************************
# Subroutine:  load_ncbi_names
# Description: load database  
#***************************************************************************
sub load_ncbi_names {

    my ($db_obj, $dbh) = @_;

    # Definition of the table
    my %names_fields = (
        tax_id                  => 'int',
        name_txt                => 'varchar',
        unique_name             => 'varchar',
        name_class              => 'varchar',
    );   
    my $ncbi_names = MySQLtable->new('ncbi_names', $dbh, \%names_fields);
    $db_obj->{ncbi_names} = $ncbi_names;
}

#***************************************************************************
# Subroutine:  load_ncbi_nodes
# Description: load database  
#***************************************************************************
sub load_ncbi_nodes {

    my ($db_obj, $dbh) = @_;

    # Definition of the table
    my %nodes_fields = (
        tax_id                        => 'int',
        parent_tax_id                 => 'int',
        rank                          => 'varchar',
        embl_code                     => 'varchar',
        division_id                   => 'int',
        inherited_div_flag            => 'int',
        inherited_gc_flag             => 'int',
        mitochondrial_genetic_code_id => 'int',
        inherited_mgc_flag            => 'int',
        genbank_hidden_flag           => 'int',
        hidden_subtree_root_flag      => 'int',
        comments                      => 'varchar',
    );
    my $ncbi_nodes = MySQLtable->new('ncbi_nodes', $dbh, \%nodes_fields);
    $db_obj->{ncbi_nodes} = $ncbi_nodes;
}

#***************************************************************************
# Subroutine:  retrieve_taxonomy 
# Description:  
#***************************************************************************
sub retrieve_taxonomy {

	my ($db_obj, $tax_id, $taxonomy, $limit_rank, $i) = @_;

	$i++;
	if ($i > 30) { return; } # Prevent deep recursion

	#print "\n\t Getting data for $tax_id";	
	my $ncbi_names = $db_obj->{ncbi_names}; 
	my $ncbi_nodes = $db_obj->{ncbi_nodes}; 
	#sleep 1;

	my %node_data;
	my @node_fields = qw [ parent_tax_id rank ]; 
	my $where  = "WHERE tax_id = '$tax_id' ";
	$ncbi_nodes->select_row(\@node_fields, \%node_data, $where);
	my $rank           = $node_data{rank};
	my $parent_tax_id   = $node_data{parent_tax_id};
	unless ($parent_tax_id and $rank) {
		print "\n\t got nothing for id $tax_id";
	}

	my %data;
	my @fields = qw [ tax_id name_txt unique_name name_class ];
	$where = "WHERE tax_id = '$tax_id' AND name_class = 'scientific name' ";
	$ncbi_names->select_row(\@fields, \%data, $where);
	my $name = $data{name_txt};

	unless ($tax_id and $name) {
		die;
		return;
	}
	else {

		unless ($rank eq 'no rank') {
			$taxonomy->{$rank} = $name;
		}

		my $taxa_string = $taxonomy->{taxonomy_string};
		unless ($taxa_string) {
			$taxa_string = $name;
		}
		else {
			$name .= ":$taxa_string";
		}
		$taxonomy->{taxonomy_string} = $name;

		#if ($rank eq $limit_rank) {
		if ($rank eq 'domain') {
			return;
		}
		elsif ($parent_tax_id) {
			#print "\n id $tax_id = $name ($rank), parent = $parent_tax_id)";
			retrieve_taxonomy($db_obj, $parent_tax_id, $taxonomy, $limit_rank, $i);
		}
	}
}

#***************************************************************************
# Subroutine:  crunch_loci
# Description:  
#***************************************************************************
sub crunch_loci {

	my ($file) = @_;

	# Read the file
	print "\n\t # Reading file '$file' for crunch_loci";
	my @file;
	$fileio->read_file($file, \@file);
	
	my $header = shift @file;
	my @header = split("\t", $header);
	my $i = 0;	
	foreach my $column_heading (@header) {
		$i++;
		print "\n\t $i: $column_heading";		
	}
	
	# Hard coded column references
	my $scaffold_i = 3;
	my $start_i = 4;
	my $end_i = 5;

	# Iterate through the file
	my $previous_scaffold;
	my $previous_start;
	my $previous_end;
	my $scaffold_high_end = undef;

	my @output;
	
	unshift (@header, 'merge?');
	my $new_header = join ("\t", @header);
	push (@output, $new_header);

	# Process file	
	foreach my $line (@file) {

		my @line = split("\t", $line);
	
		my $merge = undef;
		my $host = $line[0];
		my $current_scaffold = $line[$scaffold_i];
		my $current_start    = $line[$start_i];
		my $current_end      = $line[$end_i];
		
		unless ($previous_scaffold) {
			$previous_scaffold = $current_scaffold;
			$previous_start    = $current_start;
			$previous_end      = $current_end;
			next;
		}
				
		if ($current_scaffold eq $previous_scaffold) { # Do comparison

			if ($current_start < $scaffold_high_end) {
				$merge = "yes";
				
			}
			else {
				$merge = "no";
			}

			print "\n\n\t # HOST: $host ";
			print "\n\t # Previous scaffold: $previous_scaffold ";
			print "\n\t # Current scaffold: $current_scaffold ";
			print "\n\t # PREVIOUS: start $previous_start  : end $previous_end ";
			print "\n\t # CURRENT:  start $current_start : end $current_end  ";
			print "\n\t # Highest end: $scaffold_high_end ";
			print "\n\t # MERGE?: $merge ";

			if ($current_end > $scaffold_high_end) {
				$scaffold_high_end = $current_end;
			}

		}
		else {	
			$scaffold_high_end = $current_end;
			$merge = "NEW-SCAFFOLD";

		}

		$previous_start    = $current_start;
		$previous_end      = $current_end;
		$previous_scaffold = $current_scaffold;
		unshift (@line, $merge);
		my $new_line = join ("\t", @line);		
		push (@output, $new_line);
		$merge = undef;
		
	}

	my $original_number = scalar @file;
	my $new_number = scalar @output;
	print "\n\t Origninal lines $original_number \n\n";
	print "\n\t new_number lines $new_number \n\n";
	$fileio->write_file("crunched.txt", \@output);

	# Do merge
	do_merge(\@output);

}

#***************************************************************************
# Subroutine:  do_merge
# Description:  
#***************************************************************************
sub do_merge {

	my ($output_ref) = @_;
	
	# Hard coded column references
	my $scaffold_i = 4;
	my $start_i = 5;
	my $end_i = 6;

	# Iterate through the file
	my $previous_merge = undef;
	my $scaffold_high_end = undef;
	my $scaffold_low_start = undef;
	my $previous_line;
	my $header = shift @$output_ref;

	# Do merge
	my @merged;
	my $combine = undef;
	my $num_lines;
	foreach my $line (@$output_ref) {

		my @line = split("\t", $line);	
		
		
		my $merge = $line[0];
		my $host  = $line[1];
		my $current_scaffold = $line[$scaffold_i];
		my $current_start    = $line[$start_i];
		my $current_end      = $line[$end_i];

	
		# If its a new scaffold or new chain, record last and reset
		if ($merge eq 'NEW-SCAFFOLD' or $merge eq 'no') {		

			my $merge_line = create_merged_line($previous_line, $scaffold_low_start, $scaffold_high_end, $combine);
			$num_lines++;
			#print "\n\t Storing merged line number $num_lines '$merge'";

 			push (@merged, $merge_line);

			$combine = undef;
			$scaffold_low_start  = $current_start;
			$scaffold_high_end   = $current_end;
								
		}
		elsif ($merge eq 'yes') {
		
			# Set flag to merge
			$combine = 'yes'; 
			print "\n\t Merging line";
		
			if ($current_end > $scaffold_high_end) {
			
				print "\n\t setting end to '$current_end";
				$scaffold_high_end = $current_end;
			}
			
		}
		else { die; }
			
		# Record current line as the last
		$previous_line = $line;
	}


	my $merge_line = create_merged_line($previous_line, $scaffold_low_start, $scaffold_high_end, $combine);
	$num_lines++;
 	push (@merged, $merge_line);	
	

	my $number_after_merge = scalar @merged;
	print "\n\t Crunched! there are $number_after_merge lines now! \n\n";

	$fileio->write_file("crunched2.txt", \@merged);

}

#***************************************************************************
# Subroutine:  create_merged_line
# Description:  
#***************************************************************************
sub create_merged_line {

	my ($line, $scaffold_low_start, $scaffold_high_end, $combine) = @_;


	print "\n\t $scaffold_low_start";
	print "\n\t $scaffold_high_end";
	if ($combine) {
		print "\n\t combine";
	}
	else {
		print "\n\t don't combine";
	
	}

	my $merge_line;
	unless ($combine) {	

		my @line = split("\t", $line);			
		my $seq_len = $scaffold_high_end - $scaffold_low_start;
		$line[0] = 'No MERGED';
		$merge_line = join("\t", @line);
	}
	else {

		my @line = split("\t", $line);			
		my $seq_len = $scaffold_high_end - $scaffold_low_start;
		$line[0] = 'MERGED';
		$line[5] = $scaffold_low_start;
		$line[6] = $scaffold_high_end;
		$line[7] = $seq_len;		
		$merge_line = join("\t", @line);
	}

	return $merge_line;
}

############################################################################
# Title and help display
############################################################################

#***************************************************************************
# Subroutine:  show_title
# Description: show command line title blurb 
#***************************************************************************
sub show_title {

	$console->refresh();
	my $title       = 'digs_for_eves_tools.pl';
	my $description = 'DIGS for EVEs database PERL utility';
	my $author      = 'Robert J. Gifford';
	my $contact	    = '<robert.gifford@glasgow.ac.uk>';
	$console->show_about_box($title, $version, $description, $author, $contact);
}

#***************************************************************************
# Subroutine:  show_help_page
# Description: show help page information
#***************************************************************************
sub show_help_page {

	# Initialise usage statement to print if usage is incorrect
	my ($HELP)  = "\n\t Usage: $0 -m=[option] -i=[control file]\n";
        $HELP  .= "\n\t ### Main functions\n"; 
        $HELP  .= "\n\t -m=1  parse ncbi reference virus folder to DIGS input"; 
        $HELP  .= "\n\t -m=2  parse an ncbi 'coding sequence' download to DIGS input"; 
        $HELP  .= "\n\t -m=3  get taxonomy table from species name list"; 
        $HELP  .= "\n\t -m=4  create DIGS-for-EVEs report\n\n"; 

	print $HELP;
}

############################################################################
# End of file 
############################################################################
