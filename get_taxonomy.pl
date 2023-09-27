#!/usr/bin/perl -w
############################################################################
# Script:      get_taxonomy.pl
# Creator:     R.J. Gifford
# Description: DIGS for EVEs perl tools
# History:     Version 1.0
############################################################################

# Include the PERL module library for DIGS 
use lib ($ENV{DIGS_HOME2} . '/modules/'); 

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
use Interface::MySQLtable;   

############################################################################
# Paths & Globals
############################################################################

# Create a unique process ID for this DIGS screening process
my $pid  = $$;
my $time = time;
my $process_id  = $pid . '_' . $time;
my $version = '1.0';

# Taxonomy database connection globals
my $db_name  = 'rg_ncbi_taxonomy';
my $server   = 'localhost';

############################################################################
# Instantiations
############################################################################

# Base utilites
my $fileio     = FileIO->new();
my $console    = Console->new();
my $devtools   = DevTools->new();

my $mysql_username = ($ENV{DIGS_MYSQL_USER}); 
my $mysql_password = ($ENV{DIGS_MYSQL_PASSWORD}); 

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
			    'infile|i=s'     => \$infile,
	) or die $USAGE;

	show_title();
	
	# Set flags based on options received
	if ($mode) { 
		if ($mode eq 1) {
			get_taxonomy_from_name_list($infile);
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

############################################################################
# Basic functions (IO, console interaction etc)
############################################################################

#***************************************************************************
# Subroutine:  read_file
# Description: read an input file to an array
# Arguments:   $file: the name of the file to read
#              $array_ref: array to copy to
#***************************************************************************
sub read_file {

	my ($file, $array_ref) = @_;

	unless (-f $file) {
		if (-d $file) {
			print "\n\t Cannot open file \"$file\" - it is a directory\n\n";
			die;
		}
		else {
			print "\n\t Cannot open file \"$file\"\n\n";
			die;
		}

 	}
	unless (open(INFILE, "$file")) {
		die "\n\t Cannot open file \"$file\"\n\n";
	}
	@$array_ref = <INFILE>;
	close INFILE;

	my $lines = scalar @$array_ref;
	unless ($lines) {
		die "\n\t Unable to read infile '$file'\n\n";
	}
	elsif ($lines eq 1) {
		print "\n\t Single line of data was read from infile '$file'";
		die "\n\t Check line break formatting\n\n";
	}
	
	return 1;
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
	my $title       = 'get_taxonomy.pl';
	my $description = 'PERL utility for use with NCBI taxonomy';
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
	my ($HELP)  = "\n\t Usage: $0 -m=[option] -i=[input file]\n";
        $HELP  .= "\n\t ### Main functions\n"; 
        $HELP  .= "\n\t -m=1  get taxonomy table from species name list\n\n"; 

	print "\n\t ### The input file for this script should be a list of species names (Latin binomial)\n";

	print $HELP;
}

############################################################################
# End of file 
############################################################################
