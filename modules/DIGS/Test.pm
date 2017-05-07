#!usr/bin/perl -w
############################################################################
# Module:      Test.pm   
# Description: Tests for the DIGS tool
# History:     December  2017: Created by Robert Gifford 
############################################################################
package Test;

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

# DIGS test database connection globals
my $server   = 'localhost';
my $user     = ($ENV{DIGS_MYSQL_USER});
my $password = ($ENV{DIGS_MYSQL_PASSWORD}); 

# Test globals
my $target_path = './test/targets/';

1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create new Test 'object'
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
# TESTS
############################################################################

#***************************************************************************
# Subroutine:  run_tests
# Description:  
#***************************************************************************
sub run_tests {

	my ($self) = @_;

	my $digs_obj = $self->{digs_obj};


 	# Show title
	$digs_obj->show_title();  

	# Capture the path for output
	unless ($ENV{'DIGS_OUTPUT'}) {
		print  "\n\t Required environment variable '\$DIGS_OUTPUT' is undefined\n\n";
		exit;
	}
	
	# Load the 'digs_test' database
	$digs_obj->{mysql_server}   = $server;
	$digs_obj->{mysql_username} = $user;
	$digs_obj->{mysql_password} = $password;
	$digs_obj->{genome_use_path} = $target_path;
	my $initialise_obj = Initialise->new($digs_obj);
	$initialise_obj->initialise_screening_db($digs_obj, 'digs_test_screen');
	
	# Flush the 'digs_test' database
	my $db = $digs_obj->{db}; # Get the database reference
	$db->flush_screening_db();

	# Initialise a DIGS object for these tests
	$digs_obj->{output_path}     = $ENV{'DIGS_OUTPUT'}; 
	$initialise_obj->create_output_directories($digs_obj);

	# TODO - do we need this?
	my $loader_obj = ScreenBuilder->new($digs_obj); # Create the ScreenBuilder object
	$digs_obj->{loader_obj} = $loader_obj;

	# Do a screen using test control file and synthetic target data
	print "\n\n\t ### Running DIGS tests ~ + ~ + ~ \n";
	$self->run_test_1();
	$self->run_test_2();
	$self->run_test_3();
	$self->run_test_4();
	$self->run_test_5();
	$self->run_test_6();
	$self->run_test_7();
	die;
	#$self->run_test_10();

	# Print finished message
	print "\n\n\t ### Tests completed ~ + ~ + ~\n\n\n";

	# Remove the output directory
	my $output_dir = $digs_obj->{report_dir};
	if ($output_dir) {
		my $command1 = "rm -rf $output_dir";
		system $command1;
	}
	else { die; }
}

#***************************************************************************
# Subroutine:  run_test_1
# Description:  
#***************************************************************************
sub run_test_1 {

	my ($self) = @_;

	print "\n\t ### TEST 1: Running nucleotide screen against synthetic data ~ + ~ + ~ \n\n";
	my $digs_obj   = $self->{digs_obj};
	my $loader_obj = $digs_obj->{loader_obj};

	# Set the parameters for this test (directly rather than using control file)
	$loader_obj->{seq_length_minimum}   = 100;
	$loader_obj->{defragment_range}     = 100;
	$loader_obj->{query_na_fasta}       = "./test/probes/korv_test1.fna"; 
	$loader_obj->{reference_na_fasta}   = "./test/references/korv_test1.fna"; 
	$loader_obj->{bitscore_min_blastn}  = 100;
	$loader_obj->{seq_length_minimum}   = 50;
	my @test_targets;
	my $test1_target_path = "species/datatype/version/artificial_test1_korv.fa";
	push (@test_targets, $test1_target_path);
	$loader_obj->{target_paths} = \@test_targets;
	
	# Run this test
	$digs_obj->{defragment_mode}  = 'defragment';
	$digs_obj->{defragment_range} = 100;
 	$digs_obj->{seq_length_minimum} = $loader_obj->{seq_length_minimum};
	$digs_obj->{bitscore_minimum}   = $loader_obj->{bitscore_min_blastn};
	unless ($self->{bitscore_minimum}) {
		$digs_obj->{bitscore_minimum} = $loader_obj->{bitscore_min_tblastn};
	}
	$digs_obj->{bitscore_minimum}   = 100;

	# Initialise a DIGS object
	my $initialise_obj = Initialise->new($digs_obj);
	$initialise_obj->setup_for_a_digs_run($digs_obj);
	
	# Perform DIGS
	$digs_obj->perform_digs();
	#$devtools->print_hash($self); die;

	# Check that we got expected result
	# For this test it is two hits
	# Hit 1: start KoRV, LTR: 200   end 703   in -ve orientation 
	# Hit 2: start KoRV, LTR: 10967 end 11470 in -ve orientation 
	my $db = $digs_obj->{db}; # Get the database reference
	my $results_table = $db->{digs_results_table};
	my @data;
	my @fields = qw [ assigned_gene assigned_name extract_start extract_end ];
	my $sort = " ORDER BY scaffold, extract_start ";
	$results_table->select_rows(\@fields, \@data, $sort);
	my $result1_ref = shift @data;		
	my $result2_ref = shift @data;
	my $correct_result = 1;
	unless ($result1_ref->{assigned_gene} eq 'LTR'
	   and  $result2_ref->{assigned_gene} eq 'LTR')  { $correct_result = undef; }
	unless ($result1_ref->{assigned_name} eq 'KoRV'
	   and  $result2_ref->{assigned_name} eq 'KoRV') { $correct_result = undef; }
	unless ($result1_ref->{extract_start} eq 200
	   and  $result2_ref->{extract_start} eq 10967)  { $correct_result = undef; }
	unless ($result1_ref->{extract_end}   eq 703
	   and  $result2_ref->{extract_end}   eq 11470)  { $correct_result = undef; }
	if ($correct_result)  { print "\n\n\t  blastn screen test: ** PASSED **\n" }
	else                  { die   "\n\n\t  blastn screen test: ** FAILED **\n" }
	#$devtools->print_hash($result1_ref); $devtools->print_hash($result2_ref); die;
	sleep 1;
}

#***************************************************************************
# Subroutine:  run_test_2
# Description: Defragment results test negative (i.e. do not merge loci)
#***************************************************************************
sub run_test_2 {

	my ($self) = @_;

	my $digs_obj   = $self->{digs_obj};
	$digs_obj->{defragment_mode}  = 'defragment';
	$digs_obj->{defragment_range} = 100;
 
	## Check that defragment gives expected result (negative)
	print "\n\t ### TEST 2: Defragment results test negative  ~ + ~ + ~ \n";	

	# Construct WHERE statement
	my $where  = " WHERE organism      = 'species' ";
	$where    .= " AND target_datatype = 'datatype' ";
	$where    .= " AND target_version  = 'version' ";
	$where    .= " AND target_name     = 'artificial_test1_korv.fa' "; 
	my $path   = "species/datatype/version/artificial_test1_korv.fa";
	my $target_path  = './test/targets/'  . $path;
	my %settings;
	$settings{range}     = 100;
	$settings{start}     = 'extract_start';
	$settings{end}       = 'extract_end';
	$settings{where_sql} = $where;

	# Initialise the DIGS object
	my $initialise_obj = Initialise->new($digs_obj);
	$initialise_obj->setup_for_a_digs_run($digs_obj);

	# Get digs results ready for defragment process
	my @sorted_digs_results;
	$digs_obj->get_sorted_digs_results(\@sorted_digs_results, $where);
	my $num_hits = scalar @sorted_digs_results;
	print "\n\t\t # $num_hits digs results to defragment ";
	#$devtools->print_array(\@sorted_digs_results); die;
	$settings{defragment_loci} = \@sorted_digs_results;
	$settings{digs_obj} = $digs_obj;
 
	# Defragment results for this target file
	my $defrag_obj = Defragment->new($digs_obj);	
	my $num_new = $defrag_obj->defragment_target(\%settings, $target_path, 'digs_results');
	if ($num_new eq '0' )  { print "\n\n\t  Defragment negative test: ** PASSED **\n" }
	else                   { die   "\n\n\t  Defragment negative test: ** FAILED **\n" }
	sleep 1;
}
	
#***************************************************************************
# Subroutine:  run_test_3
# Description: Partially deleted pol peptide screen against synthetic data
#***************************************************************************
sub run_test_3 {

	my ($self) = @_;

	my $digs_obj   = $self->{digs_obj};
	$digs_obj->{defragment_mode}  = 'defragment';
	$digs_obj->{defragment_range} = 100;

	# Run a peptide screen
	print "\n\t ### TEST 3: Running partially deleted pol peptide screen against synthetic data ~ + ~ + ~ \n";
	my $test_ctl_file2 = './test/test3_erv_aa.ctl';
	my $loader_obj     = $digs_obj->{loader_obj};
	#$loader_obj->parse_control_file($test_ctl_file2, $self, 2);

	# Set the parameters for this test (directly rather than using control file)
	$loader_obj->{seq_length_minimum}   = 100;
	$loader_obj->{defragment_range}     = 100;
	$loader_obj->{query_aa_fasta}       = "./test/probes/korv_test3.faa"; 
	$loader_obj->{reference_aa_fasta}   = "./test/probes/korv_test3.faa"; 
	$loader_obj->{bitscore_min_tblastn}  = 100;
	$loader_obj->{seq_length_minimum}   = 50;
	my @test_targets;
	my $test1_target_path = "species/datatype/version/artificial_test1_korv.fa";
	push (@test_targets, $test1_target_path);
	$loader_obj->{target_paths} = \@test_targets;

	my $initialise_obj = Initialise->new($digs_obj);
	$initialise_obj->setup_for_a_digs_run($digs_obj);
	$digs_obj->perform_digs();

	my $db = $digs_obj->{db}; # Get the database reference
	my $results_table = $db->{digs_results_table};
	my @data;
	my @fields = qw [ assigned_gene assigned_name extract_start extract_end ];
	my $test3_where = " WHERE probe_type = 'ORF' ORDER BY scaffold, extract_start ";
	$results_table->select_rows(\@fields, \@data, $test3_where);
	my $result3_ref = shift @data;
	my $result4_ref = shift @data;
	unless ($result3_ref and $result4_ref) { die; };
	my $correct_result = 1;
	unless ($result3_ref->{assigned_gene} eq 'pol'
	   and  $result4_ref->{assigned_gene} eq 'pol') { $correct_result = undef; }
	unless ($result3_ref->{assigned_name} eq 'KoRV'
	   and  $result4_ref->{assigned_name} eq 'KoRV') { $correct_result = undef; }
	unless ($result3_ref->{extract_start} eq 5681
	   and  $result4_ref->{extract_start} eq 7481)   { $correct_result = undef; }
	unless ($result3_ref->{extract_end}   eq 7144
	   and  $result4_ref->{extract_end}   eq 9064)   { $correct_result = undef; }
	if ($correct_result)  { print "\n\n\t  tblastn test: ** PASSED **\n" }
	else                  { die   "\n\n\t  tblastn test: ** FAILED **\n" }
	#$devtools->print_hash($result1_ref); $devtools->print_hash($result2_ref); die;
	sleep 1;
}

#***************************************************************************
# Subroutine:  run_test_4
# Description: Defragment results test positive (i.e. do merge loci)
#***************************************************************************
sub run_test_4 {

	my ($self) = @_;

	## Check that defragment gives expected result	(should join gag and pol with range of 200)	
	print "\n\t ### TEST 4: Defragment results test positive  ~ + ~ + ~ \n";
 
 	my $digs_obj   = $self->{digs_obj};
	$digs_obj->{defragment_mode}  = 'defragment';
	$digs_obj->{defragment_range} = 500;

	# Construct WHERE statement
	my $where  = " WHERE probe_type = 'ORF' ";
	my $path   = "species/datatype/version/artificial_test1_korv.fa";
	my $target_path  = './test/targets/'  . $path;
	my %settings;
	$settings{range}     = 500;
	$settings{start}     = 'extract_start';
	$settings{end}       = 'extract_end';
	$settings{where_sql} = $where;

	# Get digs results ready for defragment process
	my @sorted_digs_results;
	$digs_obj->get_sorted_digs_results(\@sorted_digs_results, $where);
	my $num_hits = scalar @sorted_digs_results;
	print "\n\t\t # $num_hits digs results to defragment ";
	#$devtools->print_array(\@sorted_digs_results); die;
	$settings{defragment_loci} = \@sorted_digs_results;
	$settings{digs_obj} = $digs_obj;
 
	# Defragment results for this target file
	my $defrag_obj = Defragment->new($digs_obj);	
	my $num_new = $defrag_obj->defragment_target(\%settings, $target_path, 'digs_results');
	my $db = $digs_obj->{db}; # Get the database reference
	my $results_table = $db->{digs_results_table};	
	my @data;
	my @fields = qw [ assigned_gene assigned_name extract_start extract_end ];
	my $sort = " ORDER BY scaffold, extract_start ";
	$results_table->select_rows(\@fields, \@data, $sort);
	my $num_rows = scalar @data;

	my $fail = undef;
	my $result_ref = shift @data;
	unless ($result_ref->{extract_start} eq 5681 and $result_ref->{extract_start} eq 9064) { 
	   $fail = 1;
	}
	if ($num_new  eq '0' ) { 
		die   "\n\t  Defragment positive test: ** FAILED ** No merge \n";
		$fail = 1;
	}
	elsif ($num_rows ne '3')  { 
		die   "\n\t  Defragment positive test: ** FAILED ** No cleanup in digs table \n";
		$fail = 1;
	}
	else {
		print "\n\t  Defragment positive test: ** PASSED **\n";
	}
	
	sleep 2;
}

#***************************************************************************
# Subroutine:  run_test_5
# Description: gag + env peptide screen
#***************************************************************************
sub run_test_5 {

	my ($self) = @_;

	# Run the second peptide screen
	print "\n\t ### TEST 5: Running gag + env peptide screen (entails merge of result rows) ~ + ~ + ~ ";

 	my $digs_obj   = $self->{digs_obj};
	my $loader_obj = $digs_obj->{loader_obj};
	my $db = $digs_obj->{db}; # Get the database reference
	my $results_table = $db->{digs_results_table}; # Get the database reference
	$results_table->flush();
	my $test5_path = './test/tabular/test5.txt';;
	print "\n\t ### Flushed digs_results table & now uploading data from file '$test5_path'";
	$db->upload_data_to_digs_results($test5_path);
	print "\n\t ### Data uploaded, starting from point of having successfully conducted tests 1,2,3, & 4\n";

	# Set the parameters for this test (directly rather than using control file)
	$loader_obj->{seq_length_minimum}   = 100;
	$loader_obj->{defragment_range}     = 100;
	$loader_obj->{query_aa_fasta}       = "./test/probes/korv_test5.faa"; 
	$loader_obj->{reference_aa_fasta}   = "./test/probes/korv_test5.faa"; 
	$loader_obj->{bitscore_min_tblastn}  = 100;
	$loader_obj->{seq_length_minimum}   = 50;
	my @test_targets;
	my $test1_target_path = "species/datatype/version/artificial_test1_korv.fa";
	push (@test_targets, $test1_target_path);
	$loader_obj->{target_paths} = \@test_targets;

	$digs_obj->{defragment_mode}  = 'defragment';
	$digs_obj->{defragment_range} = 100;

	my $initialise_obj = Initialise->new($digs_obj);
	$initialise_obj->setup_for_a_digs_run($digs_obj);
	$digs_obj->perform_digs();

	my @data;
	my @fields = qw [ assigned_gene assigned_name extract_start extract_end ];
	my $sort = " ORDER BY scaffold, extract_start ";
	$results_table->select_rows(\@fields, \@data, $sort);
	my $num_rows = scalar @data;

	my $fail = undef;
	if ($num_rows ne '5')  { 
		$fail = 1;
		die   "\n\t  tBLASTn screen, gag + env peptides: ** FAILED ($fail) ** Wrong number of rows in digs_results table \n";
	}
	my $result1_ref = shift @data;
	my $result2_ref = shift @data;
	my $result3_ref = shift @data;
	my $result4_ref = shift @data;
	my $result5_ref = shift @data;

	unless ($result1_ref->{assigned_gene} eq 'LTR'
	   and  $result1_ref->{assigned_name} eq 'KoRV'
	   and  $result1_ref->{extract_start} eq 200  
	   and  $result1_ref->{extract_end}   eq 703)    { $fail = 2; }
	unless ($result2_ref->{assigned_gene} eq 'gag'
	   and  $result2_ref->{assigned_name} eq 'KoRV'
	   and  $result2_ref->{extract_start} eq 4001  
	   and  $result2_ref->{extract_end}   eq 5566)   { $fail = 3; }
	unless ($result3_ref->{assigned_gene} eq 'pol'
	   and  $result3_ref->{assigned_name} eq 'KoRV'
	   and  $result3_ref->{extract_start} eq 5681  
	   and  $result3_ref->{extract_end}   eq 9064)   { $fail = 4; }
	unless ($result4_ref->{assigned_gene} eq 'env'
	   and  $result4_ref->{assigned_name} eq 'KoRV'
	   and  $result4_ref->{extract_start} eq 8946  
	   and  $result4_ref->{extract_end}   eq 10925)  { $fail = 5; }
	unless ($result5_ref->{assigned_gene} eq 'LTR'
	   and  $result5_ref->{assigned_name} eq 'KoRV'
	   and  $result5_ref->{extract_start} eq 10967  
	   and  $result5_ref->{extract_end}   eq 11470)  { $fail = 6; }
	unless ($fail) {
		print "\n\n\t  tBLASTn screen, gag + env peptides: ** PASSED **\n";
	}
	else {
		die   "\n\t  tBLASTn screen, gag + env peptides: ** FAILED ($fail) **\n";
	}
	sleep 2;
}

#***************************************************************************
# Subroutine:  run_test_6
# Description: Test consolidation function
#***************************************************************************
sub run_test_6 {

	my ($self) = @_;

	my $digs_obj = $self->{digs_obj};
	
	# Test consolidation
 	print "\n\n\t ### TEST 6: Testing consolidation function ~ + ~ + ~ \n";
	$digs_obj->{consolidate_range} = 200;

	my @test_targets;
	my $test1_target_path = "species/datatype/version/artificial_test1_korv.fa";
	push (@test_targets, $test1_target_path);
	my $loader_obj = $digs_obj->{loader_obj};
	$loader_obj->{target_paths} = \@test_targets;

	# Set target sequence files for screening
	my %targets;
	my $num_targets = $loader_obj->set_targets(\%targets);

	# Show error and exit if no targets found
	unless ($num_targets) {
		$loader_obj->show_no_targets_found_error();
	}

	# Set target sequence files for screening
	my %target_groups;
	$loader_obj->set_target_groups(\%targets, \%target_groups);
	$digs_obj->{target_groups} = \%target_groups; 
	$digs_obj->{genome_use_path} = './test/'; 

	my $initialise_obj = Initialise->new($digs_obj);
	$initialise_obj->set_up_consolidate_tables($digs_obj);
	$digs_obj->{defragment_mode} = 'consolidate';

	# Get contig lengths and capture in a table
	#$digs_obj->calculate_contig_lengths();	
	# Set the parameters for consolidation

	my %consolidate_settings;
	$consolidate_settings{range} = 200;
	$consolidate_settings{start} = 'extract_start';
	$consolidate_settings{end}   = 'extract_end';
	$digs_obj->{consolidate_settings} = \%consolidate_settings;

	# Do the consolidation
	my $num_consolidated = $digs_obj->consolidate_loci();
	my $fail = undef;
	if ($num_consolidated ne '2')  { 
		$fail = 1;
		die   "\n\t  Consolidation +ve orientation test: ** FAILED ($fail) ** Wrong number of rows\n";
	}
	
	unless ($fail) {
		print "\n\n\t  Consolidation +ve orientation test: ** PASSED **\n";
	}

	sleep 2;
}

#***************************************************************************
# Subroutine:  run_test_7
# Description: Reassign test 
#***************************************************************************
sub run_test_7 {

	my ($self) = @_;

	# Reassign test
	print "\n\n\t ### TEST 7: Reassign test ~ + ~ + ~ \n";

	# Upload a data set
	my $digs_obj   = $self->{digs_obj};
	my $loader_obj = $digs_obj->{loader_obj};
	my $db = $digs_obj->{db}; # Get the database reference
	my $results_table = $db->{digs_results_table}; # Get the database reference
	$results_table->flush();
	my $test7_path = './test/tabular/test7.txt';;
	print "\n\t ### Flushed digs_results table & now uploading data from file '$test7_path'";
	$db->upload_data_to_digs_results($test7_path);
	print "\n\t ### Data uploaded\n";

	# Set the parameters for this test (directly rather than using control file)
	$loader_obj->{seq_length_minimum}   = 50;
	$loader_obj->{defragment_range}     = 200;
	$loader_obj->{query_aa_fasta}       = "./test/probes/korv_test5.faa"; 
	$loader_obj->{reference_aa_fasta}   = "./test/references/korv_test7_reassign.faa"; 
	$loader_obj->{query_na_fasta}       = "./test/probes/korv_test1.fna"; 
	$loader_obj->{reference_na_fasta}   = "./test/probes/korv_test1.fna"; 
	$loader_obj->{bitscore_min_tblastn} = 100;

	my @test_targets;
	my $test1_target_path = "species/datatype/version/artificial_test1_korv.fa";
	push (@test_targets, $test1_target_path);
	$loader_obj->{target_paths} = \@test_targets;

	$digs_obj->{defragment_mode}  = 'defragment';
	$digs_obj->{defragment_range} = 100;
	$digs_obj->{force} = 'true';

	# If we're doing a reassign, get the assigned digs_results
	my @reassign_loci;
	$digs_obj->get_sorted_digs_results(\@reassign_loci);
	$digs_obj->{reassign_loci} = \@reassign_loci;

	# Set up the reference library
	my $initialise_obj = Initialise->new($digs_obj);
	$initialise_obj->setup_for_reassign($digs_obj);
	$digs_obj->reassign();	
	
	# Get data and check reassign looks right
	my @data;
	my @fields = qw [ assigned_gene assigned_name extract_start extract_end ];
	my $sort = " WHERE probe_type = 'ORF' ORDER BY scaffold, extract_start ";
	$results_table->select_rows(\@fields, \@data, $sort);
	#$devtools->print_array(\@data); die;
	my $fail = undef;
	foreach my $result_ref (@data) {
		unless ($result_ref->{assigned_name} eq 'MLV') { $fail = 1; }
	}
	
	# Final message
	unless ($fail) { print "\n\n\t  Test 7: Reassign test  ** PASSED **\n"; }
	else           { die   "\n\n\t  Test 7: Reassign test  ** FAILED ($fail) **\n"; }
		
	sleep 2;

}

#***************************************************************************
# Subroutine:  run_test_8
# Description: Reverse complement screen 
#***************************************************************************
sub run_test_8 {

	my ($self) = @_;

	# Test reverse complemente hit
	print "\n\t ### TEST 8: Reverse complement screen ~ + ~ + ~ \n";
	sleep 2;

}

#***************************************************************************
# Subroutine:  run_test_10
# Description:  big defragment with re-extract and genotype disabled
#***************************************************************************
sub run_test_10 {

	my ($self) = @_;

	# Run the second peptide screen
	print "\n\t ### TEST 10: Running big defragment with re-extract and genotype disabled ~ + ~ + ~ ";

	# Set to defragment
	$self->{defragment_mode} = 'defragment';

	# Get digs_results table handle and flush the table
	# Upload a data set
	my $digs_obj      = $self->{digs_obj};
	my $loader_obj    = $digs_obj->{loader_obj};
	my $db            = $digs_obj->{db}; # Get the database reference
	my $results_table = $db->{digs_results_table}; # Get the database reference
	$results_table->flush();

	# Upload the data
	my $test10_path = './test/tabular/test10.txt';;
	print "\n\t ### Flushed digs_results table & now uploading data from file '$test10_path'";
	my $original_num_rows = $db->upload_data_to_digs_results($test10_path);
	print "\n\t ### EVE data uploaded\n";
	die;
	
	$digs_obj->interactive_defragment();
	
	# Check result
	my @result;
	my @fields = qw [ assigned_gene assigned_name extract_start extract_end ];
	my $sort = " ORDER BY scaffold, extract_start ";
	$results_table->select_rows(\@fields, \@result, $sort);
	my $num_rows = scalar @result;
	print "\n\t ### Compressed from $original_num_rows to $num_rows rows\n\n\n";

	exit;
}

############################################################################
# EOF
############################################################################