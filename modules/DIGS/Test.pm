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
my $user     = 'root';
my $password = 'blenat2';
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
	
	# Load the 'digs_test' database
	$digs_obj->{mysql_server}   = $server;
	$digs_obj->{mysql_username} = $user;
	$digs_obj->{mysql_password} = $password;
	$digs_obj->initialise_screening_db('digs_test_screen');
	my $db = $digs_obj->{db}; # Get the database reference
	$db->flush_screening_db();

	# Set paths for processing etc
	#$config->{process_id} = 
	#$config->{report_dir} = 
	#$config->{tmp_path} = 
	#$config->{output_path} = 
	#$config->{blast_obj} = 

	# Display current settings	
	print "\n\n\t ### Running DIGS tests ~ + ~ + ~ \n";

	# Do a live screen using test control file and synthetic target data
	$self->run_test_1();
	#$self->run_test_2();
	#$self->run_test_3();
	#$self->run_test_4();
	#$self->run_test_5();
	#$self->run_test_6();
	#$self->run_test_7();
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

	print "\n\t ### TEST 1: Running live nucleotide screen against synthetic data ~ + ~ + ~ \n\n";
	my $digs_obj = $self->{digs_obj};

	# Do a DIGS run against synthetic data (included in repo)
	my %config = ScreenBuilder->new($digs_obj);
	#$config->{tmp_path} = 
	#$config->{seq_length_minimum} = 
	#$config->{consolidate_range} = 
	#$config->{bitscore_min_tblastn} = 
	#$config->{bitscore_min_blastn} = 
	#$config->{reference_aa_fasta} = 
	#$config->{defragment_range} = 
	#$config->{reference_aa_fasta} = 
	#$config->{query_aa_fasta} =  
	#$config->{genome_use_path} =  
	
	# Arrays
	#$config->{target_paths} =  
	#$config->{skipindexing_paths} =  
	#$config->{exclude_paths} =  

	# Connection details
	#db_name => erv_homo_rt
	#mysql_username => root
	#mysql_server => localhost
	#mysql_password => blenat2

	#previously_executed_searches

	$digs_obj->{loader_obj} = \%config;
	
	$digs_obj->setup_for_digs();
	die;
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
	if ($correct_result)  { print "\n\n\t  Live blastn screen test: ** PASSED **\n" }
	else                  { die   "\n\n\t  Live blastn screen test: ** FAILED **\n" }
	#$devtools->print_hash($result1_ref); $devtools->print_hash($result2_ref); die;
	sleep 1;
}

#***************************************************************************
# Subroutine:  run_test_2
# Description: Defragment results test negative (i.e. do not merge loci)
#***************************************************************************
sub run_test_2 {

	my ($self) = @_;

	## Check that defragment gives expected result (negative)
	print "\n\t ### TEST 2: Defragment results test negative  ~ + ~ + ~ \n";	
	$self->{defragment_mode} = 'defragment';
	# Construct WHERE statement
	my $where  = " WHERE organism      = 'fake_species' ";
	$where    .= " AND target_datatype = 'fake_datatype' ";
	$where    .= " AND target_version  = 'fake_version' ";
	$where    .= " AND target_name     = 'artificial_test1_korv.fa' "; 
	my $path         = "/test/fake_species/fake_datatype/fake_version/artificial_test1_korv.fa";
	my $target_path  = $ENV{DIGS_GENOMES}  . $path;
	my %settings;
	$settings{range}     = 100;
	$settings{start}     = 'extract_start';
	$settings{end}       = 'extract_end';
	$settings{where_sql} = $where;
	
	my $num_new = $self->defragment_target(\%settings, $target_path, 'digs_results');
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

	# Run a peptide screen
	print "\n\t ### TEST 3: Running live partially deleted pol peptide screen against synthetic data ~ + ~ + ~ \n";
	my $test_ctl_file2 = './test/test3_erv_aa.ctl';
	my $loader_obj     = $self->{loader_obj};
	$loader_obj->parse_control_file($test_ctl_file2, $self, 2);
	$self->setup_for_digs();
	$self->perform_digs();

	my $db = $self->{db}; # Get the database reference
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
	if ($correct_result)  { print "\n\n\t  Live tblastn test: ** PASSED **\n" }
	else                  { die   "\n\n\t  Live tblastn test: ** FAILED **\n" }
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

	# Set to defragment
	$self->{defragment_mode} = 'defragment';

	# Construct WHERE statement
	my $where  = " WHERE probe_type = 'ORF' ";
	my $path         = "/test/fake_species/fake_datatype/fake_version/artificial_test1_korv.fa";
	my $target_path  = $ENV{DIGS_GENOMES}  . $path;
	my %settings;
	$settings{range}     = 500;
	$settings{start}     = 'extract_start';
	$settings{end}       = 'extract_end';
	$settings{where_sql} = $where;

	my $num_new = $self->defragment_target(\%settings, $target_path, 'digs_results');
	my $db = $self->{db}; # Get the database reference
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
# Description: Live gag + env peptide screen
#***************************************************************************
sub run_test_5 {

	my ($self) = @_;

	# Run the second peptide screen
	print "\n\t ### TEST 5: Running live gag + env peptide screen (entails merge of result rows) ~ + ~ + ~ ";

	my $db = $self->{db}; # Get the database reference
	my $results_table = $db->{digs_results_table}; # Get the database reference
	$results_table->flush();
	my $test5_path = './test/test5.txt';;
	print "\n\t ### Flushed digs_results table & now uploading data from file '$test5_path'";
	$db->upload_data_to_digs_results($test5_path);
	print "\n\t ### Data uploaded, starting from point of having successfully conducted tests 1,2,3, & 4\n";
		
	my $test_ctl_file2 = './test/test5_erv_aa.ctl';
	my $loader_obj     = $self->{loader_obj};
	$loader_obj->parse_control_file($test_ctl_file2, $self, 2);
	$self->setup_for_digs();
	$self->perform_digs();

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

	# Test consolidation
	print "\n\n\t ### TEST 6: Testing consolidation function ~ + ~ + ~ \n";
	$self->{consolidate_range} = 200;
	my $num_consolidated = $self->consolidate_loci();
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
	my $db = $self->{db}; # Get the database reference
	my $results_table = $db->{digs_results_table}; # Get the database reference
	$results_table->flush();
	my $test7_path = './test/test7.txt';;
	print "\n\t ### Flushed digs_results table & now uploading data from file '$test7_path'";
	$db->upload_data_to_digs_results($test7_path);
	print "\n\t ### Data uploaded\n";

	# Defragment - Construct WHERE statement
	my $where  = " WHERE organism      = 'fake_species' ";
	$where    .= " AND target_datatype = 'fake_datatype' ";
	$where    .= " AND target_version  = 'fake_version' ";
	$where    .= " AND target_name     = 'artificial_test1_korv.fa' "; 
	my $path         = "/test/fake_species/fake_datatype/fake_version/artificial_test1_korv.fa";
	my $target_path  = $ENV{DIGS_GENOMES}  . $path;
	my %settings;
	$settings{range}     = 100;
	$settings{start}     = 'extract_start';
	$settings{end}       = 'extract_end';
	$settings{where_sql} = $where;

	$self->{defragment_mode} = 'defragment';
	$self->defragment_target(\%settings, $target_path, 'digs_results');

	# Do the reassign
	print "\n\n\t ### Now reassigning the ORF hits to an MLV reference library\n";
	my $test_ctl_file7 = './test/test7_erv_aa.ctl';
	my $loader_obj     = $self->{loader_obj};
	$loader_obj->parse_control_file($test_ctl_file7, $self, 2);
	my @digs_results;
	$self->initialise_reassign(\@digs_results); # Set up 
	$self->reassign(\@digs_results);	
	
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
	my $db = $self->{db}; # Get the database reference
	my $results_table = $db->{digs_results_table}; # Get the database reference
	$results_table->flush();

	# Upload the data
	my $test10_path = './test/test10.txt';;
	print "\n\t ### Flushed digs_results table & now uploading data from file '$test10_path'";
	my $original_num_rows = $db->upload_data_to_digs_results($test10_path);
	print "\n\t ### EVE data uploaded\n";

	# Execute the defragment procedure
	$self->interactive_defragment();
	
	# Check result
	my @result;
	my @fields = qw [ assigned_gene assigned_name extract_start extract_end ];
	my $sort = " ORDER BY scaffold, extract_start ";
	$results_table->select_rows(\@fields, \@result, $sort);
	my $num_rows = scalar @result;
	print "\n\t ### Compressed from $original_num_rows to $num_rows rows\n\n\n";

	exit;
}
