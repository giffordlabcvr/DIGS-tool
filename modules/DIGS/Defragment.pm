#!usr/bin/perl -w
############################################################################
# Module:      Defragment.pm
# Description: Functions for clustering, defragmenting, consolidating loci
# History:     April  2017: Created by Robert Gifford 
############################################################################
package Defragment;

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

############################################################################
# Globals
############################################################################

# Base objects
my $fileio    = FileIO->new();
my $console   = Console->new();
my $devtools  = DevTools->new();

# Ordered fields for digs result table files
my @ordered = qw [ digs_result_id organism target_datatype target_version target_name
			  probe_type scaffold extract_start extract_end sequence_length sequence
			  assigned_name assigned_gene orientation bitscore identity evalue_num evalue_exp
			  subject_start subject_end query_start query_end align_len gap_openings mismatches ];
			  
# Maximum range for defragment
my $max = 100000000;
my $start_token = 'extract_start';
my $end_token = 'extract_end';

1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create new Defragment 'object'
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Declare empty data structures
	my %crossmatching;

	# Set member variables
	my $self = {

		# Set-up params
		defragment_range       => $parameter_ref->{defragment_range}, 
		defragment_mode        => $parameter_ref->{defragment_mode}, 
	    defragment_settings    => $parameter_ref->{defragment_settings}, 
		
		# Member classes 
		db                     => $parameter_ref->{db},  
		blast_obj              => $parameter_ref->{blast_obj},
		   
		# Paths used in defragment
		genome_use_path        => $parameter_ref->{genome_use_path},
		output_path            => $parameter_ref->{output_path},

		# Paths used in re-extract
		target_groups          => $parameter_ref->{target_groups},

		# Paths used for re-assign
		tmp_path               => $parameter_ref->{tmp_path},
		aa_reference_library   => $parameter_ref->{aa_reference_library},
		na_reference_library   => $parameter_ref->{na_reference_library},
		blast_orf_lib_path     => $parameter_ref->{blast_orf_lib_path},
		blast_utr_lib_path     => $parameter_ref->{blast_utr_lib_path},

		# Flags
		verbose                => $parameter_ref->{verbose},
		force                  => $parameter_ref->{force},

	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# Compose clusters (this is the main defragment routine)
############################################################################

#***************************************************************************
# Subroutine:  compose_clusters 
# Description: process a sorted list of loci and group into 'clusters' of
#              overlapping feature annotations
#***************************************************************************
sub compose_clusters {

	my ($self, $defragmented_ref, $loci_ref) = @_;
	
	my $i = 1;
	my %cluster_tracking_data;
	my $initialised = undef;
	my $cluster_high_end = undef; # Variable to track highest end within a cluster of loci
	
	# Iterate through loci, grouping them into clusters when they are within range
	foreach my $locus_ref (@$loci_ref)  {

		my $id = $locus_ref->{digs_result_id};
		my $sc = $locus_ref->{scaffold};
		my $end = $locus_ref->{$end_token}; 
		
		# Get end variables for new hits 
		unless ($end) { 
			$end = $locus_ref->{subject_end}; 
		}
		
		unless ($initialised) { # Is this the first locus?
			$self->initialise_cluster($defragmented_ref, $locus_ref, $i);
			$cluster_tracking_data{scaffold} = $locus_ref->{scaffold};
			$self->set_tracking_data($locus_ref, \%cluster_tracking_data);
            $initialised = 'true'; # Set flag 
			$cluster_high_end = $end;
		}
		else {

			# Test if we need to merge this one
			my $merge = $self->is_distinct_locus($locus_ref, \%cluster_tracking_data, $cluster_high_end);
						
			unless ($merge) { # Distinct locus - initialise record
				$i++; # Increment the count
				$self->initialise_cluster($defragmented_ref, $locus_ref, $i);
				$cluster_high_end = $end;
			}
			else {             
				# Extend record
				$self->extend_cluster($defragmented_ref, $locus_ref, $i);
				if ($end > $cluster_high_end) {
					$cluster_high_end = $end;
				}
			}

			# Update tracking data
			$self->set_tracking_data($locus_ref, \%cluster_tracking_data);
		}	
	}
}

#***************************************************************************
# Subroutine:  is_distinct_locus
# Description: determine whether or not to merge a locus into previous
#***************************************************************************
sub is_distinct_locus {

	my ($self, $locus_ref, $data_ref, $high_end) = @_;

	my $merge = undef;
	
	# Is it on the same scaffold?
	my $same_scaffold = $self->is_locus_on_same_scaffold($locus_ref, $data_ref);

	# Is it a match to the same gene?
	my $same_gene = $self->is_locus_assigned_to_same_gene($locus_ref, $data_ref);

	# Is it in range?
	my $is_in_range = $self->is_locus_in_range($locus_ref, $high_end);
	
	if ($same_scaffold and $same_gene and $is_in_range) {
	
		# Should locus be merged based on current rules
		$merge = $self->should_locus_be_merged($locus_ref, $data_ref);
		$merge = 1;
	
	}
	
	return $merge;
}

#***************************************************************************
# Subroutine:  is_locus_on_same_scaffold
# Description: determine if a locus is on the same scaffold as another
#***************************************************************************
sub is_locus_on_same_scaffold {

	my ($self, $locus_ref, $data_ref) = @_;

	# Are loci on the same scaffold
	my $same_scaffold = undef;
	my $scaffold = $locus_ref->{scaffold};	
	my $last_scaffold = $data_ref->{scaffold};	
	if ($scaffold eq $last_scaffold) {
		$same_scaffold = 'true';
	}
	return $same_scaffold;
}

#***************************************************************************
# Subroutine:  is_locus_assigned_to_same_gene
# Description: determine if locus is assigned to same reference gene as another
#***************************************************************************
sub is_locus_assigned_to_same_gene {

	my ($self, $locus_ref, $data_ref) = @_;

	# Are loci on the same scaffold
	my $same_gene = undef;
	my $gene = $locus_ref->{assigned_gene};	
	my $last_gene = $data_ref->{assigned_gene};	
	if ($gene eq $last_gene) {
		$same_gene = 'true';
	}
	return $same_gene;
}

#***************************************************************************
# Subroutine:  is_locus_in_range
# Description: determine if a locus is merging range of another
#***************************************************************************
sub is_locus_in_range {

	my ($self, $locus_ref, $high_end) = @_;

	my $in_range = undef;
	
	# Get data structures and values
	my $range        = $self->{defragment_range};
	my $start        = $self->get_locus_start($locus_ref);
	my $end          = $self->get_locus_end($locus_ref);

	my $buffered_start = $start - $range;
	unless ($range) { die; }	
	#print "\n\t\t\t Check if '$buffered_start' < '$high_end'";
	if ($buffered_start < $high_end) {
		#print "\n\t\t\t MERGE cos '$buffered_start' < '$high_end'";
		$in_range = 'true';
	}
	
	return $in_range;
}

#***************************************************************************
# Subroutine:  should_locus_be_merged
# Description: decide if an in-range locus is merged in current process
#***************************************************************************
sub should_locus_be_merged {

	my ($self, $locus_ref, $data_ref) = @_;

	# Get the values
	my $mode             = $self->{defragment_mode};
	my $name             = $locus_ref->{assigned_name};
	my $gene             = $locus_ref->{assigned_gene};					
	my $orientation      = $locus_ref->{orientation};					
	my $last_name        = $data_ref->{assigned_name};
	my $last_gene        = $data_ref->{assigned_gene};		
	my $last_orientation = $data_ref->{orientation};			

	# Check orientation
	if ($orientation ne $last_orientation) {
		unless ($mode eq 'consolidate') { 
			print "\n\t\t Identified pair of loci that are in range, but different orientations";
			return 0;
		}
		else {
			unless ($last_gene and $gene) { die; } # Should never get here		
		}
	}

	unless ($gene) { # If there is no assigned gene set, use probe gene	
		$gene = $locus_ref->{probe_gene};
	}			
	unless ($last_gene) {  # If there is no assigned gene set, use probe gene	
		$last_gene = $data_ref->{probe_gene};
	}
	
	# Take action depending on whether we are DEFRAGMENTING or CONSOLIDATING
	if ($mode eq 'defragment') {
		if ($gene ne $last_gene) { return 0; }  # different genes
	}
	elsif ($mode eq 'consolidate' or $mode eq 'consolidate2') { 
		# do nothing (these loci can be merged, even though different genes)
	}
	else { # Shouldn't get here
		die;
	}
	return 1;	
}

#***************************************************************************
# Subroutine:  set tracking data
# Description: 
#***************************************************************************
sub set_tracking_data {

	my ($self, $locus_ref, $cluster_tracking_data_ref) = @_;

	# Get locus data
	my $scaffold      = $locus_ref->{scaffold};
	my $target_name   = $locus_ref->{target_name};
	my $assigned_name = $locus_ref->{assigned_name};
	my $probe_gene    = $locus_ref->{probe_gene}; # Currently needed for unassigned 
	my $assigned_gene = $locus_ref->{assigned_gene};
	my $orientation   = $locus_ref->{orientation};
	my $start         = $locus_ref->{$start_token};
	my $end           = $locus_ref->{$end_token};
		
	# Get values required in this function 
	my $last_scaffold      = $cluster_tracking_data_ref->{scaffold};
	my $last_start         = $cluster_tracking_data_ref->{$start_token};
	
	# Sanity checking - are sequences in sorted order for this scaffold?
	#if ( $scaffold eq $last_scaffold and $start < $last_start) {
	#	#$devtools->print_hash($locus_ref);
	#	my $error = "\n\t Error: sequences do not appear to be correctly sorted"; 
	#	$error   .= "\n\t (i.e. by unique scaffold ID and locus start position)\n\n"; 
	#	die $error; 
	#}

	# Update cluster tracking data	
	$cluster_tracking_data_ref->{assigned_name} = $assigned_name;
	$cluster_tracking_data_ref->{assigned_gene} = $assigned_gene;
	$cluster_tracking_data_ref->{probe_gene}    = $probe_gene;
	$cluster_tracking_data_ref->{scaffold}      = $scaffold;
	$cluster_tracking_data_ref->{orientation}   = $orientation;
	$cluster_tracking_data_ref->{$start_token}  = $start;
	$cluster_tracking_data_ref->{$end_token}    = $end;

}

#***************************************************************************
# Subroutine:  reextract_and_reassign
# Description: 
#***************************************************************************
sub reextract_and_reassign {

	my ($self, $target_ref, $loci_ref) = @_;

	my $genome_use_path  = $self->{genome_use_path};
	my $target_group_ref = $self->{target_groups};
	my $copy_table_name  = $self->{copy_table_name};
	my $extract_obj      = Extract->new($self);
	my $db_ref           = $self->{db};

	unless ($genome_use_path and $target_group_ref)  { die; }

	# Defragment each target file in turn
	# Get the target details (and thus the target path)
	#my $target_path = $self->get_target_file_path($target_ref);
	my $organism        = $target_ref->{organism};
	my $target_datatype = $target_ref->{target_datatype};
	my $target_version  = $target_ref->{target_version};
	my $target_name     = $target_ref->{target_name};
	my @genome = ( $organism , $target_datatype, $target_version, $target_name );
	my $target_id       = join ('|', @genome);
	print "\n\t\t # Defragmenting hits in '$target_name'";

	# If we're re-extracting, get the path
	my $target_group = $target_group_ref->{$target_id};
	unless ($target_group) { 
		$devtools->print_hash($target_group_ref);
		print "\n\t Didn't get target group name for target file with id '$target_id'\n\n";
		die; 
	}
		
	# Construct the path to this target file
	my @path;
	push (@path, $genome_use_path);
	push (@path, $target_group);
	push (@path, $organism);
	push (@path, $target_datatype);
	push (@path, $target_version);
	push (@path, $target_name);
	my $target_path = join ('/', @path);

	# Extract newly identified or extended sequences
	my @extracted;
	$extract_obj->extract_sequences_using_blast($target_path, $loci_ref, \@extracted);	
			
	# Do the genotyping step for the newly extracted locus sequences
	my $assigned_count   = 0;
	my $crossmatch_count = 0;
	my $num_extracted = scalar @extracted;

	print "\n\t\t # Genotyping $num_extracted newly extracted sequences:";
	foreach my $hit_ref (@extracted) { # Iterate through loci	
	
		my $classify_obj = Classify->new($self);
		$classify_obj->classify_sequence_using_blast($hit_ref);
		$assigned_count++;
		my $remainder = $assigned_count % 100;
		if ($remainder eq 0) { print "\n\t\t\t # $assigned_count sequences classified "; }
		
	}
	print "\n\t\t\t # $assigned_count sequences classified";

	# Update DB
	my $num_deleted = $db_ref->update_db(\@extracted, $copy_table_name, 1);
	print "\n\t\t\t # $num_deleted rows deleted from digs_results table\n";

}

#***************************************************************************
# Subroutine:  initialise_cluster
# Description: create the first element in a cluster of associated loci
#***************************************************************************
sub initialise_cluster {

	my ($self, $defragmented_ref, $hit_ref, $count) = @_;

    # Get the current hit values
	#print "\n\t New record ($count)";
    my @array;
    my %hit = %$hit_ref;
    push (@array, \%hit);
    $defragmented_ref->{$count} = \@array;
}

#***************************************************************************
# Subroutine:  extend_cluster 
# Description: add a new element to an initialised cluster of associated loci
#***************************************************************************
sub extend_cluster {

	my ($self, $defragmented_ref, $hit_ref, $count) = @_;

    # Get the current hit values
	#print "\n\t Extending record ($count) ";
    #$devtools->print_hash($hit_ref); die;
    my $array_ref = $defragmented_ref->{$count};
    push (@$array_ref, $hit_ref);

}

############################################################################
# MERGE clusters
############################################################################

#***************************************************************************
# Subroutine:  merge_clustered_loci
# Description: resolve each of multiple clustered loci to one locus
#***************************************************************************
sub merge_clustered_loci {
	
	my ($self, $clusters_ref, $merged_ref, $singletons_ref, $delete_ref) = @_;

	# Iterate through the clusters
	my @cluster_ids  = keys %$clusters_ref;
	foreach my $id (@cluster_ids) {

		# Get data for this cluster
		my $cluster_ref = $clusters_ref->{$id};
		my $num_loci = scalar @$cluster_ref;

		if ($num_loci eq 1) {	
			$num_loci++;
			my $locus_ref = shift @$cluster_ref;
			$singletons_ref->{$id} = $locus_ref;
			#$devtools->print_array($cluster_ref);
		}
		else {		
			$self->merge_cluster($id, $cluster_ref, $merged_ref, $delete_ref);
		}
	}	
}

#***************************************************************************
# Subroutine:  merge_cluster
# Description: resolve a cluster of overlapping loci to one single locus
#***************************************************************************
sub merge_cluster {
	
	my ($self, $cluster_id, $cluster_ref, $merged_ref, $delete_ref) = @_;
	
	unless ($cluster_id and $cluster_ref) { die; }

	# Create the data structures we need
	my %previously_extracted_locus_ids;
	my $highest_end   = undef;
	my $lowest_start  = undef;
	my $previous_digs_result_id = undef;
	foreach my $locus_ref (@$cluster_ref) {
 													
		# Work out if we need to delete this row when we update
		my $id  = $locus_ref->{digs_result_id};					
		if ($id) {
			push(@$delete_ref, $id);
		}		

		# Record the start and stop parameters so we know whether or not to extend
		my $start = $self->get_locus_start($locus_ref);			
		my $end   = $self->get_locus_end($locus_ref);
					
		# Get start and end variables for new hits
		if ($lowest_start and $highest_end ) {		
			if ($start < $lowest_start) { $lowest_start = $start; }
			if ($end > $highest_end)    { $highest_end  = $end;   }
		}
		elsif ($start and $end) {
			$highest_end  = $end;
			$lowest_start = $start;
		}
		else { die; }
	
 	}

	# Create the merged copy
	my $template_ref = @$cluster_ref[0];
	#$devtools->print_hash($template_ref); die;
	my %data = %$template_ref;
	$data{extract_start} = $lowest_start;
	$data{extract_end} = $highest_end;
	$data{sequence_length} = ($highest_end - $lowest_start) + 1;
	$merged_ref->{$cluster_id} = \%data;
	
}

#***************************************************************************
# Subroutine:  get_loci_to_extract
# Description: 
#***************************************************************************
sub get_loci_to_extract {

	my ($self, $merged_ref, $singletons_ref, $to_extract_ref) = @_;

	my $verbose = $self->{verbose};

	# Identify singletons to extract
	my @s_keys = keys %$singletons_ref;
	my @singletons;
	foreach my $key (@s_keys) {	
		
		# Extract if its from a new search
		my $locus_ref = $singletons_ref->{$key};
		my $digs_result_id = $locus_ref->{digs_result_id};

		unless ($digs_result_id) {  # If "digs_result_id" is undefined its new
			my $start = $locus_ref->{subject_start};
			my $end   = $locus_ref->{subject_end};
			$locus_ref->{start} = $start;
			$locus_ref->{end} = $end;
			push (@$to_extract_ref, $locus_ref);
		}		
	}

	# Convert merged to tabular data
	my @m_keys = keys %$merged_ref;
	my @merged;
	foreach my $key (@m_keys) {

		# If its been merged it definitely needs to be extracted
		my $locus_ref = $merged_ref->{$key};
		my $start = $locus_ref->{extract_start};
		my $end   = $locus_ref->{extract_end};
		unless ($start and $end ) { die; } # Should definitely have these
		$locus_ref->{start} = $start;
		$locus_ref->{end}   = $end;
		if ($verbose) { print "\n\t\t    - Adding locus '$start - $end'"; }
		
		push (@$to_extract_ref, $locus_ref);

	}
}

#***************************************************************************
# Subroutine:  convert_locus_hash_to_locus_line
# Description: does what it says, pretty straightforward
#***************************************************************************
sub convert_locus_hash_to_locus_line {
	
	my ($self, $locus_ref, $separator) = @_;

	unless ($separator) {
		$separator = ',';
	}

	my @line;
	foreach my $field (@ordered) {
	
		my $value = $locus_ref->{$field};
		unless ($value) { $value = '0'; }
		push (@line, $value);	
	}
	my $line = join($separator, @line);
	$line .= "\n";

	return $line;
}

############################################################################
# INTERACTIVE DEFRAGMENT
############################################################################

#***************************************************************************
# Subroutine:  interactive_defragment 
# Description: interactively defragment via the console
#***************************************************************************
sub interactive_defragment {

	my ($self) = @_;

	# Get objects
	my $db        = $self->{db};
	my $dbh       = $db->{dbh};
	$db->load_digs_results_table($dbh, 'digs_results');

	# Get the sorted list of results
	my @loci;
	$db->get_sorted_digs_results(\@loci);
		
	my $choice = 1;		
	do { 
	
		# Display current config and set the range	
		$self->display_config();

		my $current_range  = $self->{defragment_range};
		my $range_question = "\n\n\t # Set the range for merging hits";
		my $t_range = $console->ask_int_with_bounds_question($range_question, $current_range, $max);		
		$self->{defragment_range} = $t_range;

	    print "\n\t # Previewing defragment result using range '$t_range'\n";
	    my $preview = 'true';
	    $self->defragment(\@loci, $preview);
	    
		print "\n\n\t\t Option 1: preview new parameters";
		print "\n\t\t Option 2: apply these parameters";
		print "\n\t\t Option 3: exit";
		my $list_question = "\n\n\t # Choose an option:";
		$choice = $console->ask_list_question($list_question, 3);
		
	} until ($choice > 1);
	if ($choice eq 2) {

		print "\n\t # Applying parameters...\n ";
	    $self->defragment(\@loci);
	}	
	elsif ($choice eq 3) {  # Exit
		print "\n"; exit;
	}
	
	else { die; } # Should never get here
}


#***************************************************************************
# Subroutine:  defragment_target_files
# Description: identify overlapping or contiguous loci in digs_results tbl
#***************************************************************************
sub defragment {

	my ($self, $loci_ref, $preview) = @_;

	# Create the various objects we will use here
	my $extract_obj  = Extract->new($self);
	my $classify_obj = Classify->new($self);	# Get objects
	my $db_ref       = $self->{db};
	my $verbose      = $self->{verbose};

	# Identify clusters
	my %clusters;
	$self->compose_clusters(\%clusters, $loci_ref);
	if ($verbose) {
		$self->show_clusters(\%clusters);
	}
	
	# Derive a non-redundant locus set (i.e. one merged locus for each cluster)
	my %merged;
	my %singletons;
	my @to_delete;
	$self->merge_clustered_loci(\%clusters, \%merged, \%singletons, \@to_delete);
    # DEBUG $devtools->print_array(\@to_delete); # die;
		
	if ($preview) {
	
		$self->display_config();
		$self->summarise_defragment_result($loci_ref, \%merged, \%singletons, \@to_delete);
		return;
	}
	
	# Get loci to extract
	my @to_extract;
	$self->get_loci_to_extract(\%merged, \%singletons, \@to_extract);
    # DEBUG $devtools->print_array(\@to_extract); # die;

	# Extract newly identified or extended sequences
	my $target_groups_ref = $self->{target_groups};
    # DEBUG $devtools->print_hash($target_groups_ref); # die;
  	my %by_target;
	my @extracted;

	foreach my $locus_ref (@to_extract) {
  	    
		# Extract sequence using BLAST
		my $organism    = $locus_ref->{organism};
		my $type        = $locus_ref->{target_datatype};
		my $version     = $locus_ref->{target_version};
		my $file        = $locus_ref->{target_name};
		my @path_parts  = ( $organism , $type, $version, $file );
		my $target_key = join ('|', @path_parts); 
		my $stem_group  = $target_groups_ref->{$target_key};
		unless ($stem_group) { die; }
		my $genome_path = join ('/', @path_parts); 
		my $use_path = $self->{genome_use_path};
		unless ($use_path) { die; }
		my $target_path  = $use_path . '/' . $stem_group . '/' . $genome_path;
			
		$locus_ref->{target_path} = $target_path;
		$extract_obj->extract_locus_sequence_using_blast($locus_ref, \@extracted);	

		# Classify using BLAST
		$classify_obj->classify_sequence_using_blast($locus_ref);	
		my %copy_locus = %$locus_ref;
		push (@extracted, \%copy_locus);
	}
	
	# Update the digs_results table	
	$db_ref->update_db(\@to_delete, \@extracted, 'digs_results_table');

}

#***************************************************************************
# Subroutine:  write_defragment_files
# Description: write output from the defragment process to tabular files
#***************************************************************************
sub write_defragment_files {

	my ($self, $merged_ref, $singletons_ref) = @_;

	# Convert singletons to tabular data
	my @s_keys = keys %$singletons_ref;
	my @singletons;
	my $header = join(',', @ordered);
	$header .= "\n";
	push (@singletons, $header);
	foreach my $key (@s_keys) {
	
		my $singleton_ref = $singletons_ref->{$key};
		my $line = $self->convert_locus_hash_to_locus_line($singleton_ref);
		push (@singletons, $line);
	}

	# Convert merged to tabular data
	my @m_keys = keys %$merged_ref;
	my @merged;
	push (@merged, $header);
	foreach my $key (@m_keys) {
	
		my $merged_ref = $merged_ref->{$key};
		my $line = $self->convert_locus_hash_to_locus_line($merged_ref);
		push (@merged, $line);
	}

	# Create file with merged loci
	my $output_path = $self->{output_path};
	my $merged_file = "$output_path/merged_file.txt";
	$fileio->write_file($merged_file, \@merged);
	my $singletons_file = "$output_path/singletons_file.txt";
	$fileio->write_file($singletons_file, \@singletons);
	my @combined = (@singletons, @merged);
	unshift (@singletons, $header);
	my $combined = "$output_path/combined_all_file.txt";
	$fileio->write_file($combined, \@combined);
	#print "\n\t FILES WRITTEN!\n\n";

}

#***************************************************************************
# Subroutine:  get_locus_start
# Description: does what it says, pretty straightforward
#***************************************************************************
sub get_locus_start {
	
	my ($self, $locus_ref) = @_;

	#  Get start and end variables for new hits
 	my $start  = $locus_ref->{$start_token};
	unless ($start) { 
		$start = $locus_ref->{subject_start}; 
	}
	unless ($start) {  die; }		
	return $start;
}

#***************************************************************************
# Subroutine:  get_locus_end
# Description: does what it says, pretty straightforward
#***************************************************************************
sub get_locus_end {
	
	my ($self, $locus_ref) = @_;

 	my $end = $locus_ref->{$end_token};
	unless ($end) { 
		$end = $locus_ref->{subject_end}; 
	}
	unless ($end) {  die; }	
	return $end;
}

############################################################################
# CONSOLE OUTPUT 
############################################################################

#***************************************************************************
# Subroutine:  show_clusters
# Description: print information about clustered loci to the console
#***************************************************************************
sub show_clusters {

	my ($self, $defragmented_ref) = @_;

	my @cluster_ids = keys %$defragmented_ref;
	my $cluster_count;
	foreach my $id (@cluster_ids) {
		$cluster_count++;
		my $cluster_ref = $defragmented_ref->{$id};
		my $cluster_size = scalar @$cluster_ref;
		#print "\n\t CLUSTER '$id' SIZE = '$cluster_size'";
		if ($cluster_size > 1) {
			$self->show_cluster($cluster_ref, $cluster_count);
		}
	}
	print "\n";

}

#***************************************************************************
# Subroutine:  show_cluster
# Description: print information about a cluster of loci to the screen
#***************************************************************************
sub show_cluster {

	my ($self, $cluster_ref, $cluster_id) = @_;

 	# Iterate through the loci in this cluster
	print "\n\t CLUSTER '$cluster_id'";
	my $i;
	foreach my $locus_ref (@$cluster_ref) {
   		
   		$i++;		
		my $organism      = $locus_ref->{organism};
		my $assigned_name = $locus_ref->{probe_name};
		my $assigned_gene = $locus_ref->{probe_gene};
		my $orientation   = $locus_ref->{orientation};
		my $track_name    = $locus_ref->{track_name};
		unless ($track_name) {
			$track_name  = $locus_ref->{source};
		}
		unless ($assigned_name) {
			$assigned_name = $locus_ref->{assigned_name};
		}
		unless ($assigned_gene) {
			$assigned_gene = $locus_ref->{assigned_gene};
		}

		my $scaffold      = $locus_ref->{scaffold};
		
		# Record the start and stop parameters so we know whether or not to extend
		my $start = $self->get_locus_start($locus_ref);			
		my $end   = $self->get_locus_end($locus_ref);

		print "\n\t\t Scaffold '$scaffold' CLUSTER locus $i: $assigned_name: $assigned_gene:  $start-$end ($orientation)";
	}			
}

#***************************************************************************
# Subroutine:  display_config
# Description: display current defragment configuration
#***************************************************************************
sub display_config {

	my ($self, $t_range) = @_;

	my $current_range = $self->{defragment_range};
	unless ($current_range )  { die; } 
	print "\n\n\t # Current settings:\n";
	print "\n\t\t   defragment range:      $current_range";
	if ($t_range) {
		print "\n\t\t   current preview range: $t_range";
	}
}

#***************************************************************************
# Subroutine:  summarise_defragment_result
# Description: summarise the extent to which digs_result table was compressed
#***************************************************************************
sub summarise_defragment_result {

	my ($self, $loci_ref, $merged_ref, $singletons_ref, $to_delete_ref) = @_;

	my $total_original= scalar @$loci_ref;
	my $total_merged  = scalar keys %$merged_ref;
	my $total_singles = scalar keys %$singletons_ref;
	my $total_delete  = scalar @$to_delete_ref;
	my $total_loci    = $total_merged + $total_singles;
	print "\n\n\t\t   TOTAL ORIGINAL:      $total_original";
	print "\n\t\t   TOTAL MERGED:        $total_merged ";
	print "\n\t\t   TOTAL SINGLES:       $total_singles";
	print "\n\t\t   TOTAL AFTER MERGE:   $total_loci ";
	
}

############################################################################
# EOF
############################################################################