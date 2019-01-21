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

# Maximum range for defragment
my $max = 100000000;
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

	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# DEFRAGMENT FUNCTIONS
############################################################################

#***************************************************************************
# Subroutine:  interactive_defragment 
# Description: interactively defragment via the console
#***************************************************************************
sub interactive_defragment {

	my ($self) = @_;

	my $choice = undef;		
	do { 
		$choice = $self->preview_defragment();
	} until ($choice > 1);
	if ($choice eq 2) {
		$self->defragment(); # Apply the changes
	}	
	elsif ($choice eq 3) {  # Exit
		print "\n"; exit;
	}
	
	else { die; } # Should never get here
}

#***************************************************************************
# Subroutine:  preview_defragment
# Description: preview a defragmentation process
#***************************************************************************
sub preview_defragment {

	my ($self) = @_;

	# Display current confi and set the range	
	$self->display_config();
	my $current_range  = $self->{defragment_range};
	my $settings_ref   = $self->{defragment_settings};
	my $targets_ref    = $settings_ref->{targets};
	my $range_question = "\n\n\t # Set the range for merging hits";
	#my $t_range = $console->ask_int_with_bounds_question($range_question, $current_range, $max);		
	my $t_range = 1050;
	$settings_ref->{range} = $t_range;

	# Defragment this set of DIGS results
	print "\n\t # Previewing defragment result using range '$t_range'\n";
	my $preview = 'true';
	my $reextract = undef;
	$self->defragment($preview, $reextract);

	# Display result and prompt for what to do next
	$self->display_config();
	$self->summarise_defragment_result();
	
	print "\n\n\t\t Option 1: preview new parameters";
	print "\n\t\t Option 2: apply these parameters";
	print "\n\t\t Option 3: exit";
	my $list_question = "\n\n\t # Choose an option:";
	my $choice = $console->ask_list_question($list_question, 3);
	return $choice;

}

#***************************************************************************
# Subroutine:  defragment_target_files
# Description: identify overlapping or contiguous loci in digs_results tbl
#***************************************************************************
sub defragment {

	my ($self, $preview, $reextract) = @_;

	# Get objects
	my $db        = $self->{db};
	my $dbh       = $db->{dbh};
	$db->load_digs_results_table($dbh, 'digs_results');

	# Get the sorted list of results
	my @loci;
	$db->get_sorted_digs_results(\@loci);

	# Get the sorted results
	my %clusters;
	$self->compose_clusters(\%clusters, \@loci);
	$self->show_clusters(\%clusters);

	# Derive a single, merged locus for clusters of loci that were grouped into a cluster
	my %merged;
	my @to_reextract;
	$self->merge_clustered_loci(\%clusters, \%merged, @to_reextract);

	# Extract & reassign
	if ($reextract) {
		#$self->reextract_and_reassign($target_ref, \%merged);
	}	
	else {  # Update DB	

		# Create file with merged loci
		my $merged_file = 'merged_file.txt';
		$self->write_merged_locus_file($merged_file, \@loci);
		print "\n\t FILE WRITTEN!\n\n";
		exit;
		
		#my $db = $self->{db};
		#$db->import_data_to_digs_results($merged_file);
	}		
}

############################################################################
# INTERNAL FUNCTIONS: base fxns for merging overlapping/adjacent loci
############################################################################

#***************************************************************************
# Subroutine:  compose_clusters 
# Description: process a sorted list of loci and group into 'clusters' of
#              overlapping feature annotations
# Argument:    
#***************************************************************************
sub compose_clusters {

	my ($self, $defragmented_ref, $loci_ref) = @_;
	
	my $i = 1;
	my %cluster_tracking_data;
	my $initialised = undef;
	my $cluster_high_end = undef; # Variable to track highest end within a cluster of loci
	my $settings_ref = $self->{defragment_settings};
	my $end_token = $settings_ref->{end};
	
	# Iterate through loci, grouping them into clusters when they are within range
	foreach my $locus_ref (@$loci_ref)  {

		my $id = $locus_ref->{digs_result_id};
		my $sc = $locus_ref->{scaffold};
		my $end = $locus_ref->{$end_token}; 
		
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
# Subroutine:  merge_clustered_loci
# Description: resolve each of multiple clustered loci to one locus
#***************************************************************************
sub merge_clustered_loci {
	
	my ($self, $clusters_ref, $merged_ref, $to_extract_ref) = @_;

	# Iterate through the clusters
	my $extend_count = '0';
	my @cluster_ids  = keys %$clusters_ref;
	foreach my $id (@cluster_ids) {

		# Get data for this cluster
		my $cluster_ref = $clusters_ref->{$id};
		unless ($cluster_ref) { die; }
		my $extended = $self->merge_cluster($id, $cluster_ref, $to_extract_ref);
		if ($extended) { $extend_count = $extend_count + $extended; }
		
	}	
	return $extend_count;	
}

#***************************************************************************
# Subroutine:  merge_cluster
# Description: resolve a cluster of overlapping loci to one single locus
#***************************************************************************
sub merge_cluster {
	
	my ($self, $cluster_id, $cluster_ref, $to_extract_ref) = @_;
	
	my $num_loci = '0';
	foreach my $locus_ref (@$cluster_ref) {
		
		$num_loci++;			
							
		# Check if this is a a previously extracted locus

		# Store this BLAST result as part of a chain if it is new

		# Record the start and stop parameters so we know whether or not to extend

	}

	# Determine whether or not we need to extract sequences for this cluster
	# Extract if cluster is composed entirely of new loci (no previous result IDs)
 	
	# Is this a merge of multiple previously extracted loci?

	# If it includes a single extracted locus, does this locus need to be extended?	

	# If the locus needs to be extracted record the details

}


#***************************************************************************
# Subroutine:  is_distinct_locus
# Description: 
#***************************************************************************
sub is_distinct_locus {


}

#***************************************************************************
# Subroutine:  is_distinct_locus
# Description: 
#***************************************************************************
sub extract_locus {


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

	# Is it in range?
	my $is_in_range;
	if ($same_scaffold) {
		$is_in_range = $self->is_locus_in_range($locus_ref, $high_end);
	}
	if ($same_scaffold and $is_in_range) {
	
		# Should locus be merged based on current rules
		#$merge = $self->should_locus_be_merged($locus_ref, $data_ref);
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
# Subroutine:  is_locus_in_range
# Description: determine if a locus is merging range of another
#***************************************************************************
sub is_locus_in_range {

	my ($self, $locus_ref, $high_end) = @_;

	my $in_range = undef;
	
	# Get data structures and values
	my $settings_ref = $self->{defragment_settings};
	my $range        = $settings_ref->{range};
	my $start_token  = $settings_ref->{start};
	my $end_token    = $settings_ref->{end};
	my $start        = $locus_ref->{$start_token};
	my $end          = $locus_ref->{$end_token};
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

	# Get settings
	my $settings_ref = $self->{defragment_settings};
	my $start_token = $settings_ref->{start};
	my $end_token   = $settings_ref->{end};
	
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
# Subroutine:  write_merged_locus_file
# Description: write a set of merged loci to a file (can be imported to DB)
# Arguments: loci_ref: empty array to store the merged loci
#***************************************************************************
sub write_merged_locus_file {

    my ($self, $merged_file, $loci_ref) = @_;

	unless ($merged_file and $loci_ref) { die; }
	#$devtools->print_array($loci_ref) ; die;
	
	my @merged;
	foreach my $locus_row_ref (@$loci_ref) {

		my @line;
		#$devtools->print_hash($locus_row_ref);

		push (@line, $locus_row_ref->{record_id});
		push (@line ,$locus_row_ref->{organism});
		push (@line ,$locus_row_ref->{target_datatype});
		push (@line ,$locus_row_ref->{target_version});
		push (@line ,$locus_row_ref->{target_name});
		push (@line ,$locus_row_ref->{probe_type});
		push (@line ,$locus_row_ref->{scaffold});
		push (@line ,$locus_row_ref->{extract_start});
		push (@line ,$locus_row_ref->{extract_end});
		push (@line ,$locus_row_ref->{sequence_length});
		push (@line ,$locus_row_ref->{sequence});
		push (@line ,$locus_row_ref->{assigned_name});
		push (@line ,$locus_row_ref->{assigned_gene});
		push (@line ,$locus_row_ref->{orientation});
		push (@line ,$locus_row_ref->{bitscore});
		push (@line ,$locus_row_ref->{identity});
		push (@line ,$locus_row_ref->{evalue_num});
		push (@line ,$locus_row_ref->{evalue_exp});
		push (@line ,$locus_row_ref->{subject_start});
		push (@line ,$locus_row_ref->{subject_end});
		push (@line ,$locus_row_ref->{query_start});
		push (@line ,$locus_row_ref->{query_end});
		push (@line ,$locus_row_ref->{align_len});
		push (@line ,$locus_row_ref->{gap_openings});
		push (@line ,$locus_row_ref->{mismatches});

		# Deal with empty rows 
		foreach my $value (@line) {		
			unless ($value) {
				$value = '0';
			} 
		}

		my $line = join("\t", @line);
		push (@merged, "$line\n");		

	}
	
	$fileio->write_file($merged_file, \@merged);
}

#***************************************************************************
# Subroutine:  reextract_and_reassign
# Description: 
#***************************************************************************
sub reextract_and_reassign {

	my ($self, $target_ref, $loci_ref) = @_;

	my $genome_use_path  = $self->{genome_use_path};
	my $target_group_ref = $self->{target_groups};
	my $settings_ref     = $self->{defragment_settings};
	my $copy_table_name  = $settings_ref->{table_name};
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
		print "\n\t CLUSTER '$id' SIZE = '$cluster_size'";
		if ($cluster_size > 1) {
			$self->show_cluster($cluster_ref, $cluster_count);
		}
	}
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
		my $start         = $locus_ref->{extract_start};
		my $end           = $locus_ref->{extract_end};
		unless ($start)   { $start = $locus_ref->{subject_start}; }
		unless ($end)     { $end   = $locus_ref->{subject_end};	  }
		my $digs_result_id    = $locus_ref->{digs_result_id};

		print "\n\t\t CLUSTER locus $i: $assigned_name: $assigned_gene: $scaffold $start-$end ($orientation)";
	}			
}

#***************************************************************************
# Subroutine:  display_config
# Description: display current defragment configuration
#***************************************************************************
sub display_config {

	my ($self) = @_;

	my $settings_ref = $self->{defragment_settings};
	my $t_range      = $settings_ref->{range};

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

	my ($self) = @_;

	my $settings_ref   = $self->{defragment_settings};
	my $total_loci     = $settings_ref->{total_loci};
	my $total_clusters = $settings_ref->{total_clusters};
	print "\n\t\t   TOTAL LOCI:            $total_loci";
	print "\n\t\t   TOTAL CLUSTERS:        $total_clusters ";
	
}

############################################################################
# EOF
############################################################################
