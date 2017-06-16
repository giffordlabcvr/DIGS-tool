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
	    consolidate_settings   => $parameter_ref->{consolidate_settings}, 
		
		# Member classes 
		db                     => $parameter_ref->{db},  
		blast_obj              => $parameter_ref->{blast_obj},
		   
		# Paths used in defragment & consolidate processes
		genome_use_path        => $parameter_ref->{genome_use_path},
		output_path            => $parameter_ref->{output_path},
		target_groups          => $parameter_ref->{target_groups}

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

	# Declarations etc
	my $choice       = undef;
	my $settings_ref = $self->{defragment_settings};
 
	do { # Preview changes 
	
		$choice = $self->preview_defragment();
		
	} until ($choice > 1);

	if ($choice eq 2) {  # Apply the changes
		
		$self->defragment();
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

	# Display current settings	
	$self->display_config();

	# Get the range	
	my $current_range  = $self->{defragment_range};
	my $settings_ref   = $self->{defragment_settings};
	my $targets_ref    = $settings_ref->{targets};
	my $range_question = "\n\n\t # Set the range for merging hits";
	my $t_range = $console->ask_int_with_bounds_question($range_question, $current_range, $max);		

	my $settings_ref = $self->{defragment_settings};
	$settings_ref->{range} = $t_range;

	# Defragment this set of DIGS results
	print "\n\t # Previewing defragment'\n";
	# Summarise the results of defragment process
	$self->defragment();

	# Prompt for what to do next
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

	my ($self) = @_;

	# Initialise parameters
	my $total_loci     = '0';
	my $total_clusters = '0';
	my $settings_ref   = $self->{defragment_settings};
	my $targets_ref    = $settings_ref->{targets};
	$settings_ref->{total_clusters} = $total_clusters;
	$settings_ref->{total_loci}     = $total_loci;

	# Get objects
	my $db        = $self->{db};
	my $reextract = $settings_ref->{reextract};
    my $digs_obj  = $settings_ref->{digs_obj};
	
	# Create a copy of the digs_results table (changes will be applied to copy)
	my $copy_name = $db->backup_digs_results_table();
	my $copy_table_name = $copy_name . '_table';
	$settings_ref->{copy_table_name} = $copy_table_name;

	# Defragment each target file in turn
	print "\n\t # Merging loci in table '$copy_name'\n";
	foreach my $target_ref (@$targets_ref) {
	
		my $organism        = $target_ref->{organism};
		my $target_name     = $target_ref->{target_name};
		my $target_datatype = $target_ref->{target_datatype};
		my $target_version  = $target_ref->{target_version};
		my $t_range         = $settings_ref->{range};
		unless ($t_range) { die;}
		print "\n\n\t # Defragmenting '$target_name' using range '$t_range'";
			
		# Create the relevant set of previously extracted loci
		my @loci;
		my $where  = " WHERE organism      = '$organism' ";
		$where    .= " AND target_datatype = '$target_datatype' ";
		$where    .= " AND target_version  = '$target_version' ";
		$where    .= " AND target_name     = '$target_name' "; 

		$db->get_sorted_digs_results(\@loci, $where);
		my $num_hits = scalar @loci;

		#$devtools->print_hash($settings_ref); die;
		$self->defragment_target_file($target_name, \@loci);
	}
}

#***************************************************************************
# Subroutine:  defragment_target_file
# Description: identify overlapping or contiguous loci in a target file
#***************************************************************************
sub defragment_target_file {

    my ($self, $target_name, $loci_ref) = @_;

	# Set flags
	my $verbose      = $self->{verbose};
	my $settings_ref = $self->{defragment_settings};
    my $reextract    = $settings_ref->{reextract};
  
	# Compose clusters of overlapping/adjacent BLAST hits and extracted loci
	my $total_loci     = $settings_ref->{total_loci};
	my $total_clusters = $settings_ref->{total_clusters};
	my $num_hits = scalar @$loci_ref;
	
	my %target_defragmented;
	$self->compose_clusters(\%target_defragmented, $loci_ref, $settings_ref);
		
	# Get number of clusters for this target file, and increment totals
	my @cluster_ids  = keys %target_defragmented;
	my $num_clusters = scalar @cluster_ids;
	
	$total_loci     = $total_loci + $num_hits;
	$total_clusters = $total_clusters + $num_clusters;

	# Show clusters if verbose flag is set
	if ($verbose) { 
		print "\n\t\t '$target_name': Compressed from '$num_hits' to '$num_clusters' loci";
	}

	# Derive a single, merged locus for clusters of loci that were grouped into a cluster
	my $extended  = $self->merge_clustered_loci(\%target_defragmented, $loci_ref, $reextract);
	print "\n\t\t # $extended extensions to previously extracted sequences ";
	my $num_new   = scalar @$loci_ref;
	print "\n\t\t # $num_new loci to extract after defragment ";
	
	# Extract & reassign
	if ($reextract and $num_new) {
		$self->reextract_and_reassign($settings_ref);
	}
	
	else {  # Update DB
		
		die;
		# Need to rename some fields
		foreach my $hit_ref (@$loci_ref) {
	
			$hit_ref->{extract_start}   = $hit_ref->{start};
			$hit_ref->{extract_end}     = $hit_ref->{end};
			$hit_ref->{sequence}        = 'NULL';
			$hit_ref->{sequence_length} = 0;
			#$devtools->print_hash($hit_ref); die;
		}
		#$digs_obj->update_db(\@loci, $copy_table_name);
	}

	$settings_ref->{total_loci}     = $total_loci;
	$settings_ref->{total_clusters} = $total_clusters;
		
}

#***************************************************************************
# Subroutine:  reextract_and_reassign
# Description: 
#***************************************************************************
sub reextract_and_reassign {

	my ($self) = @_;

	die;

}

############################################################################
# INTERNAL FUNCTIONS: base fxns for merging overlapping/adjacent loci
############################################################################

#***************************************************************************
# Subroutine:  compose_clusters 
# Description: process a sorted list of loci and group into 'clusters' of
#              overlapping feature annotations
#***************************************************************************
sub compose_clusters {

	my ($self, $defragmented_ref, $loci_ref, $settings_ref) = @_;
	
	# Get settings
	my $start_token = $settings_ref->{start};
	my $end_token   = $settings_ref->{end};
	
	# Iterate through loci, grouping them into clusters when they are within range
	my $j = 1;
	my %last_locus;
	my %name_counts;
	my $initialised = undef;
	
	foreach my $locus_ref (@$loci_ref)  {

		# Get locus data
		my $record_id     = $locus_ref->{record_id};
		my $scaffold      = $locus_ref->{scaffold};
		my $target_name   = $locus_ref->{target_name};
		my $assigned_name = $locus_ref->{assigned_name};
		my $probe_gene    = $locus_ref->{probe_gene}; # TODO: document why
		my $assigned_gene = $locus_ref->{assigned_gene};
		my $orientation   = $locus_ref->{orientation};
		my $start         = $locus_ref->{$start_token};
		my $end           = $locus_ref->{$end_token};
		
		# Get last hit values
		my $last_scaffold      = $last_locus{scaffold};
		my $last_start         = $last_locus{$start_token};
	    if ($initialised) {

			# Sanity checking - are sequences in sorted order for this scaffold?
			if ( $scaffold eq $last_scaffold and $start < $last_start) {
				$devtools->print_hash($locus_ref); die; 
			}
           
            # Work out wether to merge this hit with the last
            my $merge = $self->compare_adjacent_loci($locus_ref, \%last_locus, $settings_ref);
            
            unless ($merge) {
                 # Increment the count
                $j++;       
                # Initialise record
                $self->initialise_cluster($defragmented_ref, $locus_ref, $j);
            }
            else {             
                # Extend record
                $self->extend_cluster($defragmented_ref, $locus_ref, $j);
            }
        }
		else {
            $initialised = 'true'; # Set flag - we have seen at least one   
            $self->initialise_cluster($defragmented_ref, $locus_ref, $j);
		}

		# Update last hit data
		$last_locus{record_id}     = $record_id;
		$last_locus{assigned_name} = $assigned_name;
		$last_locus{assigned_gene} = $assigned_gene;
		$last_locus{probe_gene}    = $probe_gene;
		$last_locus{scaffold}      = $scaffold;
		$last_locus{orientation}   = $orientation;
		$last_locus{$start_token}  = $start;
		$last_locus{$end_token}    = $end;
	}	

}

#***************************************************************************
# Subroutine:  compare_adjacent_loci
# Description: compare two loci and determine whether to merge into one
# Note: this is a critical function in many respects - the logic for 
#       merging or not merging loci is implemented here
#***************************************************************************
sub compare_adjacent_loci {

	my ($self, $locus1_ref, $locus2_ref, $settings_ref) = @_;

	# Get settings
	my $range       = $settings_ref->{range};
	my $start_token = $settings_ref->{start};
	my $end_token   = $settings_ref->{end};
	my $verbose     = $self->{verbose};
	my $mode        = $self->{defragment_mode};
	unless ($mode) { die; }
	unless ($range) { die; }
	
	# Get the current hit values
	my $name             = $locus1_ref->{assigned_name};
	my $gene             = $locus1_ref->{assigned_gene};					
	unless ($gene) { # If there is no assigned gene set, use probe gene	
		$gene = $locus1_ref->{probe_gene};
	}			
	my $scaffold         = $locus1_ref->{scaffold};	
	my $start            = $locus1_ref->{$start_token};
	my $end              = $locus1_ref->{$end_token};
	my $orientation      = $locus1_ref->{orientation};			


	# Get the last hit values
	my $last_name        = $locus2_ref->{assigned_name};
	my $last_gene        = $locus2_ref->{assigned_gene};		
	unless ($last_gene) {  # If there is no assigned gene set, use probe gene	
		$last_gene = $locus2_ref->{probe_gene};
	}	
	my $last_scaffold    = $locus2_ref->{scaffold};	
	my $last_start       = $locus2_ref->{$start_token};
	my $last_end         = $locus2_ref->{$end_token};
	my $last_orientation = $locus2_ref->{orientation};			

	
	# Exclude the obvious cases
	if ($scaffold ne $last_scaffold) { return 0; }  # different scaffolds

	# Check orientation
	if ($orientation ne $last_orientation) {
		unless ($mode eq 'consolidate') { 
			if ($verbose) {
				print "\n\t\t Identified pair of loci that are in range, but different orientations";
			}
			return 0;
		}
		else {
			unless ($last_gene and $gene) { die; } # Should never get here		
		}
	}

	# Take action depending on whether we are DEFRAGMENTING or CONSOLIDATING
	if ($gene and $last_gene) { 
		if ($mode eq 'defragment') {
			if ($gene ne $last_gene) { return 0; }  # different genes
		}
		elsif ($mode eq 'consolidate' or $mode eq 'consolidate2') { 
			# do nothing (these loci can be merged, even though different genes)
		}
		else { # Shouldn't get here
			die;
		}
	}
	else { # Shouldn't get here
		print "\n\t\t ERROR genes not found: ene $gene LAST $last_gene";;
		$devtools->print_hash($locus2_ref);
		die; 
	}
	
	# If on same scaffold in same orientation, determine how far apart 
	my $gap = $start - $last_end;		
	if ($verbose) {
		#print "\n\t\t    - Defragment calculation '$scaffold': '$start'-'$last_end' = $gap";
		print "\n\t\t    - Gap between loci = $gap";
	}

	# Test whether to combine this pair of loci into a single merged locus
	if ($gap <= $range) {  # Combine
		if ($verbose) {
			if ($last_name and $last_gene) {
				print "\n\t\t      - Added pair to cluster: $last_name";
				print "($last_gene [$last_orientation]), $name ($gene [$orientation])"; 		
			}
			else { print "\n\t\t      - Added pair to cluster"; }
		}
		return 1;
	}
	else { # Don't combine
		return 0;
	}
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

#***************************************************************************
# Subroutine:  merge_clustered_loci
# Description: resolve each of multiple clustered loci to one locus
#***************************************************************************
sub merge_clustered_loci {
	
	my ($self, $defragmented_ref, $to_extract_ref, $reextract) = @_;

	#$devtools->print_hash($defragmented_ref); die;

	# Get screening database table objects
	my $db_ref           = $self->{db};
	my $searches_table   = $db_ref->{searches_table};
	my $active_set_table = $db_ref->{active_set_table};	
	my $extend_count     = '0';
	my @cluster_ids      = keys %$defragmented_ref;
	foreach my $id (@cluster_ids) {

		# Get data for this cluster
		my $cluster_ref = $defragmented_ref->{$id};
		unless ($cluster_ref) { die; }
		my $num_cluster_loci = scalar @$cluster_ref;
		my $extended = $self->merge_cluster($id, $cluster_ref, $to_extract_ref, $reextract);
		if   ($extended) { $extend_count = $extend_count + $extended; }
	}
	return $extend_count;
}

#***************************************************************************
# Subroutine:  merge_cluster
# Description: resolve a cluster of overlapping loci to one single locus
#***************************************************************************
sub merge_cluster {
	
	my ($self, $cluster_id, $cluster_ref, $to_extract_ref, $reextract) = @_;

	# Determine what to extract for this cluster
	my $verbose = $self->{verbose}; # Get 'verbose' flag setting
	my %new_blast_chains;
	my @merged_results;
	my %previous_digs_result_ids;
	my $highest_end   = undef;
	my $lowest_start  = undef;
	my $previous_digs_result_id = undef;
	my $target_name;		
	my $version;
	my $target_datatype;		
	my $scaffold;
	my $orientation;
	my $organism;
	my $probe_type;
	my $extended;
	my $extract = undef;
	my $num_cluster_loci = scalar @$cluster_ref;		
	if ($verbose) {
		if ($num_cluster_loci > 1) {
			print "\n\t\t    - Merging $num_cluster_loci loci in cluster";
			#$self->show_cluster($cluster_ref, $cluster_id);
		}
	}

	my $num_loci = '0';
	foreach my $locus_ref (@$cluster_ref) {
		
		$num_loci++;			
		my $record_id      = $locus_ref->{record_id};					
		my $digs_result_id = $locus_ref->{digs_result_id};					
		my $start          = $locus_ref->{extract_start};			
		my $end            = $locus_ref->{extract_end};
		$target_name       = $locus_ref->{target_name};					
		$target_datatype   = $locus_ref->{target_datatype};			
		$version           = $locus_ref->{target_version};
		$scaffold          = $locus_ref->{scaffold};			
		$orientation       = $locus_ref->{orientation};
		$organism          = $locus_ref->{organism};
		$probe_type        = $locus_ref->{probe_type};
		unless ($organism) { $organism = $locus_ref->{organism};   }
		unless ($start)    { $start = $locus_ref->{subject_start}; }
		unless ($end)      { $end   = $locus_ref->{subject_end};   }
		if ($verbose) {
			print "\n\t\t    - ID = '$digs_result_id' ($scaffold $orientation $start-$end)";
		}
		#$devtools->print_hash($hit_ref);
							
		# Check if this is a a previously extracted locus
		if ($digs_result_id) {
			$previous_digs_result_ids{$digs_result_id} = $locus_ref;
			$previous_digs_result_id = $digs_result_id;		
		}
		
		# Store this BLAST result as part of a chain if it is new
		else {
			my %data = %$locus_ref;
			$new_blast_chains{$record_id} = \%data; 				
		}
			
		# Record the start and stop parameters so we know whether or not to extend
		if ($lowest_start and $highest_end ) {		
			if ($start < $lowest_start) { $lowest_start = $start; }
			if ($end > $highest_end)    { $highest_end  = $end;   }
		}
		elsif ($start and $end) {
			$highest_end  = $end;
			$lowest_start = $start;
		}
		else { die; } # should never get here			
	}

	# Determine whether or not we need to extract sequences for this cluster
	# Extract if cluster is composed entirely of new loci (no previous result IDs)
	unless ($reextract) {
		unless ($previous_digs_result_id) { 
			if ($verbose) {
				print "\n\t\t # Cluster $cluster_id is comprised entirely of new loci ";
			}
			$extract = 'true';	
		}
	}
	
	# Is this a merge of multiple previously extracted loci?
	my @previous_digs_result_ids = keys %previous_digs_result_ids;
	my $num_previously_extracted_loci_in_cluster = scalar @previous_digs_result_ids;

	if ($num_previously_extracted_loci_in_cluster > 1) {
		my $combined = join (',', @previous_digs_result_ids);
		if ($verbose) {
			print "\n\t\t    - Merging previously extracted loci: ($combined)";
		}
		$extract = 'true';							
	}
	# If it includes a single extracted locus, does this locus need to be extended?	
	elsif ($num_previously_extracted_loci_in_cluster eq 1 and $num_loci > 1) {
			
		# get the indexed query 
		my $ref_id   = shift @previous_digs_result_ids;
		my $data_ref = $previous_digs_result_ids{$ref_id};
		my $start    = $data_ref->{extract_start};			
		my $end      = $data_ref->{extract_end};
		unless ($start)    { $start = $data_ref->{subject_start}; }
		unless ($end)      { $end   = $data_ref->{subject_end};   }
		unless ($lowest_start >= $start and $highest_end <= $end) {	
			$extended++;
			$extract = 'true';
			if ($verbose) {
				print "\n\t\t    - Extending locus: $start, $end: ($lowest_start-$highest_end) ";
			}
		}			
	}

	# If the locus needs to be extracted record the details
	if ($extract) {
		my %extract; # Set the extract params for this cluster
		$extract{target_name}     = $target_name;
		$extract{target_datatype} = $target_datatype;
		$extract{target_version}  = $version;
		$extract{organism}        = $organism;
		$extract{probe_type}      = $probe_type;
		$extract{digs_result_id}  = $previous_digs_result_id;
		$extract{target_name}     = $target_name;
		$extract{scaffold}        = $scaffold;
		$extract{start}           = $lowest_start;	
		$extract{end}             = $highest_end;
		$extract{orientation}     = $orientation;
		$extract{digs_result_ids} = \@previous_digs_result_ids;
		my $num_chains = scalar keys %new_blast_chains;
		if ($num_chains) { $extract{blast_chains} = \%new_blast_chains; }
		push (@$to_extract_ref, \%extract);	
	}
	
	return $extended;
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

	#$devtools->print_hash($defragmented_ref); die;

	my @cluster_ids = keys %$defragmented_ref;
	my $cluster_count;
	foreach my $id (@cluster_ids) {
		$cluster_count++;
		my $cluster_ref = $defragmented_ref->{$id};
		my $cluster_size = scalar @$cluster_ref;
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

 	#$devtools->print_array($cluster_ref); die;	
	#print "\n";
	
	foreach my $locus_ref (@$cluster_ref) {
   		
   		#$devtools->print_hash($hit_ref); die;
		my $organism      = $locus_ref->{organism};
		my $assigned_name = $locus_ref->{probe_name};
		my $assigned_gene = $locus_ref->{probe_gene};
		my $orientation   = $locus_ref->{orientation};
		my $track_name    = $locus_ref->{track_name};
		
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

		print "\n\t\t $cluster_id $organism: ";
		if ($track_name) {
			print "TRACK '$track_name' ";
		}
		print "$assigned_name: $assigned_gene: $scaffold $start-$end ($orientation)";
		if ($digs_result_id) {
			print " (extract ID: $digs_result_id)";
		}
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
	print "\n\t\t\t TOTAL LOCI:     $total_loci";
	print "\n\t\t\t TOTAL CLUSTERS: $total_clusters ";
	
}

#***************************************************************************
# Subroutine:  display_config
# Description: diaplay current defragment configuration
#***************************************************************************
sub display_config {

	my ($self) = @_;

	#my $settings_ref = $self->{defragment_settings};

	my $current_range = $self->{defragment_range};
	unless ($current_range )  { die; } 
	print "\n\n\t\t Current settings (based on control file)";
	print "\n\t\t defragment range: $current_range";
}

############################################################################
# EOF
############################################################################