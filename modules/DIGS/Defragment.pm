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
my $maximum   = 100000000;
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
		
		# Member classes 
		db                     => $parameter_ref->{db},  
		   
		# Paths used in DIGS process
		genome_use_path        => $parameter_ref->{genome_use_path},
		output_path            => $parameter_ref->{output_path},

	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# TOP LEVEL ENTRY 
############################################################################

#***************************************************************************
# Subroutine:  interactive_defragment 
# Description: interactively defragment via the console
#***************************************************************************
sub interactive_defragment {

	my ($self) = @_;

	# Get a list of all the target files from the screening DB
	$self->{defragment_mode} = 'defragment';
	my $db = $self->{db};
	my $digs_results_table = $db->{digs_results_table};
	my @fields = qw [ organism target_datatype target_version target_name ];
	my @targets;
	$digs_results_table->select_distinct(\@fields, \@targets);

	# Settings for clustering
	my %settings;
	$settings{total_loci}     = '0';
	$settings{total_clusters} = '0';
	$settings{range}          = undef;
	$settings{reextract}      = 1;

	# Preview changes 
	my $choice = undef;
	do    { 
		$choice = $self->preview_defragment(\@targets, \%settings);
	}   until ($choice > 1);

    # Apply the changes	if option is chosen
	if    ($choice eq 2) { 
		$self->defragment_target_files(\@targets, \%settings);
	}
	elsif ($choice eq 3) { 
		print "\n"; exit;
	}
	else { die; } # Should never get here
}

#***************************************************************************
# Subroutine:  consolidate_loci
# Description: assemble digs_results rows into higher-order loci 
#***************************************************************************
sub consolidate_loci {

	my ($self, $sorted_loci_ref) = @_;
	
	# Compose clusters of overlapping/adjacent BLAST hits and extracted loci
	my $total_loci = scalar @$sorted_loci_ref;
	print "\n\t  Consolidating assigned extracted sequences into loci";
	print "\n\t  $total_loci loci in the digs_results table prior to consolidation'";
	my $settings_ref = $self->{consolidate_settings};
	unless ($settings_ref) { die; }
	my %consolidated;
	$self->compose_clusters(\%consolidated, $sorted_loci_ref, $settings_ref);
	
	# Check the output
	my @cluster_ids  = keys %consolidated;
	my $num_clusters = scalar @cluster_ids;
	if ($total_loci > $num_clusters) {
		my $range = $settings_ref->{range};
		print "\n\t  $num_clusters clusters of loci within '$range' bp of one another ";
	}
	
	# Update locus data based on consolidated results
	$self->derive_locus_table_from_clustered_digs_results(\%consolidated);
	
	# Return the number of clusters
	return $num_clusters;
}

############################################################################
# MAIN FXNS
############################################################################

#***************************************************************************
# Subroutine:  derive_locus_table_from_clustered_digs_results
# Description: compile locus information  and update the locus tables
#***************************************************************************
sub derive_locus_table_from_clustered_digs_results {

	my ($self, $consolidated_ref) = @_;

	unless ($consolidated_ref) { die; }
	
	# Get parameters and data structures
	my $verbose           = $self->{verbose};
	my $db_ref            = $self->{db};
	my $loci_table        = $db_ref->{loci_table};
	my $loci_chains_table = $db_ref->{loci_chains_table};
	
	# Flags for how to handle
	#my $reextract = undef;
	my $reextract = 'true';
	#my $annotate_ends = undef;
	my $annotate_ends = 'true';

	# Iterate through the clusters	
	my @cluster_ids  = keys %$consolidated_ref;
	foreach my $cluster_id (@cluster_ids) {

		# Get the loci in this cluster
		my $cluster_ref = $consolidated_ref->{$cluster_id};

		# Turn this cluster into an annotated locus
		my %locus;
		$self->derive_locus_structure(\%locus, $cluster_ref);
		#$devtools->print_hash(\%locus); die;

		# Extract the consolidate locus if the flag is set
		if ($reextract) {
			$self->extract_consolidated_locus(\%locus);
		}

		# Do the annotation for truncated versus non-truncated 	matches
		if ($annotate_ends) {
			#$self->annotate_consolidated_locus_flanks(\%locus);
		}

		# Insert the consolidated locus information
		my $locus_array = $locus{locus_array};
		my $locus_structure = join('-', @$locus_array);
		$locus{locus_structure} = $locus_structure;
	
		# Insert the data	
		my $locus_id  = $loci_table->insert_row(\%locus);
				
		# Create the links between the loci and digs_results tables
		foreach my $digs_result_ref (@$cluster_ref) {
			my $digs_result_id = $digs_result_ref->{record_id};
			my %chain_data;
			$chain_data{digs_result_id} = $digs_result_id;
			$chain_data{locus_id}       = $locus_id;
			$loci_chains_table->insert_row(\%chain_data);
		}		
	}
}

#***************************************************************************
# Subroutine:  derive_locus_structure
# Description: derive locus structure based on clustered digs results 
#***************************************************************************
sub derive_locus_structure {

	my ($self, $consolidated_ref, $cluster_ref) = @_;

	my $annotate_flanks = undef;
	my $initialised = undef;
	my $organism;
	my $version;
	my $target_name;
	my $datatype;
	my $assigned_name;
	my $feature;
	my $lowest;
	my $highest;
	my $scaffold;
	my $orientation;
	my $target_datatype;
	my $target_version;
	my @locus_structure;
	my $target_id;
	my $multiple_orientations = undef;
	my $last_element = undef;
	foreach my $element_ref (@$cluster_ref) {
		
		# Capture values from the previous iterations	
		my $last_feature     = $feature;
		my $last_orientation = $orientation;
		#my $last_scaffold    = $scaffold;
		
		# Get the data about this digs_results table row
		my $start       = $element_ref->{extract_start};
		my $end         = $element_ref->{extract_end};

		$feature        = $element_ref->{assigned_gene};
		$assigned_name  = $element_ref->{assigned_name};
		$version        = $element_ref->{target_version};
		$datatype       = $element_ref->{target_datatype};
		$orientation    = $element_ref->{orientation};
		$scaffold       = $element_ref->{scaffold};
		$organism       = $element_ref->{organism};
		$target_name    = $element_ref->{target_name};
		unless ($feature and $orientation) { die; } # Sanity checking
		my $record      = "$feature($orientation)";
		#print "\n\t RECORD $record";

		# Create a target key
		$organism        = $element_ref->{organism};
		$target_name     = $element_ref->{target_name};
		$target_datatype = $element_ref->{target_datatype};
		$target_version  = $element_ref->{target_version};
		my @genome = ( $organism , $target_datatype, $target_version );
		my $this_target_id = join ('|', @genome);
		if ($target_id) {
			unless ($this_target_id eq $target_id) { 
				print "\n\t WHY??? $target_id NE $this_target_id\n\n";
				#die; 
			} 
		}
		$target_id = $this_target_id;

		
		# Deal with first locus in a cluster
		unless ($initialised) {
			$highest      = $end;
			$lowest       = $start;
			$last_element = $element_ref;
			$initialised  = 'true';								
			push(@locus_structure, $record);
			next;		
		}


		# Capture information about coordinates 			
		if ($end > $highest) {
			$highest = $end;
		}
		if ($start < $lowest) {
			$lowest = $start;					
		}

		# Deal with loci that follow at least one previous locus
		if ($orientation eq $last_orientation
		and $feature eq $last_feature) {
			next;
		}
		if ($orientation ne $last_orientation) {
			$multiple_orientations = 'true';
		}

		# Add this element to the start or end of the locus array, based on orientation
		if ($multiple_orientations) { # If multiple orientations, base it on last element
			if ($start >= $last_element->{extract_start}) {
				push(@locus_structure, $record);
			}
			else {
				unshift(@locus_structure, $record);	
			}
		}
		elsif ($orientation eq '+') {
			#my $record = $feature;
			push(@locus_structure, $record);
		}
		elsif ($orientation eq '-') {
			unshift(@locus_structure, $record);			
		}
		$last_element = $element_ref;				
	}

	# Store the data
	$consolidated_ref->{organism}        = $organism;
	$consolidated_ref->{target_version}  = $version;
	$consolidated_ref->{target_name}     = $target_name;
	$consolidated_ref->{target_datatype} = $datatype;
	$consolidated_ref->{scaffold}        = $scaffold;
	$consolidated_ref->{target_id}       = $target_id;
	$consolidated_ref->{orientation}     = $orientation;
	$consolidated_ref->{start}           = $lowest;
	$consolidated_ref->{end}             = $highest;
	$consolidated_ref->{extract_start}   = $lowest;
	$consolidated_ref->{extract_end}     = $highest;
	$consolidated_ref->{assigned_name}   = $assigned_name;
	$consolidated_ref->{assigned_name}   = $assigned_name;
	$consolidated_ref->{locus_array}     = \@locus_structure;
	
}

#***************************************************************************
# Subroutine:  extract_consolidated_locus
# Description: 
#***************************************************************************
sub extract_consolidated_locus {

	my ($self, $consolidated_ref) = @_;

	my $db_ref    = $self->{db};
	my $verbose   = $self->{verbose};
	my $blast_obj = $self->{blast_obj};
	my $seq_len   = 0;
	
	my $genome_use_path  = $self->{genome_use_path};
	my $target_group_ref = $self->{target_groups};
	#my $target_path = $self->get_target_file_path($target_ref);

	my $organism        = $consolidated_ref->{organism};
	my $target_version  = $consolidated_ref->{target_version};
	my $target_datatype = $consolidated_ref->{target_datatype};
	my $target_name     = $consolidated_ref->{target_name};
	my $target_id       = $consolidated_ref->{target_id};
	my $lowest          = $consolidated_ref->{start};
	my $highest         = $consolidated_ref->{end};

	my $full_id = $target_id . '|' . $target_name;
	my $target_group = $target_group_ref->{$full_id};
	unless ($target_group) {
		print " \n\t No target group found for TARGET ID $full_id\n\n"; 
		#$devtools->print_hash($target_group_ref);
        sleep 1;
		return 0;
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

	# Extract the sequence
	#print "\n\t\t    # TARGET: '$target_path'";
	my $sequence   = $blast_obj->extract_sequence($target_path, $consolidated_ref);
	my $seq_length = length $sequence; # Set sequence length
	if ($sequence) {
		
		# If we extracted a sequence, update the data for this locus
		if ($verbose) { print "\n\t\t    - Re-extracted sequence: $seq_length nucleotides "; }
		$consolidated_ref->{sequence}        = $sequence;
		$consolidated_ref->{sequence_length} = $seq_length;
	}
	elsif ($verbose) { 
		print "\n\t\t    # Sequence extraction failed ";
	}
}

#***************************************************************************
# Subroutine:  annotate_consolidated_locus_flanks
# Description: 
#***************************************************************************
sub annotate_consolidated_locus_flanks {

	my ($self, $consolidated_ref) = @_;

	my $db_ref = $self->{db};
	my $contigs_table = $db_ref->{contigs_table};
	my $lowest   = $consolidated_ref->{start};
	my $highest  = $consolidated_ref->{end};
	my $scaffold = $consolidated_ref->{scaffold};

	# Get the length of this contig
	my %data;
	my @fields = qw [ contig_id seq_length ];
	my $where = " WHERE contig_id = '$scaffold'";
	$contigs_table->select_row(\@fields, \%data, $where);
	my $contig_length = $data{seq_length};
	unless ($contig_length) { die; }

	# Check the start of the match
	my @locus_structure;
	if ($lowest eq 1) {  
		unshift(@locus_structure, 'T');
	}
	else {
		unshift(@locus_structure, 'X');
	}
	# Check the end of the match
	if ($highest eq $contig_length) { 
		push(@locus_structure, 'T');
	}
	else {
		push(@locus_structure, 'X');
	}
}

############################################################################
# INTERNAL FUNCTIONS: defragmenting results
############################################################################

#***************************************************************************
# Subroutine:  preview_defragment
# Description: preview a defragmentation process
#***************************************************************************
sub preview_defragment {

	my ($self, $targets_ref, $settings_ref) = @_;

	# Display current settings	
	my $defragment_range = $self->{defragment_range};
	unless ($defragment_range )  { die; } 
	print "\n\n\t\t Current settings (based on control file)";
	print "\n\t\t defragment range: $defragment_range";

	# Get the range	
	my $question1 = "\n\n\t # Set the range for merging hits";
	my $t_range = $console->ask_int_with_bounds_question($question1, $defragment_range, $maximum);		
	#my $t_range = 1000;

	# Preview this defragment
	$self->defragment_digs_results($targets_ref, $settings_ref, $t_range);

	# Summarise results
	my $total_loci     = $settings_ref->{total_loci};
	my $total_clusters = $settings_ref->{total_clusters};
	print "\n\t\t\t TOTAL LOCI:     $total_loci";
	print "\n\t\t\t TOTAL CLUSTERS: $total_clusters ";
	print "\n\n\t\t Option 1: preview new parameters";
	print "\n\t\t Option 2: apply these parameters";
	print "\n\t\t Option 3: exit";

	# Prompt for what to do next
	my $list_question = "\n\n\t # Choose an option:";
	my $choice = $console->ask_list_question($list_question, 3);
	$settings_ref->{range} = $t_range;
	#my $choice = 2;

	return $choice;
}

#***************************************************************************
# Subroutine:  defragment_digs_results
# Description: preview results of a defragment process (for interactive defragment)
#***************************************************************************
sub defragment_digs_results {

    my ($self, $targets_ref, $cluster_params, $t_range) = @_;
   
	# Apply the settings
	my $verbose = $self->{verbose};
	my $total_loci = '0';
	my $total_clusters = '0';
	foreach my $target_ref (@$targets_ref) {

		my $organism        = $target_ref->{organism};
		my $target_name     = $target_ref->{target_name};
		my $target_datatype = $target_ref->{target_datatype};
		my $target_version  = $target_ref->{target_version};
			
		# Create the relevant set of previously extracted loci
		my @loci;
		my $where  = " WHERE organism      = '$organism' ";
		$where    .= " AND target_datatype = '$target_datatype' ";
		$where    .= " AND target_version  = '$target_version' ";
		$where    .= " AND target_name     = '$target_name' "; 

		$self->get_sorted_digs_results(\@loci, $where);
		my $num_hits = scalar @loci;
		
		# Compose clusters of overlapping/adjacent BLAST hits and extracted loci
		my %settings;
		my %target_defragmented;
		$settings{range} = $t_range;
		$settings{start} = 'extract_start';
		$settings{end}   = 'extract_end';
		$self->compose_clusters(\%target_defragmented, \@loci, \%settings);
		
		# Get number of clusters
		my @cluster_ids  = keys %target_defragmented;
		my $num_clusters = scalar @cluster_ids;

		# Show clusters if verbose flag is set
		if ($verbose) { print "\n\n\t\t Interactive defrag: $num_hits hits in target $target_name"; }
		$total_loci     = $total_loci + $num_hits;
		$total_clusters = $total_clusters + $num_clusters;
	}

	$cluster_params->{total_loci}     = $total_loci;
	$cluster_params->{total_clusters} = $total_clusters;
		
}

#***************************************************************************
# Subroutine:  defragment_target_files
# Description: implement a defragmentation process for a set of target files
#***************************************************************************
sub defragment_target_files {

	my ($self, $targets_ref,  $settings_ref) = @_;

	my $genome_use_path  = $self->{genome_use_path};
	my $target_group_ref = $self->{target_groups};
	my $db               = $self->{db};
	my $t_range          = $settings_ref->{range};
	my $reextract        = $settings_ref->{reextract};
	if ($reextract) {
		unless ($genome_use_path and $target_group_ref)  { die; }
	}

	# Create a copy of the digs_results table (changes will be applied to copy)
	# TODO: fix this
	print "\n\t # Defragmenting using range '$t_range'\n";
	my $copy_name = $db->backup_digs_results_table();
	print "\n\t # Copied DIGS results to '$copy_name'\n";
	my $dbh = $db->{dbh};
	$db->load_digs_results_table($dbh, 'digs_results');	
	unless ($db->{digs_results_table}) { die; }
	
	# Iterate through the target files, applying the defragment process to each		
	foreach my $target_ref (@$targets_ref) {

		# Get the target details (and thus the target path)
		#my $target_path = $self->get_target_file_path($target_ref);
		my $organism        = $target_ref->{organism};
		my $target_datatype = $target_ref->{target_datatype};
		my $target_version  = $target_ref->{target_version};
		my $target_name     = $target_ref->{target_name};
		my @genome = ( $organism , $target_datatype, $target_version, $target_name );
		my $target_id       = join ('|', @genome);
		print "\n\t\t # Defragmenting hits in '$target_name'";

		my $target_path = 'NULL';
		if ($reextract) {
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
			$target_path = join ('/', @path);
		}

		# Construct WHERE statement
		my $where  = " WHERE organism      = '$organism' ";
		$where    .= " AND target_datatype = '$target_datatype' ";
		$where    .= " AND target_version  = '$target_version' ";
		$where    .= " AND target_name     = '$target_name' "; 
		$settings_ref->{start}     = 'extract_start';
		$settings_ref->{end}       = 'extract_end';
		$settings_ref->{where_sql} = $where;
		$self->defragment_target($settings_ref, $target_path, 'digs_results');
	}
}

#***************************************************************************
# Subroutine:  defragment_target 
# Description: implement a defragmentation process for a single target file
# Note: similar to compile_nonredundant_locus_set fxn, diff details
#***************************************************************************
sub defragment_target {

	my ($self, $settings_ref, $target_path, $copy_name) = @_;
	
	# Create the relevant set of previously extracted loci
	my %target_defragmented;
	my $loci_ref = $settings_ref->{defragment_loci};
 	my $digs_obj = $settings_ref->{digs_obj};
 	my $num_hits = scalar @$loci_ref;
			
	# Compose clusters of overlapping/adjacent BLAST hits and extracted loci
	$self->compose_clusters(\%target_defragmented, $loci_ref, $settings_ref);
	my @cluster_ids  = keys %target_defragmented;
	my $num_clusters = scalar @$loci_ref;
	if ($num_clusters < $num_hits) {
		#$self->show_clusters(\%target_defragmented);  # Show clusters
		my $range = $settings_ref->{range};
		print "...compressed to $num_clusters overlapping/contiguous clusters within '$range' bp of one another";
	}
	
	# Determine what to extract, and extract it
	my @loci;
	my $reextract = $settings_ref->{reextract};
	my $extended  = $self->merge_clustered_loci(\%target_defragmented, \@loci, $reextract);
	print "\n\t\t # $extended extensions to previously extracted sequences ";
	my $num_new   = scalar @loci;
	print "\n\t\t # $num_new loci to extract after defragment ";
    my $copy_table_name = $copy_name . '_table';
	if ($reextract and $num_new) {

		# Extract newly identified or extended sequences
		my @extracted;
		$self->extract_sequences_from_target_file($target_path, \@loci, \@extracted);
		
		# Do the genotyping step for the newly extracted locus sequences
		my $assigned_count   = 0;
		my $crossmatch_count = 0;
		my $num_extracted = scalar @extracted;
		print "\n\t\t # Genotyping $num_extracted newly extracted sequences:";
		foreach my $hit_ref (@extracted) { # Iterate through loci		
			my $classify_obj = Classify->new($digs_obj);
			$classify_obj->classify_sequence_using_blast($hit_ref);
			$assigned_count++;
			my $remainder = $assigned_count % 100;
			if ($remainder eq 0) { print "\n\t\t\t # $assigned_count sequences classified "; }
		}
		print "\n\t\t\t # $assigned_count sequences classified";

		# Update DB
		my $num_deleted = $self->update_db(\@extracted, $copy_table_name, 1);
		print "\n\t\t\t # $num_deleted rows deleted from digs_results table\n";
	}
	else { # DEBUG
		# Update DB
		$self->prepare_locus_update(\@loci);
		$digs_obj->update_db(\@loci, $copy_table_name);
	}


	return $num_new;
}

############################################################################
# INTERNAL FUNCTIONS: clustering/merging overlapping/adjacent loci
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

#***************************************************************************
# Subroutine:  show_clusters
# Description: print information about clustered loci to the screen
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

############################################################################
# Development
############################################################################

#***************************************************************************
# Subroutine:  prepare_locus_update 
# Description: 
#***************************************************************************
sub prepare_locus_update {

	my ($self, $loci_ref) = @_;

	# Get parameters from self
	foreach my $hit_ref (@$loci_ref) {
	
		$hit_ref->{extract_start}   = $hit_ref->{start};
		$hit_ref->{extract_end}     = $hit_ref->{end};
		$hit_ref->{sequence}        = 'NULL';
		$hit_ref->{sequence_length} = 0;
		#$devtools->print_hash($hit_ref); die;
	}
}

############################################################################
# EOF
############################################################################
