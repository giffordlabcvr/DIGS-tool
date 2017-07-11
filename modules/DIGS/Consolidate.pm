#!usr/bin/perl -w
############################################################################
# Module:      Consolidate.pm
# Description: Functions for clustering merging matches to different probes
#              into higher order structures
# History:     April  2017: Created by Robert Gifford 
############################################################################
package Consolidate;

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

# Maximum range for Consolidate
my $max = 100000000;
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create new Consolidate 'object'
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Declare empty data structures
	my %crossmatching;

	# Set member variables
	my $self = {

		# Set-up params
		consolidate_range       => $parameter_ref->{consolidate_range}, 
		consolidate_mode        => $parameter_ref->{consolidate_mode}, 
	    consolidate_settings    => $parameter_ref->{consolidate_settings}, 
		
		# Member classes 
		db                     => $parameter_ref->{db},  
		blast_obj              => $parameter_ref->{blast_obj},
		   
		# Paths used in consolidate processes
		genome_use_path        => $parameter_ref->{genome_use_path},
		output_path            => $parameter_ref->{output_path},
		target_groups          => $parameter_ref->{target_groups},
		tmp_path               => $parameter_ref->{tmp_path},

	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# CONSOLIDATION FXNS
############################################################################

#***************************************************************************
# Subroutine:  consolidate_loci
# Description: assemble digs_results rows into higher-order loci 
#***************************************************************************
sub consolidate_loci {

	my ($self) = @_;

	# Get the digs results sorted by scaffold and extract start
	my $db = $self->{db};
	my @sorted;
	$db->get_sorted_digs_results(\@sorted);	
	
	# Compose clusters of overlapping/adjacent BLAST hits and extracted loci
	my $total_loci = scalar @sorted;
	print "\n\t  Consolidating assigned extracted sequences into loci";
	print "\n\t  $total_loci loci in the digs_results table prior to consolidation'";
	my $settings_ref = $self->{consolidate_settings};
	unless ($settings_ref) { die; }
	my %consolidated;
	$self->{defragment_mode} = 'consolidate';
	my $defragment_obj = Defragment->new($self);
	$defragment_obj->compose_clusters(\%consolidated, \@sorted, $settings_ref);
	
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
	my $reextract = 'true';
	#my $reextract = undef;
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
	my $last_assigned_name;

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
		my $record      = "$feature\[$assigned_name($orientation)\]";
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
	unless ($blast_obj) { die; }
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
		print " \n\t Defreag: No target group found for TARGET ID $full_id\n\n"; 
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
# EOF
############################################################################