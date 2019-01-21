#***************************************************************************
# Subroutine:  merge_cluster
# Description: resolve a cluster of overlapping loci to one single locus
#***************************************************************************
sub merge_cluster {
	
	my ($self, $cluster_id, $cluster_ref, $to_extract_ref) = @_;

	# Determine what to extract for this cluster
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
		#print "\n\t\t    - ID = '$organism' ($scaffold $orientation $start-$end)";
							
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
	#print "\n\t\t # Cluster $cluster_id is comprised entirely of new loci ";
	#$extract = 'true';	
 	
	# Is this a merge of multiple previously extracted loci?
	my @previous_digs_result_ids = keys %previous_digs_result_ids;
	my $num_previously_extracted_loci_in_cluster = scalar @previous_digs_result_ids;

	if ($num_previously_extracted_loci_in_cluster > 1) {
		my $combined = join (',', @previous_digs_result_ids);
		print "\n\t\t    - Merging previously extracted loci: ($combined)";
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
			print "\n\t\t    - Extending locus: $start, $end: ($lowest_start-$highest_end) ";
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

