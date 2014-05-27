
	# LOAD the query files we got
	if ($query_vglue) {
		print "\n\t ### LOADING GLUE QUERY\n";
		$self->load_glue_query($pipeline_obj, $queries_ref);
	}
	elsif ($query_fasta) {
		print "\n\t ### LOADING FASTA QUERY\n";
		$self->load_fasta_query($pipeline_obj, $queries_ref);
		#$devtools->print_hash($queries_ref); exit;
	}
	else { 
		if ($mode) { die "\n\t NO reference sequences found for screen\n\n\n"; }
		else       { die "<br>Error<br>";   }
	}
	
	# Create an indexed DB for assign if we're not using default
	my $ref_fasta_nt = $pipeline_obj->{reference_nt_fasta};
	my $ref_fasta_aa = $pipeline_obj->{reference_aa_fasta};
