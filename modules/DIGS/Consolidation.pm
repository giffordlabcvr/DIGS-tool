#!/usr/bin/perl -w
############################################################################
# Module:      Consolidation.pm
# Description: Consolidate Extracted hits into Loci
# History:     August 2014: Created by Daniel Blanco-Melo 
############################################################################
package Consolidation;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::FileIO;
use Base::DevTools;
use Base::Console;

# Program components
use DIGS::ScreenBuilder; # Functions to set up screen
#### DBM TESTING
use DIGS::ScreeningDB; # Functions to insert Loci data

############################################################################
# Globals
############################################################################

# Default minimum length of BLAST hit to extract
my $default_min_seqlen = 100; # minimum sequence length of BLAST hit to extract
my $process_dir_warn   = 100; # folders in process directory to trigger warning

# Base objects
my $fileio    = FileIO->new();
my $devtools  = DevTools->new();
my $console   = Console->new();
my $verbose   = undef;
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create new Pipeline 'object'
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {
	
		# Flags
		process_id             => $parameter_ref->{process_id},

		# Paths and member variables
		blast_bin_path         => $parameter_ref->{blast_bin_path},
		genome_use_path        => $parameter_ref->{genome_use_path},
		output_path            => $parameter_ref->{output_path},

		# Database variables
		db_name                => $parameter_ref->{db_name},   
		mysql_server           => $parameter_ref->{server},   
		mysql_username         => $parameter_ref->{username},   
		mysql_password         => $parameter_ref->{password},   

		# Member classes
		blast_obj              => $parameter_ref->{blast_obj},
		db                     => $parameter_ref->{db},

		# Consolidation Variables
		group					=> $parameter_ref->{group},
		length_threshold		=> $parameter_ref->{length_threshold},
		bitscore_min_tblastn	=> $parameter_ref->{bitscore_min_tblastn},
		genome_structure		=> $parameter_ref->{genome_structure}
	};

	bless ($self, $class);
	return $self;

}

#***************************************************************************
# Subroutine: consolidate
# Description: execute the backbone for consolidation and table retrival
#***************************************************************************
sub consolidate {

	my ($self, $scaffold, $orientation, $fam) = @_;
	
	# get parameters for consolidation that were set in the control file
	my $db_name   		 = $self->{db_name};
	my $length_threshold = $self->{length_threshold};
	my $bitscore 		 = $self->{bitscore_min_tblastn};
	my $genome_struc 	 = $self->{genome_structure};
	unless ($db_name)          { die; }
	unless ($length_threshold) { die; }
	unless ($bitscore)         { die; }
	unless ($genome_struc)     { die; }

	my $evalue_exp = 4;
	my %feature_structure;
	$self->make_feature_structure_hash($genome_struc, \%feature_structure);
	#$devtools->print_hash(\%feature_structure); #die;

	my @data_ref;
	$self->load_extract_data($evalue_exp, $bitscore, $db_name, \@data_ref, $scaffold, $orientation, $fam);
	my @organism_array = split(/\//,$scaffold);
	my @scaffold_array = split(/\^/,$organism_array[3]);
	my $head_scaff	=	$scaffold_array[1];

	if(!@data_ref){
		print "No hits on $scaffold $orientation\n";
		return('nohits');
	}

	# if the query executes you'll get an array of hashes as results
	# each hash corresponds to a database row, with the values
	# indexed by the database field names
	my $flagfirsthit=undef;
	my $counter4structure;
	my $flag_contiguous_gene;
	my $previous_hit;       #stores the hash ref of previous hit
	my @previous_assigned_to; #will store the previous assigned to hit (in two pieces Assigned_to = family, Assigned_to_gene = gene)
	my @current_assigned_to; #will store the current assigned to hit
    my @ORFtable;
	my @ORFcons_results;
	my @consolidated_features;
	my $chunk_path = $self->{genome_use_path} . '/' . $self->{group} . '/' . $data_ref[0]->{organism};
	$chunk_path .= '/' . $data_ref[0]->{data_type} . '/' . $data_ref[0]->{version} . '/' . $data_ref[0]->{target_name};

	foreach my $current_hit (@data_ref){ #$current_hit stores the hash ref of the current hit

		#print "$current_hit->{'Extracted.record_id'}\n";
		unless ($flagfirsthit) {
			$previous_hit = $current_hit;
			$flagfirsthit='true';
			next;
		}
		$previous_assigned_to[0]=$previous_hit->{'assigned_name'};
		$previous_assigned_to[1]=$previous_hit->{'assigned_gene'};
		$current_assigned_to[0]=$current_hit->{'assigned_name'};
		$current_assigned_to[1]=$current_hit->{'assigned_gene'};
		#$chunk_path = $self->{genome_use_path} . '/' . $self->{group} . '/' . $previous_hit->{organism};
		#$chunk_path .= '/' . $previous_hit->{data_type} . '/' . $previous_hit->{version} . '/' . $previous_hit->{target_name}; 

		# If both hits are from the same gene
		if($previous_assigned_to[1] eq $current_assigned_to[1]) { #try to consolidate the orf

			#ORFcons returns 1|0 if consolidate and a reference to the hash with new pvl
			@ORFcons_results = $self->ORFcons($orientation, $length_threshold, $previous_hit, $current_hit);
			if($ORFcons_results[0]==0){ #if the two hits were NOT consolidated 

				push(@ORFtable, $previous_hit); #add the previous orf to the ORFtable and try to consolidate a loci
				#consolidate_features retrurns reference to a hash containing the row describing the feature
				push(@consolidated_features,$self->consolidate_features($orientation, $head_scaff, \@ORFtable,$chunk_path));
				@ORFtable=();
				$previous_hit = $ORFcons_results[1];#this is the current_hit
			}
			else{ #if the hits were consolidated
				
				$previous_hit = $ORFcons_results[1]; #this is consolidated previous_hit
				#continue consolidating
			}
		}

		# If the two hits belong to the same family not same gene
		else {
			$flag_contiguous_gene = 0;
			push(@ORFtable, $previous_hit); #if the current hit is from a different gene than the previous
			
			#add the previous orf to the ORFtable
			if(($previous_assigned_to[1] eq 'ltr') && (scalar(@ORFtable) > 1)){
				#SPECIAL CASE: if the previous hit is an LTR 
				#and is not the first one (i.e. is the end of the loci) 
				#prevent to consolidate more members
				$flag_contiguous_gene = 0;
			}
			else { #if either the previous hit is not an LTR or is the first LTR
				for($counter4structure=0;$feature_structure{$previous_assigned_to[1]}[$counter4structure];$counter4structure++){
				#go through the %feature_structure to see if the next gene is the consecutive
					if($current_assigned_to[1] eq $feature_structure{$previous_assigned_to[1]}[$counter4structure]){#if consecuitive gene
						$flag_contiguous_gene = 1;      #turn the flag of consecutive gene on
						
						if(abs(($previous_hit->{'extract_end'})-($current_hit->{'extract_start'}))<= $length_threshold){
							#if the distance between the end of the previous and the begining of the actual hit
							#is smaller than $length_threshold
							$previous_hit = $current_hit;
							#continue consolidating
						}
						else{#if bigger than $length_threshold close and consolidate the loci
							push(@consolidated_features,$self->consolidate_features($orientation, $head_scaff,\@ORFtable,$chunk_path));
							@ORFtable=();
							$previous_hit = $current_hit;
						}
						last;
					}
				}
			}
			if($flag_contiguous_gene == 0){# if the hit was NOT the consecutive gene
				#close and consolidate the loci
				push(@consolidated_features,$self->consolidate_features($orientation, $head_scaff, \@ORFtable,$chunk_path));
				@ORFtable=();
				$previous_hit = $current_hit;
			}
		}
	 }

	#consolidate the last hit
	push(@ORFtable, $previous_hit);
	push(@consolidated_features,$self->consolidate_features($orientation, $head_scaff,  \@ORFtable, $chunk_path));
	return(\@consolidated_features);
	#returns the reference of an array with references to hashes that contain each Feature

}

#***************************************************************************
# Subroutine: make_feature_structure_hash
# Description: creates a hash that set up the order of the elements (ORFs) that consitute the feature
#***************************************************************************
sub make_feature_structure_hash {

	my ($self, $genome_struc, $feature_structure_ref) = @_;

	unless ($genome_struc) { die; }

	my $i;
	my @order = split(/-/, $genome_struc);
	my $first_eq_last = 0;
	my $last = pop(@order);
	my $lclast = lc($last);
	my $first = shift(@order);
	my $lcfirst = lc($first);
	my $lcelement;
	if ($lcfirst eq $lclast){   #if the first element is equal to the last 
                                #(If the feature is a retrovirus the genome structure is falnked by identical LTRs)
		my @array =();
		foreach my $element (@order) {
			$lcelement = lc($element);
			push(@array, $lcelement);
		}
		$feature_structure_ref->{$lcfirst} = \@array;
		$first_eq_last = 1;
	}
	push(@order, $lclast);

	if($first_eq_last == 0){
		unshift(@order,$lcfirst);
	}

	my $number_of_elements = scalar(@order);
	my $one_before_last = ($number_of_elements - 1);
	for($i=0;$i<$one_before_last;$i++){
		my $key = shift(@order);
		my $lckey = lc($key);
		my @array =();
		foreach my $element (@order) {
			$lcelement = lc($element);
			push(@array, $lcelement);
		}	
		$feature_structure_ref->{$lckey} = \@array;
	}
}

#***************************************************************************
# Subroutine:  load_extract_data
# Description: execute the backbone for consolidation and table retrival
#***************************************************************************
sub load_extract_data {
	
	my ($self, $evalue_exp, $bitscore, $db_name, $data_ref, $scaffold, $orientation, $fam) = @_;

	# Connect to DB
	my $data_source = '';
	my $server 	= $self->{mysql_server};
	my $user 	= $self->{mysql_username};
	my $password 	= $self->{mysql_password};
	$data_source = 'dbi:mysql:' . $db_name . ':' . $server;
	my $conexion = DBI->connect($data_source, $user, $password);
	if(!$conexion) {
			print "I cannot connect to data source\n$DBI::errstr\n";
			exit(-2);
	}
	
	my @organism_array = split(/\//,$scaffold);
	my @scaffold_array = split(/\^/,$organism_array[3]);
	# GET ALL THE HITS IN THIS ORIENTATION FROM EXTRACTED
	my @fields = qw (record_id extract_start extract_end subject_start subject_end
			assigned_name assigned_gene organism scaffold target_name data_type version);
	my $fieldssr = join (', ',@fields);
	my $query =	"SELECT $fieldssr 
				FROM Extracted 
				WHERE Orientation = '$orientation'  
				AND Scaffold = '$scaffold_array[1]'
				AND Bit_score > $bitscore 
				AND Organism = '$organism_array[0]'
				AND Data_type = '$organism_array[1]'
				AND Version = '$organism_array[2]'
				AND Target_name = '$scaffold_array[0]'";
	#This Change is to allow for a mixed consolidation
	if($fam){
		$query .= " AND Assigned_name = '$fam'";
	}
	#AND (evalue_exp >= 100 OR evalue_exp = 1)
	#AND Assigned_name = '$fam'
	#AND (evalue_exp >= $evalue_exp OR evalue_exp = 1)" ;
	#AND Bit_score > $bitscore
	# INCORPRATE PARAMETERS FOR THIS SEARCH INTO QUERY
	if ($orientation eq '+') {
			$query .= " ORDER BY extract_start";
	} else{
			$query .= " ORDER BY extract_start DESC";
	}

	# Execute the query
	my $sth = $conexion->prepare($query);
	unless ($sth->execute()) { print $query; exit; }


	my $row_count = 0;
	my ($row,$i,$value, $field);
	if ($orientation eq '+'){
		while ($row = $sth->fetchrow_arrayref) {
			$row_count++;
			$i = 0;
			my %row;
			foreach $field (@fields) {
				$value = @$row[$i];
				if ($field eq 'assigned_gene'){
					#change the name of the genes to be consistent with gene array
					$row{$field} = lc($value);
					$i++;
					next;
				}
				$row{$field} = $value;
				$i++;
			}
			push (@$data_ref, \%row);
		}
	}
	else {
	
		while ($row = $sth->fetchrow_arrayref) {
			
			$row_count++;
			my %row;
			#when -ve the coordinates of start and end at the nucleotide level are swapped
			$row{'record_id'} = @$row[0];
			$row{'extract_end'} = @$row[1];
			$row{'extract_start'} = @$row[2];
			$i = 3;
			foreach $field (@fields) {
				if (($field eq 'record_id')||($field eq 'extract_start')||($field eq 'extract_end')){
						next;
				}
				$value = @$row[$i];
			
				if ($field eq 'assigned_gene'){
					#change the name of the genes to be consistent with gene array
					$row{$field} = lc($value);
					$i++;
					next;

				}
				$row{$field} = $value;
				$i++;
			}
			push (@$data_ref, \%row);
		}
	}
	$conexion->disconnect() || warn "Fail to disconnect from data source\n$DBI::errstr\n";
}

#***************************************************************************
# Subroutine: ORFcons
# Description: execute the consolidation of each ORF and LTR
#If consolidates returns 1 and new consolidated previous_hit
#if not retrives 0 and the current_hit
#returns 1|0 if consolidate and a reference to the hash with new pvl
#***************************************************************************
sub ORFcons {

	my ($self, $orientation, $length_threshold, $previous_hit, $current_hit) = @_;
	my $subject_end         = 'subject_end';
	my $subject_start       = 'subject_start';
	my $realex_end          = 'extract_end';
	my $realex_start        = 'extract_start';
	my ($end_previous_NT, $start_current_NT, $end_current_NT, $comparison,$comparison2, $out);
	#store the hit coordinates at the aminoacid level
	my ($end_previous_AA, $start_current_AA, $end_current_AA);
	$end_previous_AA        = $previous_hit->{$subject_end};
	$start_current_AA       = $current_hit->{$subject_start};
	$end_current_AA         = $current_hit->{$subject_end};
	#store the hit coordinates at the nucleotide level
	$end_previous_NT        = $previous_hit->{$realex_end};
	$start_current_NT       = $current_hit->{$realex_start};
	$end_current_NT         = $current_hit->{$realex_end};

	#SET UP THE COMPARISON DEPENDING OF THE ORIENTATION
	#IN ORDER TO STATE IF ONE HIT GOES AFTER THE OTHER
	#AND NOT OVERLAPPING
	if($orientation eq '-') {
		$comparison = ($start_current_NT < $end_previous_NT);
		$comparison2 = ($end_current_NT < $end_previous_NT);
	}
	else {
		 $comparison = ($start_current_NT > $end_previous_NT);
		$comparison2 = ($end_current_NT > $end_previous_NT);
	}

	if($start_current_AA > $end_previous_AA) { #IF THE CURRENT_START GOES AFTER PREVIOUS_END AT THE AMINOACID LEVEL

		if($comparison) { #IF THE CURRENT_START GOES AFTER PREVIOUS_END AT THE NUCLEOTIDE LEVEL
			
			if(abs($end_previous_NT-$start_current_NT)<= $length_threshold) {       
				#IF THE LENGTH BETWEEN THE TWO HITS IS MINOR THAN $length_threshold
				#CONSOLIDATE ORF (AND STORE NEW CONSOLIDATED HIT IN $previous_hit       
				$previous_hit->{'record_id'} .= '-' . $current_hit->{'record_id'};
				$previous_hit->{'assigned_name'} .= '_' . $current_hit->{'assigned_name'};
				$previous_hit->{$realex_end} = $end_current_NT;
				$previous_hit->{$subject_end} = $end_current_AA;
				$out = 1; #FLAG OF SUCCESFUL CONSOLIDATION
			}
			else { #IF THE LENGTH BETWEEN THE TWO HITS NOT MINOR THAN $length_threshold #NOT CONSOLIDATE
				$previous_hit = $current_hit; #RETURN THE UNCONSOLIDATED CURRENT HIT AS THE PREVIOUS
				$out = 0;
			}

		}
		else { #IF THE CURRENT_START GOES AFTER PREVIOUS_END AT THE AMINOACID LEVEL
			#BUT THE CURRENT_START GOES BEFORE PREVIOUS_END AT THE NUCLEOTIDE LEVEL (THEY ARE OVERLAPPING)
			#CONSOLIDATE ORF (AND STORE NEW CONSOLIDATED HIT IN $previous_hit 
			$previous_hit->{'record_id'} .= '-' . $current_hit->{'record_id'};
			$previous_hit->{'assigned_name'} .= '_' . $current_hit->{'assigned_name'};
			$previous_hit->{$realex_end} = $end_current_NT;
			$previous_hit->{$subject_end} = $end_current_AA;
			$out = 1;
		}
	}
	else {   #IF THE CURRENT_START GOES BEFORE PREVIOUS_END AT THE AMINOACID LEVEL (THEY ARE OVERLAPPING)
		
		if($comparison){ #IF THE CURRENT_START GOES AFTER PREVIOUS_END AT THE NUCLEOTIDE LEVEL
			
			if(abs($end_previous_NT-$start_current_NT)<= $length_threshold){ 
				#IF THE LENGTH BETWEEN THE TWO HITS IS MINOR THAN $length_threshold
				#CONSOLIDATE ORF (AND STORE NEW CONSOLIDATED HIT IN $previous_hit
				$previous_hit->{'record_id'} .= '-' . $current_hit->{'record_id'};
				$previous_hit->{'assigned_name'} .= '_' . $current_hit->{'assigned_name'};
				$previous_hit->{$realex_end} = $end_current_NT;
				$previous_hit->{$subject_end} = $end_current_AA;
				$out = 1;
			}
			else{#IF THE LENGTH BETWEEN THE TWO HITS NOT MINOR THAN $length_threshold
				  #NOT CONSOLIDATE
				 $previous_hit = $current_hit; #RETURN THE UNCONSOLIDATED CURRENT HIT AS THE PREVIOUS
				 $out = 0;
			}
		}
		else{  #IF THE CURRENT_START GOES BEFORE PREVIOUS_END AT THE AMINOACID LEVEL (THEY ARE OVERLAPPING)
				#AND THE CURRENT_START GOES BEFORE PREVIOUS_END AT THE NUCLEOTIDE LEVEL (THEY ARE OVERLAPPING)
			if(($comparison2)||($end_current_NT = $end_previous_NT)){         
				
				#IF THE CURRENT_END IS BIGGER THAN THE PREVIOUS_END OR IF THEY DO FINISH AT THE SME POINT 
				#CONSOLIDATE ORF (AND STORE NEW CONSOLIDATED HIT IN $previous_hit
				$previous_hit->{'record_id'} .= '-' . $current_hit->{'record_id'};
				$previous_hit->{'assigned_name'} .= '_' . $current_hit->{'assigned_name'};
				$previous_hit->{$realex_end} = $end_current_NT;
				$previous_hit->{$subject_end} = $end_current_AA;
				$out = 1;
			}
			else { #IF THE ACTL HIT DO NOT EXPAND THE ORF DISCARD IT
				$out = 1;
			}
		}
	}

	return($out, $previous_hit);#returns 1|0 if consolidate or not and a hash with the new pvl
}

#***************************************************************************
# Subroutine: consolidate_features
# Description: execute the consolidation of Features
##returns reference to a hash containing the row describing the Feature
#push(@consolidated_features,consolidate_features($orientation, \@ORFtable,\%Featurefams));
#Recombination or not will be called by consolidate_features
#***************************************************************************
sub consolidate_features {

	my ($self, $ori, $scaffold, $ORFtable, $chunk_path) = @_;

	my ($reference_to_loci, $flag_1st_member, $lastORF,$flag_family, $counter4names, $ORF, $gene, $current_name, $previous_name);
	my $feature_name;
	my $threshold_chr_end = 200;
	$lastORF = scalar(@{$ORFtable}); #get the number of members in the loci to be consolidated
	$lastORF--; #because the first position of the array is [0]
	my @Feature;
	my @names=();
	my @feature_names;
	my @counter_name_array=();
	my @THEnames=();
	my $scaffold_len;
	my @paths;
	my $gi;
	my $name;

	if($scaffold =~ /\|(\d+)\|/){
			$gi = $1;
	}
	else{ $gi = $scaffold; }
	
	$Feature[0]='';
	$Feature[6]='';
	if($lastORF == 0){ #if the loci is only composed of one member
		
		#format the information of the member into a Feature loci
		$reference_to_loci=$ORFtable->[$lastORF];
		$Feature[0]=$reference_to_loci->{'record_id'};
		if ($ori eq '+'){
				$Feature[1]=$reference_to_loci->{'extract_start'};
				$Feature[2]=$reference_to_loci->{'extract_end'};
		}else{
				$Feature[2]=$reference_to_loci->{'extract_start'};
				$Feature[1]=$reference_to_loci->{'extract_end'};
		}

		$gene=$reference_to_loci->{'assigned_gene'};
		#The following is to get the names off all the Sequences that compose a Locus
		#When The consolidation mode is set on MIXED
		$feature_name=$reference_to_loci->{'assigned_name'};
		if($feature_name=~ /_/){
			@feature_names = ();
			@feature_names = split(/_/,$feature_name);
			#find out what family is represented the most and it will be chosen as the loci family
        	$flag_family = 0;
        	@THEnames=();
			foreach $name (sort(@feature_names)){
        	    #get the first family name
        	    if($flag_family == 0){ 
        	        push(@THEnames, $name);
        	        $previous_name = $name;
        	        $counter4names=1;
        	        $flag_family =1;
        	        next;
        	    }else{
        	        #count the abundance of the families that nake the loci
        	        #and store it into an array representing the abundance
        	        if($previous_name eq $name){
        	            $counter4names++;
        	            $counter_name_array[$counter4names]=$previous_name;
        	        }
        	        else{
        	            $counter_name_array[$counter4names]=$previous_name;
        	            push(@THEnames, $name);
       	           		$previous_name=$name;
                   		$counter4names=1;
       	                $counter_name_array[$counter4names]=$previous_name;
           		    }
	            }
	        }
        	#$Feature[3]=pop(@counter_name_array); #get the last family name (i.e. the most abundant)
        	$current_name=pop(@counter_name_array);
        	$Feature[4]='';
        	#$Feature[5]=$ori;
        	$Feature[4] .= join('/', @THEnames); #concatenate the different family names that compose the loci 
	
		}else{
			$current_name=$feature_name;
			$Feature[4]=$current_name;
		}
		#$current_name=$reference_to_loci->{'assigned_name'};
		$Feature[3]=$current_name;
		#$Feature[4]=$current_name;
		$Feature[5]=$ori;
		$Feature[6]=$gene;
		$Feature[7]=$reference_to_loci->{'target_name'};
		$Feature[8]=$reference_to_loci->{'organism'};
		$Feature[9]=$reference_to_loci->{'data_type'};
		$Feature[10]=$reference_to_loci->{'version'};
	}
	else{ #if the loci contains more than one member
		$flag_1st_member = 0;
		@names=();
		foreach $ORF (@{$ORFtable}){ #go through the members
			if($flag_1st_member==0){ #get the chunk name and the organim only from the first member
				$Feature[7]=$ORF->{'target_name'};
				$Feature[8]=$ORF->{'organism'};
				$Feature[9]=$ORF->{'data_type'};
				$Feature[10]=$ORF->{'version'};
				$flag_1st_member=1;
			}
			$Feature[0] .= $ORF->{'record_id'}; #concatenate a list of the Record_ID of each member
			$gene     = $ORF->{'assigned_gene'};
			$current_name  = $ORF->{'assigned_name'};
			if($current_name=~ /_/){
				@feature_names=();
				@feature_names = split(/_/,$current_name);
				push(@names,@feature_names);
			}else{
				push(@names, $current_name); #put in an array the names of the families that make the loci
			}
			$Feature[6] .= $gene; #concatenate the name of the members
			#push(@names, $current_name); #put in an array the names of the families that make the loci
			if ($ORFtable->[$lastORF] != $ORF){ #separate the record_id and the member names by a '-'
					$Feature[0] .= '-';
					$Feature[6] .= '-';
			}
		}
		#arrange the proper coordinates acording to the orientation
		if ($ori eq '+'){
			$reference_to_loci=$ORFtable->[0];
			$Feature[1]= $reference_to_loci->{'extract_start'}; #get the start from the first member
			$reference_to_loci=$ORFtable->[$lastORF];
			$Feature[2]= $reference_to_loci->{'extract_end'}; #get the end from the last member
		}
		else{
			$reference_to_loci=$ORFtable->[0];
			$Feature[2]= $reference_to_loci->{'extract_start'};
			$reference_to_loci=$ORFtable->[$lastORF];
			$Feature[1]= $reference_to_loci->{'extract_end'};
		}
		#find out what family is represented the most and it will be chosen as the loci family
		$flag_family = 0;
		@THEnames=();
		foreach $name (sort(@names)){
			#get the first family name
			if($flag_family == 0){
				push(@THEnames, $name);
				$previous_name = $name;
				$counter4names=1;
				$flag_family =1;
				next;
			}else{
				#count the abundance of the families that nake the loci
				#and store it into an array representing the abundance
				if($previous_name eq $name){
					$counter4names++;
					$counter_name_array[$counter4names]=$previous_name;
				}
				else{
					$counter_name_array[$counter4names]=$previous_name;
					push(@THEnames, $name);
					$previous_name=$name;
					$counter4names=1;
					$counter_name_array[$counter4names]=$previous_name;
				}
			}
		}
		$Feature[3]=pop(@counter_name_array); #get the last family name (i.e. the most abundant)
		$current_name=$Feature[3];
		$Feature[4]='';
		$Feature[5]=$ori;
		$Feature[4] .= join('/', @THEnames); #concatenate the different family names that compose the loci 
	}

	#format the consolidated loci into the hash results
	my $gen_st = $Feature[6];
	#print "locate $Feature[7]\n";
	#@paths = `locate $Feature[7]`;
	#foreach my $line (@paths){
	#    if (($line =~ /$Feature[8]/) && ($line =~ /fa$/)){
    #        chomp($line);
    #        $chunk_path = $line;
    #    }    
    #}
	#my $paleotools_obj = Paleotools->new();
	#my %chunks_paths;
	#$paleotools_obj->get_path_to_chunk(\%chunks_paths);
	#my $chunk_key = $Feature[8] . '_' . $Feature[7];
	#$chunk_path = $chunks_paths_ref->{$chunk_key};
	#print "./bin/blast/blastdbcmd -db $chunk_path -entry $gi -outfmt %l\n";    
	my $blast_exe = $self->{blast_bin_path} . "blastdbcmd";
	$scaffold_len = `$blast_exe -db $chunk_path -entry $gi -outfmt %l`;
	chomp($scaffold_len);
	unless($scaffold_len){
		print "\n\t$blast_exe -db $chunk_path -entry $gi -outfmt %l\n";
	}
	#print "blastdbcmd $gi = $scaffold_len\n";
	#die;
	if ($ori eq '+'){
		if ($Feature[1] >= $threshold_chr_end){
				$Feature[6] = 'X-' . $gen_st;
		}else{
				$Feature[6] = 'T-' . $gen_st;
		}
	}
	else{
		if ($Feature[2] <= ($scaffold_len - $threshold_chr_end)){
			$Feature[6] = 'X-' . $gen_st;
		}else{
			$Feature[6] = 'T-' . $gen_st;
		}
	}
	if ($ori eq '+'){
		if ($Feature[2] <= ($scaffold_len - $threshold_chr_end)){
			$Feature[6] .= '-X';
		}else{
				$Feature[6] .= '-T';
		}
	}
	else{
		if ($Feature[1] >= $threshold_chr_end){
			$Feature[6] .= '-X';
		}else{
			$Feature[6] .= '-T';
		}
	}

	 my %results = (
		'extracted_id'     => $Feature[0],
		'orientation'      => $Feature[5],
		'extract_start'    => $Feature[1],
		'extract_end'      => $Feature[2],
		'assigned_name'    => $Feature[3],
		'assigned_notes'   => $Feature[4], #the different family names that compose the loci
		'scaffold'         => $scaffold,
		'target_name'      => $Feature[7],
		'organism'         => $Feature[8],
		'genome_structure' => $Feature[6],
		'data_type'        => $Feature[9],
		'version'          => $Feature[10]
	);
	return(\%results); # return a reference to a hash
	
}

#***************************************************************************
# Subroutine:  insert_loci_data 
# Description: 
#***************************************************************************
sub insert_loci_data {
	
	#I add the consolidation file that will store the identities of 
	#the sequences that composed a Locus (when consolidation mode is set to MIXED)
	my ($self, $db_name, $Featuresmain_ref, $consolidation_file) = @_;

	my $db_obj = $self->{db};
	unless ($db_obj) { die; }
	if($consolidation_file){
		open(TAB, ">>$consolidation_file") || die "\n\tCannot open $consolidation_file\n";
	}
	# GET A LIST OF UNIQUE SCAFFOLDS, associated with chunk info
	my $loci_table = $db_obj->{loci_table};
	my $loci_link_table = $db_obj->{loci_link_table};
	unless ($loci_table and $loci_link_table) { die; }
	#print "\n\tI entered here!!!!!!!\n";	
	my ($reftoarray, $reftohash);
	my $insert_clause;
	foreach $reftoarray (@$Featuresmain_ref){ #each entry of Featuresmain is an aarray of hashes for each scaffold

		# Skip if no hits
		if ($reftoarray eq 'nohits') { next; }
		# Iterate through the consolidated hits
		foreach $reftohash (@{$reftoarray}){
			# each entry in the array is a hash with every loci
			#$devtools->print_hash($reftohash); die;
			my $locus_id = $loci_table->insert_row($reftohash);
			my $id_string = $reftohash->{'extracted_id'};
			my @extract_ids = split('-', $id_string);
			foreach my $extract_id (@extract_ids) {
				my %link_data;
				print "\n\t ID $id_string";
				print "\n\t locus ID $locus_id extract $extract_id";
				$link_data{locus_id}      = $locus_id;
				$link_data{extracted_id}  = $extract_id;
				#$devtools->print_hash(\%link_data); die;
				$loci_link_table->insert_row(\%link_data);
			}
			my $names_string = $reftohash->{'assigned_notes'};
            my @assigned_names = split('/', $names_string);
			my $assigned_name = $reftohash->{'assigned_name'};
			if($consolidation_file){
				print TAB "$locus_id\t$assigned_name";
            	foreach my $name (@assigned_names) {
                	print TAB "\t$name";
				
            	}
				print TAB "\n";
			}   
		}
	}

	if($consolidation_file){close(TAB);}
}

