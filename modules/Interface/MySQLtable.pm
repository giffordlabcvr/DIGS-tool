#!/usr/bin/perl
############################################################################
# Module:      MySQLtable.pm
# Description: 
# History:     (RG) March  2011: Creation 
############################################################################
package MySQLtable;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

############################################################################
# Globals
############################################################################
my $devtools = DevTools->new();
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Create a new table object
#***************************************************************************
sub new {

	my ($invocant, $name, $dbh, $fields_ref) = @_;
	
	my $class = ref($invocant) || $invocant;

	# Member variables
	my $self = {
        name      => $name,
		dbh       => $dbh,
		fields    => $fields_ref,  
	};

	bless ($self, $class);
	
	return $self;
}

############################################################################
# Member Functions
############################################################################

#***************************************************************************
# Subroutine:  insert_row
# Description: generic fxn to submit a row of data to the table
# Arguments:   $sample_data_ref: hash reference containing field values
#***************************************************************************
sub insert_row {

	my ($self, $data_ref) = @_;

	
	my $dbh = $self->{dbh};
	my $field_ref = $self->{fields};
	#$devtools->print_hash($field_ref);
	#$devtools->print_hash($data_ref); die;
	my @fields = keys %$field_ref;
	my $fields = join(',', @fields);	
	my $insert_clause = "INSERT INTO $self->{name} \( $fields\)";
	my @values;
	foreach my $field (@fields) {

		#print "\n\t $field";
		my $value = $data_ref->{$field};
		chomp $value;        # remove newline
		$value =~ s/\s+$//;  # remove trailing whitespace
		my $type  = $field_ref->{$field};
		
		# format the value according to type
		my $type_value = $value;
		if ($type eq 'varchar' or $type eq 'text' or $type eq 'date') {
			$type_value = "'$value'"; # string values should be in quotes
		}
		push (@values, $type_value);
	}

	my $values = join(',', @values);	
	
	my $value_clause = " VALUES \($values\)";
	my $insert = "$insert_clause $value_clause";
	my $sth = $dbh->prepare($insert);
	# DEBUG
	#print "\n\n\t$insert\n\n";
	unless ($sth->execute()) { print $insert; exit; }	
	

	# Get the sample ID (generated via autoincrement)
	my $db_id = $sth->{mysql_insertid};

	return $db_id;
}	

#***************************************************************************
# Subroutine:  select
# Description: Generic select fxn using an array to store rows of data 
# Arguments:   $field_ref: reference to an array with fields to fetch 
#              $data_ref:  reference to an array to store the rows
#              $where:     the where clause of the select as a string
#***************************************************************************
sub select {

	my ($self, $field_ref, $row_ref, $where) = @_;

	my $dbh = $self->{dbh};
	my $fields = join(',', @$field_ref);
	my $query = "SELECT $fields FROM $self->{name} $where";
	my $sth = $dbh->prepare($query);
	unless ($sth->execute()) { print $query; exit;}	
	while (my $row = $sth->fetchrow_arrayref) {
		push(@$row_ref, (join("\t", @$row)));
	}
}

#***************************************************************************
# Subroutine:  select row
# Description: Generic select fxn for a selecting a single row and storing 
#              data in a hash 
# Arguments:   $field_ref: reference to an array with fields to fetch 
#              $data_ref:  reference to a hash to store the data
#              $where:     the where clause of the select as a string
#***************************************************************************
sub select_row {

	my ($self, $field_ref, $data_ref, $where) = @_;

	my $dbh = $self->{dbh};
	my $fields = join(',', @$field_ref);
	my $query = "SELECT $fields FROM $self->{name} $where";
	my $sth = $dbh->prepare($query);
	#print "\n\t QUERY $query\n\n\n";
	$sth->execute();	

	# get the values into a hash
	my $i = 0;
	my @row = $sth->fetchrow_array;
	foreach my $field (@$field_ref) {
		
		$data_ref->{$field} = $row[$i];
		$i++;
	}
}

#***************************************************************************
# Subroutine:  select rows
# Description: Generic select fxn for a one or more rows and storing the
#              data in an array of hashes (one row per hash)
#***************************************************************************
sub select_rows {

	my ($self, $field_ref, $data_ref, $where) = @_;

	my $dbh       = $self->{dbh};
	my $fields = join(',', @$field_ref);	
	
	# remove any quotations
	$fields =~ s/"//g;
	
	my $query = "SELECT $fields FROM $self->{name}";
	if ($where) { $query .= " $where"; }

	my $sth = $dbh->prepare($query);
	unless ($sth->execute()) { print $query; exit; }
	
	# DEBUG
	#print "\n\n\t$query\n\n"; 
	my $row_count = 0;
	while (my $row = $sth->fetchrow_arrayref) {
		
		$row_count++;
		my $i = 0;
		my %row;
		foreach my $field (@$field_ref) {
			
			# Deal with an aliased field
			if ($field =~ m/ AS /) {
              	my @field = split(/\'/,$field);
				my $alias = pop @field; 
                #print "\n\t $field aliased to $alias"; die;
				$field = $alias;
			}
			
			my $value = @$row[$i];
			$row{$field} = $value;
			$i++;
		}
		push (@$data_ref, \%row);
	}
}

#***************************************************************************
# Subroutine:  select distinct
# Description: Generic select distinct 
# Arguments:   $field:     the field to do the 'select distinct' on
#              $data_ref:  reference to an array to store the rows
#              $where:     the where clause of the select as a string
#***************************************************************************
sub select_distinct {
	
	my ($self, $field_ref, $data_ref, $where) = @_;

	unless ($where) { $where = ''; }

	my $dbh = $self->{dbh};
	my $fields = join(',', @$field_ref);	
	
	my $query = "SELECT DISTINCT $fields FROM $self->{name}";
	if ($where) { $query .= " $where"; }
	#print "\n\t #### SELECT DISTINCT QUERY:\n\n $query \n\n";
	my $sth = $dbh->prepare($query);
	$sth->execute();	
	while (my $row = $sth->fetchrow_arrayref) {
		my $i = 0;
		my %row;
		foreach my $field (@$field_ref) {
			my $value = @$row[$i];
			$row{$field} = $value;
			$i++;
		}
		push (@$data_ref, \%row);
	}
}

#***************************************************************************
# Subroutine:  select distinct single field
# Description: Generic select distinct 
# Arguments:   $field:     the field to do the 'select distinct' on
#              $data_ref:  reference to an array to store the rows
#              $where:     the where clause of the select as a string
#***************************************************************************
sub select_distinct_single_field {
	
	my ($self, $field, $data_ref, $where) = @_;

	unless ($where) { $where = ''; }

	my $dbh = $self->{dbh};
	my $query = "SELECT DISTINCT $field FROM $self->{name}";
	if ($where) { $query .= " $where"; }
	#print "\n\t #### SELECT DISTINCT QUERY:\n\n $query \n\n";
	my $sth = $dbh->prepare($query);
	$sth->execute();	
	while (my $row = $sth->fetchrow_arrayref) {
		my $i = 0;
		my %row;
		my $value = @$row[$i];
		push (@$data_ref, $value);
	}
}

#***************************************************************************
# Subroutine:  select value
# Description: use when only a single value returned - e.g. with COUNT(*) 
#***************************************************************************
sub select_value {
	
	my ($self, $field, $data_ref, $where) = @_;

	unless ($where) { $where = ''; }

	my $dbh = $self->{dbh};
	my $query = "SELECT $field FROM $self->{name} $where";
	#print $query;

    my $sth = $dbh->prepare($query);
	$sth->execute();	
	
	while (my $row = $sth->fetchrow_arrayref) {
		my $value = @$row[0];
        push(@$data_ref, $value);
	}
}

#***************************************************************************
# Subroutine:  update
# Description: generic update fxn using hashes 
# Arguments:   $set:   reference to hash with the new values
#              $where: reference to hash with the where relationships
#***************************************************************************
sub update {

	my ($self, $set, $where) = @_;
	
	my $dbh          = $self->{dbh};
	my $fields_ref   = $self->{fields};
	#$devtools->print_hash($fields_ref);die;
	
	my @set_clause;
	my @fields = keys %$set;
	foreach my $field (@fields) {
		my $value = $set->{$field};
		my $type = $fields_ref->{$field};
		unless ($type) { die "no type for field '$field'"; }
		
		my $f_value = $value;
		if ($type eq 'varchar' or $type eq 'text') {
			$f_value = "'$value'";
		}
		
		my $subclause = "$field = $f_value";
		push (@set_clause, $subclause);
	}
	my $set_clause = join (',', @set_clause);

	#my $where_clause = $self->hash_to_delimited_ist($where, ' AND ');
	my $query = "UPDATE $self->{name}
	             SET    $set_clause 
				 $where";
	#print "\n\n\t $query";
	my $sth = $dbh->prepare($query);
	$sth->execute();	
	#exit;
}

#***************************************************************************
# Subroutine:  delete_rows
# Description: generic fxn to delete rows (safer than flush because requires
#              a 'WHERE" statement)
#***************************************************************************
sub delete_rows {

	my ($self, $where) = @_;	
	
	my $dbh = $self->{dbh};
	unless ($where) { die; }
	my $query = "DELETE from $self->{name} $where";
	my $sth = $dbh->prepare($query);
	$sth->execute();	
}

#***************************************************************************
# Subroutine:  flush
# Description: generic fxn to empty table of all data 
#***************************************************************************
sub flush {

	my ($self, $where) = @_;	
	
	my $dbh = $self->{dbh};
	#print "\n\t*** Emptying Table $self->{name} ***\n";
	my $query = "DELETE from $self->{name}";
	if ($where) {
		$query .= $where;
	}
	my $sth = $dbh->prepare($query);
	$sth->execute();	
}

#***************************************************************************
# Subroutine:  reset_primary_key
# Description:  
#***************************************************************************
sub reset_primary_keys {

	my ($self) = @_;
	
	my $dbh = $self->{dbh};
	#print "\n\t*** Resetting Primary Key in $self->{name} ***\n";
	my $alter = "ALTER TABLE $self->{name} AUTO_INCREMENT=1";
	my $sth = $dbh->prepare($alter);
	$sth->execute();	
}

############################################################################
# EOF
############################################################################
