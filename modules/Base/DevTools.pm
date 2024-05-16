#!/usr/bin/perl -w
############################################################################
# Module:       DevTools.pm 
# Description:  Functions for viewing the contents of PERL data structures
# History:      Rob Gifford, Novemeber 2006: Creation
############################################################################
package DevTools;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::FileIO;

############################################################################
# Globals
############################################################################

# Create base objects
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Parameters
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Member variables
	my $self = {
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# Public Member Functions
############################################################################

#***************************************************************************
# Subroutine:  print_hash
# Description: prints the contents of a hash, including other arrays and
#              hashes.
# Arguments:   $hash_ref:  reference to the hash being printed
#              $formatting: additional formatting for recursive fxn call
#***************************************************************************
sub print_hash {

	my ($self, $hash_ref, $formatting) = @_;
	
	# don't show this line if it's a recursive call
	unless ($formatting) {
		print "\n\t#~o~ ~o~ ~o~ ~o~ ~o~ ~o~# Showing Hash Contents \n";
	}

	while ( my ($key, $value) = each(%$hash_ref) ) {
		
		unless ($value) {
			$value = "\t\t\t\t>>> UNDEFINED <<<";
		}

		if (ref($value) eq 'HASH' ) {
			print "\n\t Hash '$key' in hash:";
			$self->print_hash($value, "\t");
		}	
		elsif (ref($value) eq 'ARRAY') { 
			print "\n\t Array '$key' in hash:";
			$self->print_array($value, "\t");
		}
		else {
			if ($formatting) { print "\n\t$formatting $key => $value"; }
			else             { print "\n\t$key => $value"; }
		}
	}	
	print "\n";
}

#***************************************************************************
# Subroutine:  print_array
# Description: prints the contents of an array, including other arrays and
#              hashes.
# Arguments:   $array_ref:  reference to the hash being printed
#              $formatting: additional formatting for recursive fxn call
#***************************************************************************
sub print_array {

	my ($self, $array_ref, $formatting) = @_;	

	# don't show this line if it's a recursive call
	unless ($formatting) {
		print "\n\t#~o~ ~o~ ~o~ ~o~ ~o~ ~o~# Showing Array Contents \n";
	}

	foreach my $item(@$array_ref) {
		
		if (ref($item) eq 'ARRAY') { 
			print "\n\t Array in array:";
			$self->print_array($item, "\t");
		}
		elsif (ref($item) eq 'HASH' ) {
			print "\n\t Hash in array:";
			$self->print_hash($item, "\t");
		}	
		else {
			chomp $item;
			if ($formatting) { print "\n\t$formatting $item"; }
			else             { print "\n\t$item"; }
		}
	}
	print "\n";
}
#***************************************************************************
# Subroutine:  print_web_hash
# Description: prints the contents of a hash, including other arrays and
#              hashes.
# Arguments:   $hash_ref:  reference to the hash being printed
#              $formatting: additional formatting for recursive fxn call
#***************************************************************************
sub print_web_hash {

	my ($self, $hash_ref, $formatting) = @_;
	
	# don't show this line if it's a recursive call
	unless ($formatting) {
		print "<BR>&nbsp;#~o~ ~o~ ~o~ ~o~ ~o~ ~o~# Showing Hash Contents \n";
	}

	while ( my ($key, $value) = each(%$hash_ref) ) {
		
		unless ($value) {
			$value = "&nbsp;&nbsp;&nbsp;>>> UNDEFINED <<<";
		}

		if (ref($value) eq 'HASH' ) {
			print "<BR>&nbsp; Hash '$key' in hash:";
			$self->print_web_hash($value, "&nbsp;");
		}	
		elsif (ref($value) eq 'ARRAY') { 
			print "<BR>&nbsp; Array '$key' in hash:";
			$self->print_web_array($value, "&nbsp;");
		}
		else {
			if ($formatting) { print "<BR>&nbsp;$formatting $key => $value"; }
			else             { print "<BR>&nbsp;$key => $value"; }
		}
	}	
	print "<BR>";
}

