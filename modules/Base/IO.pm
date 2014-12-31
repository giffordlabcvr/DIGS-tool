#!/usr/bin/perl -w
############################################################################
# Script:       IO.pm 
# Description:  
# History:      Rob Gifford, Novemeber 2006: Creation
############################################################################
package IO;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base classes
use Base::DevTools;

############################################################################
# Globals
############################################################################
	
# For displaying alignment progress
my $increment = 50;

my $devtools = DevTools->new();
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
# Subroutine:   show_output_message 
# Description:  
#***************************************************************************
sub show_output_message {

	my ($self, $message, $output_type) = @_;

	my $f_message;
	unless ($output_type) { 
		$output_type = 'command_line';
	}
	if ($output_type eq 'html') {
		$f_message = "\n\t<BR>$message";
	}
	else {
		$f_message = "\n\t$message";
	}
	print $f_message;
}

#***************************************************************************
# Subroutine:  close_output_section
# Description:  
#***************************************************************************
sub close_output_section {

	my ($self, $output_type) = @_;

	if ($output_type eq 'html') { 
		print "</ul></font>"; 

	}
	else { 
		print "\n\n"; 
	}

}

############################################################################
# EOF
############################################################################
