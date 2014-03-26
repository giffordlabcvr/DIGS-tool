#!/usr/bin/perl -w
############################################################################
# Module:       Console.pm 
# Description:  Functions for text console programs 
# History:      Rob Gifford, November 2006: Creation
############################################################################
package Console;

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
my $fileio   = FileIO->new();
my $console_width = 70; # assumed width of the console
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create new Console.pm object 
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
# Subroutine:  refresh_console
# Description: Clear the screen by sending the 'clear' command via 'system'
#***************************************************************************
sub refresh {

	my ($self) = @_;
	my $command = 'clear';
	system $command;
}

#***************************************************************************
# Subroutine:  show_about_box 
# Description: show a formatted title box for a console application
# Arguments:   the program description as a series of strings:
#              - $title, $version, $description, $author, $contact
#***************************************************************************
sub show_about_box {

	my ($self, $title, $version, $description, $author, $contact) = @_;

	my $solid_line  = "\n\t" . '#' x $console_width;
	my $border_line = "\n\t" . '#' . (' ' x ($console_width - 2)) . "#";

	# Format the text
	my $title_version   = $title . ' ' . $version; 
	
	my $f_title_version = enclose_box_text($title_version);
	my $f_description   = enclose_box_text($description);
	my $f_author        = enclose_box_text($author);
	my $f_contact       = enclose_box_text($contact);

	# Print the box
	print "\n\n";
	print $solid_line;
	print $border_line;
	print $f_title_version;
	print $f_description;
	print $f_author;
	print $f_contact;
	print $border_line;
	print $solid_line;
	print "\n\n"; 
}

############################################################################
# Getting user input via the console 
############################################################################

#***************************************************************************
# Subroutine:  ask_yes_no_question
# Description: ask a question and accept only 'y' or 'n' as an answer 
# Arguments:   $question: the question to ask 
# Returns:     $answer: the integer value entered by the user
#***************************************************************************
sub ask_yes_no_question {

	my ($self, $question) = @_;
	
	my $answer;
	do {
		print "$question " . "\(y\/n\): ";
		$answer = <STDIN>;
		chomp $answer; 
	} until ($answer eq 'n' or $answer eq 'y');

	return $answer;
}

############################################################################
# Private Member Functions
############################################################################

#***************************************************************************
# Subroutine:  enclose_box_text
# Description: Format text for an about box by centering it within a box 
#***************************************************************************
sub enclose_box_text {

	my ($text) = @_;

	my $f_text;
	my $left_spacing;
	my $right_spacing;
	my $text_length = length $text;
	
	if ($text_length > ($console_width - 4)) {
		die ("\n\t Title field was more than max length");
	
	}
	else {
		# calculate total white space
		my $space = ($console_width - ($text_length + 2));
		
		# use this value to centre text
		$left_spacing = $space / 2;
		my $adjust_for_uneven = $space % 2;
		$right_spacing = ($space / 2) + $adjust_for_uneven;
	}

	$f_text  = "\n\t#" . (' ' x $left_spacing);
	$f_text .= $text;
	$f_text .= (' ' x $right_spacing) . "#";

	return $f_text;
}

############################################################################
# EOF
############################################################################
