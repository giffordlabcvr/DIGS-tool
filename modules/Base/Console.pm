#!/usr/bin/perl -w
############################################################################
# Script:       Console.pm 
# Description:  
# History:      Rob Gifford, Novemeber 2006: Creation
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
# Subroutine:  
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


#***************************************************************************
# Subroutine:  show_menu_title
# Description: does what it says 
#***************************************************************************
sub show_menu_title {
	
	my ($self, $title) = @_;
	print "\n\t#=======================# $title \n";
}

############################################################################
# Getting user input via the console 
############################################################################

#***************************************************************************
# Subroutine:  ask_input_file_question
# Description: get the name of an input file and copy the contents to an 
#              array.
# Arguments:   $question: the question to ask 
#              $array_ref: array to copy to
# Returns:     $file: the file name entered by the user
#***************************************************************************
sub ask_input_file_question {

	my ($self, $question, $array_ref) = @_;
	
	print "$question : ";
	my $file = <STDIN>;
	chomp $file;
	unless (open(INFILE, $file)) {
		print "\n\t Cannot open file \"$file\"\n\n";
		return;
	}
	@$array_ref = <INFILE>;
	close INFILE;

	return $file;
}

#***************************************************************************
# Subroutine:  ask_question
# Description: just ask a question and return the input
# Arguments:   $question: the question to ask 
# Returns:     $answer: the response of the user
#***************************************************************************
sub ask_question {

	my ($self, $question) = @_;

	unless ($question) { die; }	
	print "$question : ";
	my $answer = <STDIN>;
	chomp $answer; 
	return $answer;
}


#***************************************************************************
# Subroutine:  ask_question
# Description: ask a question with a deafult option and return the input
# Arguments:   $question: the question to ask 
# Returns:     $answer: the response of the user
#***************************************************************************
sub ask_default_question {

	my ($self, $question, $default) = @_;
	
	print "\n\t Default = '$default'";
	print "$question: (return=default)";
	my $answer = <STDIN>;
	if ($answer =~ /^\s*$/) { 
		$answer = $default;
	} 
	
	chomp $answer; 
	return $answer;
}

#***************************************************************************
# Subroutine:  ask_strict_question
# Description: ask a question and return the input, accepting only answers
#              that conform to a certain format
# Arguments:   $question: the question to ask
#              $answer_format: regular expression string to which the
#                              answer must conform
# Returns:     $answer: the response of the user
# NOTE:        this actually supercedes some of the fxns below, but I have
#              left the others in for now anyway
#***************************************************************************
sub ask_strict_question {

	my ($self, $question, $answer_format) = @_;
	
	my $answer;
	do {
		print "$question : ";
		$answer = <STDIN>;
		chomp $answer; 
	} until ($answer =~ /$answer_format/); 
	return $answer;
}

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

#***************************************************************************
# Subroutine:  ask_simple_choice_question
# Description: ask a question and accept only a range of possible answers 
# Arguments:   $question: the question to ask 
#              $choice_ref: reference to an array with the possible choices
# Returns:     $answer: the integer value entered by the user
#***************************************************************************
sub ask_simple_choice_question {

	my ($self, $question, $choice_ref) = @_;

	# convert the choices array to a scalar
	my $choice_string = join('/', @$choice_ref);
	
	my $answer;
	my $answer_in_set = undef;
	do {
		print "$question " . "\($choice_string\): ";
		$answer = <STDIN>;
		chomp $answer; 
		foreach my $choice (@$choice_ref) {
			#print "\n\t '$choice' '$answer'";
			if ($answer eq $choice) {
				$answer_in_set = 'true';
			}
		}
		
	} until ($answer_in_set);

	return $answer;
}

#***************************************************************************
# Subroutine:  ask_either_or_question
# Description: ask user to choose from one of two possibilities 
# Arguments:   $question: the question to ask 
#              $choice_one, $choice_two: the options (suggest one char each)
# Returns:     $answer: the choice entered by the user
#***************************************************************************
sub ask_either_or_question {

	my ($self, $question, $choice_one, $choice_two) = @_;
	
	my $answer;
	do {
		print "$question ($choice_one/$choice_two): ";
		$answer = <STDIN>;
		chomp $answer; 
	} until ($answer eq $choice_one or $answer eq $choice_two);
	
	return $answer;
}

#***************************************************************************
# Subroutine:  ask_list_question
# Description: ask user to choose an option from a numbered list 
# Arguments:   $question: the question to ask 
#              $list_length: the number of options in the list
# Returns:     $answer: the integer value entered by the user
#***************************************************************************
sub ask_list_question {

	my ($self, $question, $list_length) = @_;

	my $answer;
	my $return = undef;
	do {
		print "$question \(1-$list_length\): ";
		$answer = <STDIN>;
		chomp $answer; 
		if ($answer =~ /^-?\d/) {
			if ($answer <= $list_length) {
				$return = 1;
			}
		} 
	} until ($return eq 1);
	
	return $answer;
}

#***************************************************************************
# Subroutine:  ask_int_question
# Description: ask a question and accept only an integer as a response
# Arguments:   $question: the question to ask 
# Returns:     $answer: the integer value entered by the user
#***************************************************************************
sub ask_int_question {

	my ($self, $question) = @_;
	
	my $answer;
	do {
		print "$question : ";
		$answer = <STDIN>;
		chomp $answer; 
	} until ($answer =~ /\d/); # To do: this isn't strict enough
	return $answer;
}

#***************************************************************************
# Subroutine:  ask_int_with_bounds_question
# Description: ask a question and accept only an integer that falls within
#              a defined range as a response
# Arguments:   $question: the question to ask 
#              $lower_bound, $upper_bound: the specified bounds
# Returns:     $answer: the integer value entered by the user
# TODO:       finish: presently this doesn't discriminate ints and floats
#***************************************************************************
sub ask_int_with_bounds_question {

	my ($self, $question, $lower_bound, $upper_bound) = @_;
	
	my $answer;
	do {
		print "$question \($lower_bound-$upper_bound\): ";
		$answer = <STDIN>;
		chomp $answer; 
	} until ($answer >= $lower_bound and $answer <= $upper_bound);
	return $answer;
}

#***************************************************************************
# Subroutine:  ask_float_question
# Description: ask a question and accept only a floating point value  as a 
#              response.
# Arguments:   $question: the question to ask 
# Returns:     $answer: the integer value entered by the user
# TODO:       finish: presently this doesn't discriminate ints and floats
#***************************************************************************
sub ask_float_question {

	my ($self, $question) = @_;
	
	my $answer;
	do {
		print "$question : ";
		$answer = <STDIN>;
		chomp $answer; 
	} until ($answer =~ /\d/); 
	return $answer;
}

#***************************************************************************
# Subroutine:  ask_float_with_bounds_question
# Description: ask a question and accept only a floating point value
#			   that falls within a defined range as a response
# Arguments:   $question: the question to ask 
#              $lower_bound, $upper_bound: the specified bounds
# Returns:     $answer: the integer value entered by the user
# TODO:       finish: presently this doesn't discriminate ints and floats
#***************************************************************************
sub ask_float_with_bounds_question {

	my ($self, $question, $lower_bound, $upper_bound) = @_;
	
	my $answer;
	do {
		print "$question \($lower_bound-$upper_bound\): ";
		$answer = <STDIN>;
		chomp $answer; 
	} until ($answer >= $lower_bound and $answer <= $upper_bound);
	return $answer;
}

#***************************************************************************
# Subroutine:  show_message
# Description: show a formatted message
# Arguments:   $message: the message to display
#***************************************************************************
sub show_message {

	my ($self, $message) = @_;
	
	my $solid_line  = "\n\t" . '=' x $console_width;
	my $border_line = "\n\t" . '=' . (' ' x ($console_width - 2)) . "=";

	# Format the text
	my $f_message = enclose_box_text($message);

	# Print the box
	print "\n";
	print $solid_line;
	print $f_message;
	print $solid_line;
	print "\n"; 

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
