#!usr/bin/perl -w
############################################################################
# Module:      CrossMatch.pm 
# Description: Capture information about cross-matching during DIGS
# History:     May 2017: Created by Robert Gifford 
############################################################################
package CrossMatch;

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
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: create new CrossMatch 'object'
#***************************************************************************
sub new {

	my ($invocant, $parameter_ref) = @_;
	my $class = ref($invocant) || $invocant;

	# Declare empty data structures
	my %crossmatching;

	# Set member variables
	my $self = {
		
		# Global settings
		process_id             => $parameter_ref->{process_id},
		program_version        => $parameter_ref->{program_version},
		
		# Flags
		verbose                => $parameter_ref->{verbose},
		force                  => $parameter_ref->{force},
		
		# Data structures
		crossmatching          => \%crossmatching,

	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# INTERNAL FUNCTIONS: recording cross-matching
###########################################################################

#***************************************************************************
# Subroutine:  update_cross_matching
# Description: update a hash to record cross-matches
#***************************************************************************
sub update_cross_matching {

	my ($self, $probe_key, $assigned) = @_;
	
	my $crossmatch_ref = $self->{crossmatching};
	
	if ($crossmatch_ref->{$probe_key}) {
		my $cross_matches_ref = $crossmatch_ref->{$probe_key};
		if ($cross_matches_ref->{$assigned}) {
			$cross_matches_ref->{$assigned}++;
		}
		else {
			$cross_matches_ref->{$assigned} = 1;
		}
	}
	else {
		my %crossmatch;
		$crossmatch{$assigned} = 1;
		$crossmatch_ref->{$probe_key} = \%crossmatch;
	}
}

#***************************************************************************
# Subroutine:  show_cross_matching
# Description: show contents of hash that records cross-matches
#***************************************************************************
sub show_cross_matching {

	my ($self) = @_;

	print "\n\n\t  Summary of cross-matching";   
	my $crossmatch_ref = $self->{crossmatching};
	my @probe_names = keys 	%$crossmatch_ref;
	foreach my $probe_name (@probe_names) {
		
		my $cross_matches_ref = $crossmatch_ref->{$probe_name};
		my @cross_matches = keys %$cross_matches_ref;
		foreach my $cross_match (@cross_matches) {
			my $count = $cross_matches_ref->{$cross_match};
			print "\n\t\t #   $count x $probe_name to $cross_match";
		}
	}
}

############################################################################
# EOF
############################################################################
