#!/usr/bin/perl -w
############################################################################
# Module:      Tasks.pm
# Description: Scheduled tasks  
# History:     November 2011: Created by Robert Gifford 
############################################################################
package Tasks;

############################################################################
# Import statements/packages (externally developed packages)
############################################################################
use strict;

############################################################################
# Import statements/packages (internally developed packages)
############################################################################

# Base Modules
use Base::Console;
use Base::DevTools;
use Base::FileIO;
use Base::BioIO;
use Base::SeqIO;

# Component Modules
use Component::GLUE::Sequence;
use Component::GLUE::RefSeq;
use Component::GLUE::RefSeqParser;
use Component::GLUE::RefSeqLibrary;

############################################################################
# Globals
############################################################################

# Base objects
my $fileio    = FileIO->new();
my $devtools  = DevTools->new();
my $seqio     = SeqIO->new();
my $seq_obj   = Sequence->new();
my $console   = Console->new();
my $writer    = HTML_Utilities->new();
1;

############################################################################
# LIFECYCLE
############################################################################

#***************************************************************************
# Subroutine:  new
# Description: Parameters
#***************************************************************************
sub new {

	my ($invocant, $parameters) = @_;
	my $class = ref($invocant) || $invocant;

	# Set member variables
	my $self = {
		
		# Member variables
		output_path           => $parameters->{output_path},
		microlineage_path     => $parameters->{microlineage_path},
		refseq_lib_path       => $parameters->{refseq_lib_path},
		refseq_use_path       => $parameters->{refseq_use_path},
		blast_db_path         => $parameters->{blast_db_path},
		blast_orf_lib_path    => $parameters->{blast_orf_lib_path},
		blast_nt_orf_lib_path => $parameters->{blast_nt_orf_lib_path},
		blast_utr_lib_path    => $parameters->{blast_utr_lib_path},
		blast_genome_lib_path => $parameters->{blast_genome_lib_path},
		
	
	};
	
	bless ($self, $class);
	return $self;
}

############################################################################
# VIEW OPTIONS 
############################################################################

#***************************************************************************
# Subroutine:  run_task
# Description: top level handler for scheduled tasks
#***************************************************************************
sub run_task {

	my ($self, $task, $logpath) = @_;
	
   	my $USAGE .= "\n\t  v = refresh:\n";
   	$USAGE    .= "\n\t     1 = daily"; 
   	$USAGE    .= "\n\t     2 = weekly";
   	$USAGE    .= "\n\t     3 = monthly"; 
   	$USAGE    .= "\n\n";
	
	if ($task eq 1) { 
		$self->run_daily_tasks();
	}
	elsif ($task eq 2) { 
		$self->run_weekly_tasks();
	}
	elsif ($task eq 3) { 
		$self->run_monthly_tasks();
	}
	else { die $USAGE; }	
}

#***************************************************************************
# Subroutine:  run_daily_tasks 
# Description: 
#***************************************************************************
sub run_daily_tasks {

	my ($self, $task, $logpath) = @_;
	

}

#***************************************************************************
# Subroutine:  run_weekly_tasks 
# Description: 
#***************************************************************************
sub run_weekly_tasks {

	my ($self, $task, $logpath) = @_;
	
	#print "\n\t # Cleanup site\n\n";
	#$self->clean_up();
	#$self->test_phylogeny($run_params_ref);    
	#$self->test_alignment_analysis($run_params_ref);    
	#$self->test_reanalyze($run_params_ref);    
	#$self->test_thirdparty_loading($run_params_ref);    
	#$self->test_preprocessing($run_params_ref);    
	#$self->run_pipeline_tests(); 


}

#***************************************************************************
# Subroutine:  run_monthly_tasks 
# Description: 
#***************************************************************************
sub run_monthly_tasks {

	my ($self, $task, $logpath) = @_;
	

}

#***************************************************************************
# Subroutine:  run_pipeline_tests	
# Description: 
#***************************************************************************
sub run_pipeline_tests {

	my ($self) = @_;

	# RUN an VGLUE screen
	#my $test_file1 = 'ctl_files/screening/TEST_1_HML2_HsapCh15.ctl';
	my $test_file1 = 'ctl_files/screening/REV.ctl';
	my $command = "$0 -s=$test_file1";
	print "\n\t Executing Batchfile: $command";
	system $command;
}

#***************************************************************************
# Subroutine:  clean_up 
# Description:
#***************************************************************************
sub clean_up {

	my ($self) = @_;
	
	# Clean up the report directory
	die;
	my $reports_command = "rm -rf ./site/reports/*";
	system $reports_command;

	# Clean up BLAST DBs
	my $blast_command = "rm -rf ./db/blast/*";
	#system $blast_command;
	
	# Clean up tmp path
	my $tmp_path = $self->{tmp_path};
	my $tmp_command = "rm -rf ./db/blast/*";
	system $tmp_command;

	# Clean up log files
	my $formatlog_command = "rm -rf ./formatdb.log*";
	system $formatlog_command;
	#my $errorlog_command = "rm -rf ./errordb.log*";
	#system $errorlog_command;
}

############################################################################
# EOF 
############################################################################
