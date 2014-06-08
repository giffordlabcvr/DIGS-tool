############################################################################
# README for Database-Integrated Genome Screening Tool
############################################################################

This package comprises a framework for performing bi-directional BLAST-based
screening of biological sequence databases.

The database-integrated genome screening (DIGS) tool is designed for use with
locally stored genome data. It requires the BLAST+ program suite to perform
sequence similarity searches, and captures data in MySQL database tables.
It also requires PERL with the PERL DBI (database interface) module
installed. 

############################################################################
# Getting Started
############################################################################

############################################################################
#*** STEP ONE: set essential paths and variables:

These instructions assume that PERL is installed locally, and that the
user's intention is to screen a set of FASTA files containing genome contigs
or assembled chromosomes. 

These instructions also assume that the user is reasonably familiar with a
linux/unix programming environment.

The PIPELINE screening framework is controlled using the 'pipeline.pl' script
You must initially set a small number of paths and variables in the header
section of 'pipeline.pl' before it can be used.

If the BLAST+ executables are installed in your local path you do not need
to set a path for BLAST.

If you want to run blast executables from another directory (e.g. a subdirectory
of this directory, set the path to that directory under '$blast_bin_path'.

# NOTE The BLAST+ executables can be obtained from the following URL:
# blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
# ...as of 10th December 2013.

Finally, the path to the top level of a genome data directory structured as 
described below should be set under the environment variable $GENOMES:

(i.e. open ~/.bash_profile and add the following line :

export GENOMES=[full path to your genomes directory here];


############################################################################
#*** STEP TWO: set up sequence libraries

The fundamental components of an analysis using this framework are:

i)   a reference sequence library comprised of aset of nucleotide or amino 
     acid sequences in FASTA format.

ii)  a set of nucleotide or amino acid sequence 'probes', in FASTA format

iii) a set of target sequences to be screened. If the targets are a set of 
     whole genome sequences, these should be contained in a directory
     conforming to the following  structure:

     ---> Genomes [top level directory, can have any name]

        ---> Mammalia [host group level directory, can have any name]

           ---> Homo_sapiens [species directory]

              ---> assembly [datatype directory, e.g. 'wgs', 'assembly']
    
                 ---> ncbi-celera_36.3_march_2008 [version directory]

                      version directory contains the sequence files to be screened.
                      
                      e.g.     hs_alt_Celera_chrX.fa
                               hs_alt_Celera_chrY.fa
                               hs_alt_Celera_chr1.fa
                               
                               ...etc.


############################################################################
#*** STEP THREE: create a control file for your analysis

To start screening genomes, a control file must first be defined for use
with  'pipeline.pl'.

The control file contains information about how to perform screening. It
specifies which sequences to use as probes, which sequence files to screen
in the initial 'search' round of reciprocal BLAST, and which reference 
sequence library to use in the second 'classify' round.

An example control file can be found in this package, under:

./example/ctl/bornavirus_screen.ctl


############################################################################
#*** STEP FOUR: run the screen

To start screening run pipeline.pl as follows

./pipeline.pl -s=[path to your control file]


############################################################################
# CONTENTS OF THE REPOSITORY

./pipeline.pl         # Main process control script

### modules directory -  PERL modules
modules/Base/         # Some base functions (IO etc), used by everything
modules/Interface/    # Interfaces to BLAST and MySQL
modules/DIGS/         # Program-specific functionality

### example directory -  worked examples for DIGS (see DIGS 1.0 user guide)


############################################################################
# LICENSE
############################################################################

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

############################################################################
# END OF FILE
############################################################################
