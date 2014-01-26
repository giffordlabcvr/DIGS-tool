############################################################################
# README for Database-Integrated Genome Screening Tool
############################################################################

This package comprises a framework - implemented in PERL - for performing
so-called 'reciprocal' or 'bi-directional' similarity search-based
screening of biological sequence databases.

In a reciprocal similarity search-based screening, a set of 'probe'
sequences are used to search one or more 'target' sequence files (typically
large files, such as chromosome assemblies, or collections of whole genome
shotgun (WGS) contigs). Sequences in target files that disclose
above-threshold similarity to probes are extracted and 'genotyped' in a
second similarity search, wherein they are compared to a set of 'reference'
sequences representing genetic diversity among the sequences under
investigation.

The database-guided genome screening (DIGS) tool uses the BLAST+ program
suite to perform sequence similarity searches, and captures data in MySQL
database tables.

It requires PERL 5.8 or higher, with the PERL DBI (database interface)
module installed, and is designed for use with locally stored genome data.

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

# Structure of the control file

The control file follows the format of NEXUS files, with blocks defined by
'BEGIN' and 'ENDBLOCK' tokens.

BEGIN SCREENDB;
db_name=Bornavirus_EVEs; # Name of screening DB (avoid special chars)
mysql_server=localhost;  # MySQL server to use (e.g. 'localhost' if using local database) 
mysql_username=;         # Name of user with CREATE, DROP, SELECT, ALTER, DELETE privileges
mysql_password=;         # $mysql_password: Password of the above user
ENDBLOCK;

BEGIN SCREENSETS;
query_aa_fasta=db/probes.fa;       # Path to file with amino acid probe sequences
query_nt_fasta=db/probes.fa;       # Path to file with nucleic acid probe sequences
reference_aa_fasta=db/refseqs.fa;  # Path to file with amino acid reference sequences
reference_nt_fasta=db/refseqs.fa;  # Path to file with nucleic acid reference sequences
bit_score_min_tblastn=100;         # Minimum bit score of tBLASTn hit to extract
bit_score_min_blastn=80;           # Minimum bit score of BLASTn hit to extract
seq_length_minimum=100;            # Minimum sequence length of BLAST hit to extract
ENDBLOCK;

BEGIN PARAMS;           
ENDBLOCK;

BEGIN TARGETS;                     # List of target databases
Mammalia/Homo_sapiens/
Mammalia/Pan_troglodytes/
ENDBLOCK;

BEGIN SCREENSQL;

ENDBLOCK;

Note that not all the parameters shown above need to be defined.

############################################################################
#*** STEP FOUR: run the screen

To start screening run pipeline.pl as follows

./pipeline.pl -s=[path to your control file]


############################################################################
# CONTENTS OF THE REPOSITORY

./pipeline.pl         # Main process control script

### db directory 
screenset/            # defaukt directory to store probe and reference sequences

### process directory
process/              # contains files generated during screening process

### modules directory -  PERL modules
modules/Base/         # Some base functions (IO etc), used by everything
modules/Interface/    # Interfaces to BLAST and MySQL
modules/DIGS/         # Program-specific functionality

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
