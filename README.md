**The Database-Integrated Genome Screening (DIGS) Tool, Version 1.0**
------------------------------------------------------------------------------------

### **Database-integrated genome-screening (DIGS)**

Molecular sequence data are highly information rich, and are now being generated much faster than they can be analysed. Consequently, the immense quantities of genome, transcriptome, and metagenome data accumulating in public databases contain multitudes of gene, pseudogene, transposon, virus, and non-coding element sequences that have not been yet been identified, or are only poorly described.

Database-integrated genome screening (DIGS) is an approach for exploring sequence data.
In DIGS, the output of sequence similarity-based searches is captured in a relational database.
This allows for the interrogation and manipulation of output data using structured query language (SQL).
In addition, it provides all the benefits of a relational database management system (RDBMS)
with respect to features such as data recoverability, multi-user support and network access. 

### **The DIGS tool**

The DIGS tool is a PERL program for implementing DIGS.
It is designed for use with assembled sequence data (i.e. not short read data).
It uses the [Basic Local Alignment Search Tool (BLAST)](http://blast.ncbi.nlm.nih.gov/Blast.cgi) to perform sequence similarity searches,
and the [MySQL](https://www.mysql.com/) to capture their output. 

The DIGS tool can be used to systematically search for sequences of interest, and to support investigations of their distribution, diversity and evolution.

### QUICK INSTALLATION GUIDE

DIGS requires PERL (with DBD::MySQL package), BLAST and MySQL.

To install DIGS:

- Install Perl DBI and DBD::MySQL packages (if they are not already installed)
- Set $DIGS_HOME and $DIGS_GENOMES environment variables
- $DIGS_HOME = path to this repository
- $DIGS_GENOMES = path to the top level of the [target genomes directory](https://github.com/giffordlabcvr/DIGS-tool/wiki/Genome-data)
- Create a MySQL user for DIGS
- Set $DIGS_MYSQL_USER and $DIGS_MYSQL_PASSWORD environment variables

To see options for screening: 

```
./digs_tool.pl -h
```

Detailed instructions for installing and running the DIGS tool can be found on [these pages](https://github.com/robjgiff/DIGS-tool/wiki/Installation-and-Setup)

### Authors and Contributors

**Lead Developer**:

Robert J. Gifford (robert.gifford@glasgow.ac.uk).


**Contributors**: 

Daniel Blanco-Melo (daniel.blancomelo@mssm.edu) 

Henan Zhu (h.zhu.1@research.gla.ac.uk)


**DISCLAIMER**

We do not accept responsibility for any problems that arise from use of this software.

