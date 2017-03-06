**The Database-Integrated Genome Screening (DIGS) Tool, Version 1.0**
------------------------------------------------------------------------------------

### **Database-integrated genome-screening (DIGS)**

Molecular sequence data are highly information rich, and are now being generated much faster than they can be analysed. Consequently, the immense quantities of genome, transcriptome, and metagenome data accumulating in public databases contain multitudes of gene, pseudogene, transposon, virus, and non-coding element sequences that have not been yet been identified, or are only poorly described.

Database-integrated genome screening (DIGS) is an approach for exploring sequence data without relying on previous annotations. DIGS can be used to systematically search for sequences of interest, and to support investigations of their distribution, diversity and evolution.

In DIGS, the output of sequence similarity-based searches is captured in a relational database. This allows for the interrogation and manipulation of output data using structured query language (SQL). In addition, it provides all the benefits of a relational database management system (RDBMS) with respect to features such as data recoverability, multi-user support and network access. 

### **The DIGS tool**

The DIGS tool is a PERL program for implementing DIGS with assembled sequence data (not short read data). It uses the [Basic Local Alignment Search Tool (BLAST)](http://blast.ncbi.nlm.nih.gov/Blast.cgi) to perform sequence similarity searches, and the [MySQL](https://www.mysql.com/) RDBMS to capture their output. 

Instructions for installing and running the DIGS tool can be found on [these pages](https://github.com/robjgiff/DIGS-tool/wiki/Installation-and-Setup)

### QUICK INSTALLATION GUIDE

DIGS requires PERL, BLAST and MySQL.

After downloading this repository, you will next need to:

- Install Perl DBI and DBD::MySQL packages (if they are not already installed)
- Set $DIGS_HOME and $DIGS_GENOMES environment variables
- $DIGS_HOME = path to this repository
- $DIGS_GENOMES = path to the top level of the [target genomes directory](https://github.com/giffordlabcvr/DIGS-tool/wiki/Genome-data)
- Create a MySQL user for DIGS 
- Set up target genomes for screening
- Create probe and reference sequence libraries for screening
- Configure a DIGS control file (a template can be found under 'ctl/template.ctl')

You can then start exploring genomes using the DIGS tool. Execute the pipeline.pl script as follows:

To see options for screening: 

```
./digs_tool.pl -h
```


### Authors and Contributors

**Main developer**: Robert J. Gifford (robert.gifford@glasgow.ac.uk).


**Contributors**: 

Josh Singer (josh.singer@glasgow.ac.uk)

Henan Zhu (h.zhu.1@research.gla.ac.uk)

Tristan Dennis (t.dennis.1@research.gla.ac.uk) 



**DISCLAIMER**

This program may contain bugs, both apparent and less apparent ones. I do not accept responsibility for any problems that arise from use of this software. Use entirely at your own risk.  

