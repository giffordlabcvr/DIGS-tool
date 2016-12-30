**The Database-Integrated Genome Screening (DIGS) Tool, Version 1.0**
------------------------------------------------------------------------------------

### **Database-integrated genome-screening (DIGS)**

Much of the content of genomes consists of poorly characterised '**dark matter**', such as [transposons](https://www.broadinstitute.org/education/glossary/transposable-elements), [pseudogenes](http://pseudogene.org/background.php), [endogenous viral elements (EVEs)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1001191) and obscure [non-coding DNA elements](http://www.pbs.org/newshour/rundown/junk-dna/). These sequences, even when non-functional, contain a wealth of **useful biological information** that can be explored by using [sequence similarity searches](http://www.ebi.ac.uk/Tools/sss/) in combination with strategically chosen reference sequence datasets.

In **database-integrated genome-screening (DIGS)**, the output of sequence similarity search-based genome 'screens' is captured in a [relational database](https://docs.oracle.com/javase/tutorial/jdbc/overview/database.html). This facilitates the implementation of automated screens that can be performed on a large scale. In addition, it allows for the interrogation and manipulation of output data using [structured query language (SQL)](http://www.w3schools.com/sql/sql_intro.asp), and provides all the benefits of a relational database management system (RDBMS) with respect to features such as data recoverability, multi-user support and network access.

### **The DIGS tool**

The DIGS tool provides a computational framework for implementing DIGS. The tool is written in [PERL](https://www.perl.org/). It uses the [Basic Local Alignment Search Tool (BLAST)](http://blast.ncbi.nlm.nih.gov/Blast.cgi) to perform sequence similarity searches, and the [MySQL](https://www.mysql.com/) RDBMS to capture results and track progress.

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

To initially create your screening database:
```
./digs_tool.pl -m=1 -i=[path to your control file]
```

To initiate DIGS:

```
./digs_tool.pl -m=2 -i=[path to your control file]
```




### Authors and Contributors

**Main developer**: Robert J. Gifford (robert.gifford@glasgow.ac.uk).


**Contributors**: 

Josh Singer (josh.singer@glasgow.ac.uk)

Henan Zhu (h.zhu.1@research.gla.ac.uk)

Tristan Dennis (t.dennis.1@research.gla.ac.uk) 



**DISCLAIMER**

This program may contain bugs, both apparent and less apparent ones. I do not accept responsibility for any problems that arise from use of this software. Use entirely at your own risk.  

