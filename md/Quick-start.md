## Quick start

Steps involved in installing the DIGS tool and using it to perform DIGS:

```
  1. Install and configure DIGS 

     - Download the DIGS tool
     - Install PERL, BLAST, and MySQL
     - Set $DIGS_HOME and $DIGS_GENOMES environment variables
     - Install Perl DBI.pm if not installed
     - Create a MySQL user for DIGS

  2. Create reference sequence library and set up target sequence databases

  3. Create control file for a DIGS project

  4. Run the DIGS screen based on the control file

  5. Interrogate the output of DIGS 

  6. Update reference libraries and repeat steps 4+5 using updated information 

```


**Step 1** and its sub-components are one-offs associated with initial set-up of the DIGS tool. 

**Steps 2-3** refer to the set-up of individual DIGS projects, and will need to be repeated for each distinct screen.

**Steps 4-6** encapsulate the actual DIGS process. **Step 5** can entail analysis within the screening database (i.e. using SQL), but may also entail the analysis of DIGS output in external programs (e.g. phylogeny packages, statistical analysis programs). Iterating on a DIGS project (**Step 6**) is optional. However, it is anticipated that many DIGS projects will be heuristic in nature, and these will commonly require iteration.

&nbsp;

**How long will it take to set up DIGS and run a screen?**

*Step 1: installing the required components for DIGS*

In principle this step should be straightforward, but there is some uncertainty as it depends on your platform and operating system. All of the components of DIGS (PERL, BLAST, MySQL) and well-established and widely used within bioinformatics, and in theory should be straightforward to install. A typical bioinformatics server will usually have most if not all of these programs installed already. If installations are required, they should only take a few minutes, especially for experienced bioinformaticians working on LINUX/UNIX operating systems.

Away from LINUX, things are less predicable. On Macintosh computers one PERL library (DBD::MySQL) does not come as standard, and may not work 'out-of-the-box'). We have created a [wiki page dedicated to addressing this issue](https://github.com/giffordlabcvr/DIGS-tool/wiki/Installing-the-DIGS-tool-on-your-mac). Hopefully, these instructions will help resolve such issues, should they arise, within a few minutes.

So far, we have not attempted to install DIGS on a Windows computer, though this is in theory possible, since all the required components including [BLAST](http://www.ncbi.nlm.nih.gov/books/NBK52637/), [MySQL](https://www.mysql.com/why-mysql/windows/) and [PERL](http://learn.perl.org/installing/windows.html) are available for this platform. If you have attempted this, we'd like to hear about your experience and support you if possible - send us an email!

Before setting up your own, bespoke screening pipeline, you may want to test-run DIGS using some of the examples included in this repo.

*Step 2: setting up your DIGS screen*

Setting up a screen means choosing probes, references and targets, and formatting these for use within DIGS. As much as you may want to get going with screening, it is vital to take some time here. At this point you will need to frame the question that you hope to answer through screening. Which target genomes and why? What kind of sequences are you looking for? Do you expect there to be cross-matching to other kinds of sequences that your are not interested in? 

*Step 3: creating a control file*

This should only take a few minutes. You can use the template example ('template.ctl') included in this repo.

*Step 4: running the screen*

Similarity search based screening is a computationally intensive procedure, and can take hours to days to complete. The DIGS output will show how far along the screen is in terms of queries executed, but bear in mind that the length of time each query takes will depend on the size of the target file and the scarcity/abundance of sequences matching the probe in the target file. 

