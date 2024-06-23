# Database-Integrated Genome Screening (DIGS) Tool

<img src="md/logo_digs.png" align="right" alt="" width="220" />


Welcome to the GitHub repository for the **DIGS Tool**!

**Systematic, sequence similarity search-based genome screening** is a powerful approach for identifying and characterising genome features in silico. This approach extends the basic [sequence similarity search](https://blast.ncbi.nlm.nih.gov/) by: 

 1. Performing multiple searches systematically, involving various **query sequences** and/or **target databases**.
 2. Classifying “**hits**” (matching sequences) via comparison to a **reference sequence library** curated by the investigator.

**Database-integrated genome screening (DIGS)** is a form of systematic genome screening in which a sequence similarity search-based screening pipeline is linked to a **relational database management system** ([RDBMS](https://www.w3schools.com/mysql/mysql_rdbms.asp)). This provides a robust foundation for implementing large-scale, automated screens, and enables a 'database querying' approach to investigating screening output.

**The DIGS Tool is a software framework for implementing DIGS on UNIX/LINUX and Macintosh OSX platforms**. The program is accessible through a text-based console interface. It uses the [BLAST+ ](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) program suite to perform similarity search-based screening, and the [MySQL](https://dev.mysql.com/downloads/mysql/) RDBMS to capture screen output. 

## Overview 

To run the DIGS tool requires [PERL](https://www.perl.org/), [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and [MySQL](https://dev.mysql.com/downloads/mysql/) (or a supported fork of MySQL such as MariaDB).

Steps involved in installing the DIGS tool and using it to perform DIGS are as follows:

1. Install and configure DIGS
    - [Download](https://github.com/giffordlabcvr/DIGS-tool/zipball/master) the DIGS tool
    - Install [PERL](https://www.perl.org/), [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and [MySQL](https://dev.mysql.com/downloads/mysql/)
    - Install Perl `DBI` and `DBD::MySQL` packages (if they are not already installed)
    - Set `$DIGS_HOME` and `$DIGS_GENOMES` environment variables
        - `$DIGS_HOME` = path to DIGS tool directory
        - `$DIGS_GENOMES` = path to the top level of the target database (tDB) directory (see below)
    - Create a MySQL user for DIGS
    - Set `$DIGS_MYSQL_USER` and `$DIGS_MYSQL_PASSWORD` environment variables

3. Create reference sequence library (RSL) and set up target database

4. Create a [control file](https://github.com/giffordlabcvr/DIGS-tool/wiki/DIGS-Tool-Control-File) for a DIGS project

5. Run the DIGS screen based on the control file

6. Interrogate the output of DIGS 

7. Update reference libraries and repeat steps 4+5 using updated information 

**Step 1** and its sub-components are one-offs associated with initial set-up of the DIGS tool. 

**Steps 2-3** refer to the set-up of individual DIGS projects, and will need to be repeated for each distinct screen.

**Steps 4-6** encapsulate the actual DIGS process. **Step 5** can entail analysis within the screening database (i.e. using [SQL](https://github.com/giffordlabcvr/DIGS-tool/wiki/Example-SQL), but may also entail the analysis of DIGS output in external programs (e.g. phylogeny packages, statistical analysis programs). Iterating on a DIGS project (**Step 6**) is optional. However, it is anticipated that many DIGS projects will be heuristic in nature, and these will commonly require iteration.

Please see the [User Guide](https://github.com/giffordlabcvr/DIGS-tool/wiki) for more details.

## Quick Start

To see options for screening: 

`./digs_tool.pl -h` 

### Input Data Components

1. **Target Database (TDb):** A collection of whole genome sequence or transcriptome assemblies serving as the target for similarity searches.
2. **Query Sequences (Probes):** Input sequences for similarity searches of the Target Database.
3. **Reference Sequence Library (RSL):** Represents the genetic diversity associated with the genome feature(s) under investigation.

Before running a screen for the first time, you will need to index the target database (TDb) for BLAST searching:

`./digs_tool.pl –m=1 –i=[path to control file]`

### Running a sceen

Once the target database has been indexed, a screen can be executed as follows:

`./digs_tool.pl –m=2 –i=[path to control file]`

Progress is written to the terminal, and can also be monitored by issuing SQL queries against the relevant screening database. A screen can be stopped at any time. The next time the tool is restarted, it will initiate screening at the point it left off.


## Contributing

The DIGS tool team is very open to further development of this software by the open source bioinformatics community. It is probably worth raising any ideas you have with the team before embarking on development. 

If contributing to the DIGS tool, please review our [Contribution Guidelines](./md/CONTRIBUTING.md).

[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](./md/code_of_conduct.md) 

## Contact

For questions, issues, or feedback, please contact us at [digstool@gmail.com](mailto:digstool@gmail.com) or open an [issue](https://github.com/giffordlabcvr/DIGS-tool/issues).

## Credits

The DIGS tool was written by Robert J. Gifford.

## License

The project is licensed under the [GNU Affero General Public License v. 3.0](https://www.gnu.org/licenses/agpl-3.0.en.html)
