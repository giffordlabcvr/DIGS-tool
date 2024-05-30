# Database-Integrated Genome Screening (DIGS) Tool

## Overview

Systematic in silico genome screening is a powerful approach for identifying genome features within the dark genome. It extends the basic concept of a sequence similarity search in two ways: (i) inclusion of multiple query sequences and/or target databases and (ii) similarity-based classification of matching sequences (“hits”) via comparison to a reference sequence library 

In DIGS, a pipeline for performing systematic in silico genome screening is linked to a relational database, providing a structured and efficient approach.

The DIGS tool is a software framework for implementing DIGS.

## Key Features

- **Efficient In Silico Genome Screening:** The DIGS tool implements a robust, database-integrated approach for systematic genome screening.
- **Relational Database Management System (RDBMS):** Integration with the MySQL RDBMS allows for efficient data exploration and structured querying.
- **Iterative Discovery Process:** The screening process is iterative, with the ability to incorporate novel diversity into the reference sequence library for subsequent screens.

## Implementation and Compatibility

To run the DIGS tool requires PERL, BLAST+, and MySQL (or a supported fork of MySQL such as MariaDB).

The program is accessible through a text-based console interface.

To run the DIGS tool and see options type: `./digs_tool.pl`

## Quick Start

To install DIGS:


To see options for screening: 

```
./digs_tool.pl -h
```

## Installation

Steps involved in installing the DIGS tool and using it to perform DIGS:

1. Install and configure DIGS
    - Download the DIGS tool
    - Install PERL, BLAST, and MySQL
    - Install Perl `DBI` and `DBD::MySQL` packages (if they are not already installed)
    - Set `$DIGS_HOME` and `$DIGS_GENOMES` environment variables
        - `$DIGS_HOME` = path to DIGS tool directory
        - `$DIGS_GENOMES` = path to the top level of the target genomes directory
    - Create a MySQL user for DIGS
    - Set `$DIGS_MYSQL_USER` and `$DIGS_MYSQL_PASSWORD` environment variables

3. Create reference sequence library and set up target sequence databases

4. Create control file for a DIGS project

5. Run the DIGS screen based on the control file

6. Interrogate the output of DIGS 

7. Update reference libraries and repeat steps 4+5 using updated information 

**Step 1** and its sub-components are one-offs associated with initial set-up of the DIGS tool. 

**Steps 2-3** refer to the set-up of individual DIGS projects, and will need to be repeated for each distinct screen.

**Steps 4-6** encapsulate the actual DIGS process. **Step 5** can entail analysis within the screening database (i.e. using SQL), but may also entail the analysis of DIGS output in external programs (e.g. phylogeny packages, statistical analysis programs). Iterating on a DIGS project (**Step 6**) is optional. However, it is anticipated that many DIGS projects will be heuristic in nature, and these will commonly require iteration.

### How long will it take to set up DIGS and run a screen?

#### Step 1: installing the required components for DIGS

In principle this step should be straightforward, but there is some uncertainty as it depends on your platform and operating system. All of the components of DIGS (PERL, BLAST, MySQL) and well-established and widely used within bioinformatics, and in theory should be straightforward to install. A typical bioinformatics server will usually have most if not all of these programs installed already. If installations are required, they should only take a few minutes, especially for experienced bioinformaticians working on LINUX/UNIX operating systems. Away from LINUX, things are less predicable. On Macintosh computers one PERL library (DBD::MySQL) does not come as standard, and may not work 'out-of-the-box'). For more information about installing DIGS on a Mac, please see [here](https://github.com/giffordlabcvr/DIGS-tool/blob/master/md/tristan-mac-installation.md).

#### Step 2: setting up your DIGS screen

Setting up a screen means choosing probes, references and targets, and formatting these for use within DIGS. As much as you may want to get going with screening, it is vital to take some time here. At this point you will need to frame the question that you hope to answer through screening. Which target genomes and why? What kind of sequences are you looking for? Do you expect there to be cross-matching to other kinds of sequences that you're not interested in? 

#### Step 3: creating a control file

This should only take a few minutes. See [here](https://github.com/giffordlabcvr/DIGS-tool/blob/master/md/control-file-structure.md) for information about the control file structure.

#### Step 4: running the screen

Similarity search based screening is a computationally intensive procedure, and can take hours to days to complete. The DIGS output will show how far along the screen is in terms of queries executed, but bear in mind that the length of time each query takes will depend on the size of the target file and the scarcity/abundance of sequences matching the probe in the target file. Where screens are being run on a server, and are expected to take several hours or days to complete, you may want to run them in the background, as described [here](https://github.com/giffordlabcvr/DIGS-tool/blob/master/md/background-screen.md), so that your DIGS run is not dependent on maintaining a connection to the server.

## Overview of the screening process

### Input Data Components

1. **Target Database (TDb):** A collection of whole genome sequence or transcriptome assemblies serving as the target for similarity searches.
2. **Query Sequences (Probes):** Input sequences for similarity searches of the Target Database.
3. **Reference Sequence Library (RSL):** Represents the genetic diversity associated with the genome feature(s) under investigation.

### Organizing genome data for screening

The DIGS tool is designed for screening local sequence databases. This 'target database' (tDB) comprises sets files containing FASTA-formatted nucleotide sequence data that has been indexed for BLAST screening (a function of the BLASTx progra suite). Files in the tDB should be stored in a directory tree with a specific structure. 

```
target_database/species_group/species_name/data_type/version/file_name.fasta
```
The top directory level should be the path to the tDB specified under the DIGS environment variable ($DIGS_GENOMES).

This directory should contain sub-directories that organise tDB files by species group (e.g. plants, animals, mammals, insects). The group names are up to the user.

Within the 'species group' directory, tDB files are grouped according to the species from which they were generated. We advise the use of Latin binomial species names (e.g. homo_sapiens).

The directory level immediately below 'species' is used for organising tDB files by type, e.g. low coverage, wgs_assembly, transcriptome.

Finally, within the 'data_type' level, tDB files should be grouped according to the assembly version.

### Creating a reference sequence library

DIGS requires a library of FASTA-formatted sequences to provide (i) **probes** for the first step pf paired BLAST, and (ii) **reference sequences** for the second step of paired BLAST. 

The DIGS tool uses a simple rule to capture data from the headers of FASTA-formatted reference (and probe) sequences; headers should be structured so as to define two hierarchical name elements; ‘name’ and ‘gene_name’, separated by an underscore. In the example shown above these are a virus name and a gene name. Other two-level hierarchical naming schemes (e.g. species & gene name, gene-subdomain name) can also be used, providing the same scheme is used consistently throughout the project. Reference sequences should be stored in a file, the path to which will be specified in the DIGS control file.

### Selecting probe sequences for BLAST-based screening

A subset of reference sequences should be selected as probes. The entire reference library can be used - in which case there is no need to create a separate file – but it is often sufficient to use only a subset of sequences from the reference sequence library, in which case, a separate file containing this subset should be created. The path to this file is specified in the DIGS control file.

**NOTE:** DIGS can utilize both polypeptide and nucleotide sequences within the same project and database, but separate probe and reference files should be created for nucleotide and polypeptide sequences, and separate screens using distinct control file parameters should be performed for each. 

### Indexing genome data for BLAST

Before running a screen for the first time, you will need to index the target database (TDb) for BLAST searching:

```
./digs_tool.pl –m=1 –i=[path to control file]
```

### Running a sceen

Once the screening database is successfully created, a screen can be executed as follows:

```
./digs_tool.pl –m=2 –i=[path to control file]
```
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
