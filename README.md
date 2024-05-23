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

The DIGS tool is implemented using the [PERL](https://www.perl.org/) scripting language and is compatible with UNIX/LINUX operating systems.

The DIGS tool uses the [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) program suite to perform similarity searches, and the [MySQL](https://dev.mysql.com/downloads/mysql/) relational database management system to capture their output.

To run the DIGS tool requires PERL, BLAST+, and MySQL (or a supported fork of MySQL such as MariaDB).

The program is accessible through a text-based console interface.

To run the DIGS tool and see options type: `./digs_tool.pl`


## Quick Start

To install DIGS:

- Install Perl `DBI` and `DBD::MySQL` packages (if they are not already installed)
- Set `$DIGS_HOME` and `$DIGS_GENOMES` environment variables
- `$DIGS_HOME` = path to this repository
- `$DIGS_GENOMES` = path to the top level of the target genomes directory
- Create a MySQL user for DIGS
- Set `$DIGS_MYSQL_USER` and `$DIGS_MYSQL_PASSWORD` environment variables

To see options for screening: 

```
./digs_tool.pl -h
```

## Input Data Components

1. **Target Database (TDb):** A collection of whole genome sequence or transcriptome assemblies serving as the target for similarity searches.
2. **Query Sequences (Probes):** Input sequences for similarity searches of the Target Database.
3. **Reference Sequence Library (RSL):** Represents the genetic diversity associated with the genome feature(s) under investigation.

## Screening Process

To initiate screening using the DIGS tool, researchers provide a project-specific command file that serves as the blueprint for the screening process. This command file specifies parameters, including the user-defined name of the screening database, and file paths to the TDb, RSL, and probe sequences.

When a screen is initiated, a project-specific database is created. The core schema of this database can be extended to include any relevant "side data," such as taxonomic information related to the species and sequences included in the screen. This extension increases the power of SQL queries to reveal informative patterns.

Systematic screening proceeds automatically until all searches have been completed. If the process is interrupted at any point or if novel probe/target sequences are incorporated into the project, screening will proceed in a non-redundant way on restarting. Thus, screening projects can readily be expanded to incorporate new TDb files (e.g., recently published WGS assemblies) or novel probe/reference sequences as they become available.

The DIGS tool console allows the reclassification of sequences held in the results table, for example, following an RSL update. To increase efficiency, this process can be tailored to specific subsets of database sequences by supplying SQL constraints via the DIGS tool console.

### Handling Fragmented Matches

BLAST algorithms emphasize local similarity and tend to fragment contiguous matches into several separate hits if similarity across some internal regions of the match is low. The DIGS tool allows screening pipelines to be configured with respect to how overlapping/adjacent hits are handled. This feature enables the screening process to be tailored to the specific needs of diverse projects.

Additionally, the DIGS tool provides a ‘consolidation’ function that concatenates, rather than merges, adjacent hits. The consolidated results, along with information about their structure, are stored in a new screening database table.

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
