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
