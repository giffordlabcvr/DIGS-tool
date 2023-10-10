**The Database-Integrated Genome Screening (DIGS) Tool, Version 2.0**
------------------------------------------------------------------------------------

The DIGS tool is a program for implementing 'database-integrated genome screening' (DIGS).

Genome screening utilises similarity search tools to systematically search DNA sequence 
databases for sequences related to a one or more query sequences. Queries may be 
nucleic acid or protein sequences. 

The DIGS tool is implemented using the PERL scripting language and is compatible with 
UNIX/LINUX operating systems.

The DIGS tool uses the BLAST+ program suite  to perform similarity searches, and the MySQL 
relational database management system to capture their output.

To run the DIGS tool requires PERL, BLAST+ and MySQL (or a supported fork of MySQL such
as MariaDB).

The program is accessible through a text-based console interface.

To run the DIGS tool and see options type: ./digs_tool.pl

For information, please visit the [DIGS website](http://giffordlabcvr.github.io/DIGS-tool/) 
