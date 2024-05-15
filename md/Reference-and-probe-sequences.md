**Creating a reference sequence library**

DIGS requires a library of FASTA-formatted sequences to provide (i) **probes** for the first step pf paired BLAST, and (ii) **reference sequences** for the second step of paired BLAST. 

[[images/FASTA-format.png]]

The DIGS tool uses a simple rule to capture data from the headers of FASTA-formatted reference (and probe) sequences; headers should be structured so as to define two hierarchical name elements; ‘name’ and ‘gene_name’, separated by an underscore. In the example shown above these are a virus name and a gene name. Other two-level hierarchical naming schemes (e.g. species & gene name, gene-subdomain name) can also be used, providing the same scheme is used consistently throughout the project. Reference sequences should be stored in a file, the path to which will be specified in the DIGS [control file](https://github.com/giffordlabcvr/DIGS-tool/wiki/Setting-up-a-control-file).

**Selecting probe sequences for BLAST-based screening**

A subset of reference sequences should be selected as probes. The entire reference library can be used - in which case there is no need to create a separate file – but it is often sufficient to use only a subset of sequences from the reference sequence library, in which case, a separate file containing this subset should be created. The path to this file is specified in the DIGS [control file](https://github.com/giffordlabcvr/DIGS-tool/wiki/Setting-up-a-control-file).

**NOTE:** DIGS can utilize both polypeptide and nucleotide sequences within the same project and database, but separate probe and reference files should be created for nucleotide and polypeptide sequences, and separate screens using distinct control file parameters should be performed for each. 
