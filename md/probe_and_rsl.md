## Structure of the DIGS Control File

Control files are structured as NEXUS-style blocks delineated by BEGIN and ENDBLOCK tokens.

Below is an annotated example of a DIGS control file for a polypeptide-based screen:

```
BEGIN SCREENDB;
    db_name=[your screening database name];
    mysql_server=[localhost];
ENDBLOCK;

BEGIN SCREENSETS;
    query_aa_fasta=[path to file with amino acid probes];
    reference_aa_fasta=[path to file with amino acid references];
    output_path=[path to output folder];
    bit_score_min_tblastn=100;
    seq_length_minimum=100;
    defragment_range=100;
ENDBLOCK;

BEGIN TARGETS;
    Mammals/Mus_musculus/complete/goldenpath_mm10/
    Mammals/Rattus_norvegicus/complete/goldenpath_rn6/
ENDBLOCK;
```

#### SCREENDB Block

The SCREENDB block defines the database connection details:

- **db_name**: The name of the screening database.
- **mysql_server**: The MySQL server (e.g., 'localhost' for local databases).

Ensure the MySQL user has the necessary privileges: CREATE, DROP, SELECT, ALTER, and DELETE.

#### SCREENSETS Block

The SCREENSETS block sets paths and parameters for BLAST-based screening:

- Paths to probe and reference libraries.
- Thresholds for BLAST significance.

For a complete description of all possible parameters, refer to the DIGS tool wiki.

#### TARGETS Block

The TARGETS block specifies the paths to target files within a directory hierarchy. The DIGS tool will include every FASTA file under the specified path. For example:

```
BEGIN TARGETS;
    Mammals/
ENDBLOCK;
```


**Creating a reference sequence library**

DIGS requires a library of FASTA-formatted sequences to provide (i) **probes** for the first step pf paired BLAST, and (ii) **reference sequences** for the second step of paired BLAST. 

The DIGS tool uses a simple rule to capture data from the headers of FASTA-formatted reference (and probe) sequences; headers should be structured so as to define two hierarchical name elements; ‘name’ and ‘gene_name’, separated by an underscore. In the example shown above these are a virus name and a gene name. Other two-level hierarchical naming schemes (e.g. species & gene name, gene-subdomain name) can also be used, providing the same scheme is used consistently throughout the project. Reference sequences should be stored in a file, the path to which will be specified in the DIGS control file.

**Selecting probe sequences for BLAST-based screening**

A subset of reference sequences should be selected as probes. The entire reference library can be used - in which case there is no need to create a separate file – but it is often sufficient to use only a subset of sequences from the reference sequence library, in which case, a separate file containing this subset should be created. The path to this file is specified in the DIGS [control file](https://github.com/giffordlabcvr/DIGS-tool/wiki/Setting-up-a-control-file).

**NOTE:** DIGS can utilize both polypeptide and nucleotide sequences within the same project and database, but separate probe and reference files should be created for nucleotide and polypeptide sequences, and separate screens using distinct control file parameters should be performed for each. 



This includes all genomes stored under the path $DIGS_GENOMES/Mammals.

To target specific files, specify the complete paths:


```
BEGIN TARGETS;
    Mammals/Homo_sapiens/Complete/ncbi_37.3_june_11/ChrX.fa
    Mammals/Homo_sapiens/Complete/ncbi_37.3_june_11/ChrY.fa
ENDBLOCK;
```
