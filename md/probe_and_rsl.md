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

This includes all genomes stored under the path $DIGS_GENOMES/Mammals.

To target specific files, specify the complete paths:


```
BEGIN TARGETS;
    Mammals/Homo_sapiens/Complete/ncbi_37.3_june_11/ChrX.fa
    Mammals/Homo_sapiens/Complete/ncbi_37.3_june_11/ChrY.fa
ENDBLOCK;
```
