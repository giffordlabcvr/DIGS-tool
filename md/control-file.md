**Parameters defined in the DIGS control file**

| Parameter | Definition |
|---|---|
| SCREENDB BLOCK |  |
| db_name | Name of the a screening database for this screening run |
| mysql_server: | MySQL server to use (e.g. 'localhost' if using local database) |
| SCREENSETS BLOCK |  |
| query_aa_fasta* | Path to FASTA file with amino acid probe sequences |
| query_na_fasta* | Path to FASTA file with nucleic acid probe sequences |
| reference_aa_fasta ** | Path to FASTA file with amino acid reference sequences |
| reference_na_fasta ** | Path to FASTA file with nucleic acid reference sequences |
| bit_score_min_tblastn* | Minimum bit score of tBLASTn hit to extract |
| bit_score_min_blastn** | Minimum bit score of BLASTn hit to extract |
| seq_length_minimum | Minimum length of sequence to extract |
| defragment_range| Range within which two BLAST hits in the target sequence will be merged  |

Footnote: * Only required for screens utilizing protein sequences. ** Only required for screens utilizing nucleotide sequences.

The DIGS control file contains parameters and paths for DIGS. Control files are structured as NEXUS style blocks  delineated by BEGIN and ENDBLOCK tokens. An annotated example of a DIGS control file is shown below.  Note that not all parameters will need to be defined for every screen. 

**Example control file for a polypeptide-based screen**

```
Begin SCREENDB;
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

**SCREENDB block**

Block used to define the DB connection details. These include the name of the screening database, and three parameters for connecting to MySQL; (i) server, (ii) user and (iii) password. Note that the user should have CREATE, DROP, SELECT, ALTER, and DELETE privileges


**SCREENSETS block**

This block is used to set paths and parameters for BLAST-based screening. These include the paths to probe and reference libraries, and thresholds for BLAST significance. See [here](https://github.com/giffordlabcvr/DIGS-tool/wiki/Appendices) for a complete description of the parameters that can be specified in this block.

**TARGETS block**

The targets block should contains the paths to the target files in a directory hierarchy as shown on [this page](https://github.com/giffordlabcvr/DIGS-tool/wiki/Creating-a-target-genome-directory). The DIGS tools will include every FASTA file below the specified path as a target. For example, if only the top two levels are specified, as follows:

```
BEGIN TARGETS;
    Mammals/
ENDBLOCK;
```

...all types and versions of all genomes stored under the path '$DIGS_GENOMES/Mammals' would be included in the screen. To target a screen to individual files, write the complete path. For example:

```
BEGIN TARGETS;
    Mammals/Homo_sapiens/Complete/ncbi_37.3_june_11/ChrX.fa
    Mammals/Homo_sapiens/Complete/ncbi_37.3_june_11/ChrY.fa
ENDBLOCK;
```

