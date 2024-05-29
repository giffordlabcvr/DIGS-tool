## Parameters Defined in the DIGS Control File

The DIGS control file contains parameters and paths for DIGS. Not all parameters are necessary for every screening configuration. 

| Parameter | Definition |
|---|---|
| SCREENDB BLOCK |  |
| db_name | Name of the screening database for this screening run |
| mysql_server: | MySQL server to use (e.g. 'localhost' for a local database) |
| SCREENSETS BLOCK |  |
| query_aa_fasta* | Path to FASTA file with amino acid probe sequences |
| query_na_fasta* | Path to FASTA file with nucleic acid probe sequences |
| reference_aa_fasta ** | Path to FASTA file with amino acid reference sequences |
| reference_na_fasta ** | Path to FASTA file with nucleic acid reference sequences |
| bit_score_min_tblastn* | Minimum bit score for recording hits (tBLASTn) |
| bit_score_min_blastn** | Minimum bit score for recording hits (BLASTn) |
| seq_length_minimum | Minimum length (nucleotides) for recording hits  |
| defragment_range| Range within which two BLAST hits in the target sequence will be merged  |

Footnote: * Only required for screens utilizing protein sequences. ** Only required for screens utilizing nucleotide sequences.

