**Parameters defined in the DIGS control file**

| Parameter | Definition |
|---|---|
| SCREENDB BLOCK |  |
| db_name | Name of the a screening database for this screening run |
| mysql_server: | MySQL server to use (e.g. 'localhost' if using local database) |
| mysql_username | Name of user with CREATE, DROP, SELECT, ALTER, DELETE privileges |
| mysql_password: | Password of the above user |
| SCREENSETS BLOCK |  |
| query_aa_fasta* | Path to FASTA file with amino acid probe sequences |
| query_na_fasta* | Path to FASTA file with nucleic acid probe sequences |
| reference_aa_fasta ** | Path to FASTA file with amino acid reference sequences |
| reference_na_fasta ** | Path to FASTA file with nucleic acid reference sequences |
| bit_score_min_tblastn* | Minimum bit score of tBLASTn hit to extract |
| bit_score_min_blastn** | Minimum bit score of BLASTn hit to extract |
| seq_length_minimum | Minimum length of sequence to extract |
| redundancy mode | A flag for setting how [overlapping, redundant and fragmented sequence hits are handled](https://github.com/giffordlabcvr/DIGS-tool/wiki/Redundant,-overlapping-&-fragmented-hits) |
| defragment_range| Range within which two BLAST hits in the target sequence will be merged  |

Footnote: * Only required for screens utilizing protein sequences. ** Only required for screens utilizing nucleotide sequences.
