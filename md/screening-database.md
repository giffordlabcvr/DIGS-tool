**DIGS screening database field definitions**

When a screen is initiated by the DIGS tool, a project-specific, relational database is created. This 'screening database' captures data generated by DIGS.

**'searches_performed' table**

Records the details of BLAST searches that have been performed in this DIGS project.

| **Field** | **Type** | **Description** |
|---|---|---|
| record_ID | INT | Automatically incremented primary key |
| probe_ID | VARCHAR | Unique iID of probe sequence |
| probe_name | VARCHAR | Name of probe sequence |
| probe_gene | VARCHAR | Name of probe sequence gene |
| target_id | VARCHAR | Unique identifier the TDb file |
| organism | VARCHAR | Name of the organism (Latin binomial) from which TDb was generated |
| target_datatype | VARCHAR | Data type of the TDb file |
| target_version | VARCHAR | Version details of the TDb file |
| target_name | VARCHAR | Name of TDb file |
| Timestamp | TIMESTAMP |  Timestamp of the table entry |


**'digs_results' table**

Contains the extracted sequences of the loci specified in BLAST_results table, and results of the second round of paired BLAST, in which extracted sequences are 'genotyped' by BLAST comparison to the reference library.

| **Field** | **Type** | **Description** |
|---|---|---|
| record_ID | INT | Automatically incremented primary key |
| organism | VARCHAR | Organism name (Latin binomial) |
| target_datatype | VARCHAR | Genome data type |
| target_version | VARCHAR | Genome build version details |
| target_name | VARCHAR | Name of genome data file containing the BLAST hit |
| probe_type | VARCHAR | Type of probe sequence (amino acid or nucleotide) |
| extract_start | INT | 5’ (start) position of reverse BLAST alignment in the RSL sequence |
| extract_end | INT | 3’ (end) position of reverse BLAST alignment in the RSL sequence  |
| scaffold | VARCHAR | Name of scaffold/contig/chromosome containing the BLAST hit |
| orientation | ENUM | Orientation of the BLAST hit relative to the probe |
| assigned_name | VARCHAR | Name of closest matching sequence in RSL |
| assigned_gene | VARCHAR | Name of gene of closest matching sequence in RSL |
| bit_score | FLOAT | Bit score of the best match from reverse BLAST |
| identity | FLOAT | Percentage identity of the best match from reverse BLAST | 
| e_value_num | FLOAT | Coefficient of the expect (e) value for the best match from reverse BLAST |
| e_value_exp | INT | Exponent (base e) of the expect (e) value for the best match from reverse BLAST |
| subject_start | INT | 5’ (start) position of reverse BLAST alignment in the RSL sequence |
| subject_end | INT | 3’ (end) position of reverse BLAST alignment in the RSL sequence  |
| query_start | INT | 5’ (start) position of reverse BLAST alignment in the probe sequence |
| query_end | INT | 3’ (end) position of reverse BLAST alignment in the probe sequence |
| mismatches | INT | Number of mismatches in alignment from reverse BLAST |
| gap_openings | INT | Number of gap openings in alignment from reverse BLAST |
| sequence_length | INT | Length of the extracted sequence |
| sequence | TEXT | Text string of the extracted sequence | 
| timestamp | TIMESTAMP | Timestamp of the table entry |
