BEGIN SCREENDB;
db_name=Bornavirus_EVEs;
ENDBLOCK;

BEGIN SCREENSETS;
query_aa_fasta=db/bornavirus_probes.fa;
reference_aa_fasta=db/bornavirus_refseq_library.fa;
bit_score_min_tblastn=100;
seq_length_minimum=100;
ENDBLOCK;

BEGIN TARGETS;
Mammalia/Homo_sapiens/
Mammalia/Pan_troglodytes/
ENDBLOCK;
