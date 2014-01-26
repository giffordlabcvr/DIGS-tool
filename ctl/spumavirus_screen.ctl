BEGIN SCREENDB;
db_name=Spumavirus_ERVs;
ENDBLOCK;

BEGIN SCREENSETS;
query_aa_fasta=screensets/spumavirus_probes.fa;
reference_aa_fasta=screensets/retrovirus_refseq_library.fa;
bit_score_min_tblastn=100;
seq_length_minimum=100;
ENDBLOCK;

BEGIN TARGETS;
Mammalia/Choloepus_hoffmanni/
ENDBLOCK;
