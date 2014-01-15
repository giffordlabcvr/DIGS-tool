BEGIN PARAMS;
db_name=Spumavirus_ERVs;
query_aa_fasta=db/spumavirus_probes.fa;
reference_aa_fasta=db/retrovirus_refseq_library.fa;
bit_score_min_tblastn=100;
seq_length_minimum=100;
ENDBLOCK;

BEGIN TARGETS;
Mammalia/Choloepus_hoffmanni/
ENDBLOCK;
