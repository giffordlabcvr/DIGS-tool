BEGIN PARAMS;
query_aa_fasta=db/parvovirus_probes.fa;
reference_aa_fasta=db/parvovirus_refseq_library.fa;
bit_score_min_tblastn=100;
seq_length_minimum=100;
ENDBLOCK;

BEGIN TARGETS;
Mammalia/Octodon_degus/
Mammalia/Chinchilla_lanigera/
ENDBLOCK;
