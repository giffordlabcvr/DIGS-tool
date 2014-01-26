BEGIN SCREENDB;
db_name=Spumavirus_ERVs;
mysql_server=localhost;
mysql_username=;      
mysql_password=; 
ENDBLOCK;

BEGIN SCREENSETS;
query_aa_fasta=screensets/retrovirus/spumavirus_probes.fa;
reference_aa_fasta=screensets/retrovirus/retrovirus_refseq_library.fa;
bit_score_min_tblastn=100;
seq_length_minimum=100;
ENDBLOCK;

BEGIN TARGETS;
Mammalia/Choloepus_hoffmanni/
ENDBLOCK;

BEGIN SCREENSQL;

ENDBLOCK;
