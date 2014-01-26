BEGIN SCREENDB;
db_name=Parvovirus_EVEs;
mysql_server=localhost;
mysql_username=root;      
mysql_password=blenat2; 
ENDBLOCK;

BEGIN SCREENSETS;
query_aa_fasta=screensets/parvovirus/parvovirus_probes.fa;
reference_aa_fasta=screensets/parvovirus/parvovirus_refseq_library.fa;
bit_score_min_tblastn=100;
seq_length_minimum=100;
ENDBLOCK;

BEGIN TARGETS;
Mammalia/Octodon_degus/
Mammalia/Chinchilla_lanigera/
ENDBLOCK;

BEGIN SCREENSQL;

ENDBLOCK;
