BEGIN SCREENDB;
db_name=Bornavirus_EVE_example;
mysql_server=localhost;
mysql_username=[your username here];      
mysql_password=[your password here]; 
ENDBLOCK;

BEGIN SCREENSETS;
query_aa_fasta=example/screensets/bornavirus/bornavirus_probes.fa;
reference_aa_fasta=example/screensets/bornavirus/bornavirus_refseq_library.fa;
#reference_aa_fasta=example/screensets/bornavirus/bornavirus_refseq_library_update.fa;
bit_score_min_tblastn=100;
seq_length_minimum=100;
ENDBLOCK;

BEGIN TARGETS;
#Mammalia/Homo_sapiens/
Mammalia/Homo_sapiens/complete/ensembl_hg19/chr9.fa
ENDBLOCK;

BEGIN SCREENSQL;
select_list:Organism, assigned_name, assigned_gene, Bit_score, Scaffold, Target_name, Subject_start, Subject_end, Sequence;
where_statement:WHERE Organism="Homo_sapiens" ORDER BY Assigned_name, Bit_score DESC;
ENDBLOCK;
