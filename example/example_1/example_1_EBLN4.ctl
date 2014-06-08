BEGIN SCREENDB;
db_name=EBLN_example;
mysql_server=localhost;
#mysql_username=[your username here];      
#mysql_password=[your password here]; 
mysql_username=root;      
mysql_password=MyCVR!2; 
ENDBLOCK;

BEGIN SCREENSETS;
query_aa_fasta=example/example_1/bornavirus_probes.fa;
#reference_aa_fasta=example/example_1/bornavirus_refseq_library.fa;
reference_aa_fasta=example/example_1/bornavirus_refseq_library_update.fa;
bit_score_min_tblastn=100;
seq_length_minimum=100;
redundancy_mode=1;
threadhit_probe_buffer=1000;
threadhit_gap_buffer=1000;
threadhit_max_gap=1000;
ENDBLOCK;

BEGIN TARGETS;
#Mammalia/Homo_sapiens/
Mammalia/Homo_sapiens/complete/ensembl_hg19/chr17.fa
ENDBLOCK;

BEGIN SCREENSQL;
select_list:Organism, assigned_name, assigned_gene, Bit_score, Scaffold, Target_name, Subject_start, Subject_end, Sequence;
where_statement:WHERE Organism="Homo_sapiens" ORDER BY Assigned_name, Bit_score DESC;
ENDBLOCK;
