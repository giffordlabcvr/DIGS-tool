Begin SCREENDB;
db_name=HERVRT_example;
mysql_server=localhost;
mysql_username=root;
mysql_password=blenat2;
#mysql_username=[your username here];      
#mysql_password=[your password here]; 
ENDBLOCK;

BEGIN SCREENSETS;
query_aa_fasta=example/example_3/retrovirus_probes.fa;
reference_aa_fasta=example/example_3/retrovirus_RT_references_AA.fas;
bit_score_min_tblastn=100;
seq_length_minimum=100;
redundancy_mode=2;
threadhit_probe_buffer=1000;
threadhit_gap_buffer=1000;
threadhit_max_gap=1000;
ENDBLOCK;

BEGIN TARGETS;
#Mammalia/Homo_sapiens/
Mammalia/Homo_sapiens/complete/ensembl_hg19/chr17.fa
ENDBLOCK;

#BEGIN SCREENSQL;
#select_list:Organism, assigned_name, assigned_gene, Bit_score, Scaffold, Target_name, Subject_start, Subject_end, Sequence;
#where_statement:WHERE Organism="Homo_sapiens" ORDER BY Assigned_name, Bit_score DESC;
#ENDBLOCK;
