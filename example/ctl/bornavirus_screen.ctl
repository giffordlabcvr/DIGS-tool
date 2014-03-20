BEGIN SCREENDB;
db_name=Bornavirus_EVEs;
mysql_server=localhost;
mysql_username=root;      
mysql_password=blenat2; 
ENDBLOCK;

BEGIN SCREENSETS;
query_aa_fasta=example/screensets/bornavirus/bornavirus_probes.fa;
#reference_aa_fasta=example/screensets/bornavirus/bornavirus_refseq_library_update.fa;
reference_aa_fasta=example/screensets/bornavirus/bornavirus_refseq_library_.fa;
bit_score_min_tblastn=100;
seq_length_minimum=100;
ENDBLOCK;

BEGIN TARGETS;
Mammalia/Homo_sapiens/
ENDBLOCK;

BEGIN SCREENSQL;
select_list:Organism, assigned_name, assigned_gene, Bit_score, Scaffold, Target_name, Subject_start, Subject_end, Sequence;
where_statement:WHERE Organism="Homo_sapiens" ORDER BY Assigned_name, Bit_score DESC;
ENDBLOCK;
