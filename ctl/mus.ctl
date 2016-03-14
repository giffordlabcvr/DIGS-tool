Begin SCREENDB;
	db_name=Mus_ERVs;
	mysql_server=localhost;
	mysql_username=root;
	mysql_password=blenat2;
ENDBLOCK;

BEGIN SCREENSETS;
	query_aa_fasta=/Users/admdart/DIGS-tool/projects/Mus_RT/mus-rt-probe.fas;
	reference_aa_fasta=/Users/admdart/DIGS-tool/projects/Mus_RT/mus-rt-probe.fas;
	output_path=./tmp/;
	bit_score_min_tblastn=100;
	seq_length_minimum=100;
	redundancy_mode=2;
	threadhit_probe_buffer=100;
	threadhit_gap_buffer=100;
	threadhit_max_gap=100;
ENDBLOCK;

BEGIN TARGETS;
Mammals/Mus_musculus/complete/goldenpath_mm10
ENDBLOCK;
