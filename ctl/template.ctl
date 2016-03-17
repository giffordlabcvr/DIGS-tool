Begin SCREENDB;
    db_name=[your screening database name];
    mysql_server=[localhost];
	mysql_username=[your username];
	mysql_password=[your password];
ENDBLOCK;

BEGIN SCREENSETS;
	query_aa_fasta=[path to probes];
	reference_aa_fasta=[path to references];
	output_path=tmp/;
	redundancy_mode=2;
	seq_length_minimum=30;
	bit_score_min_tblastn=100;
	threadhit_probe_buffer=100;
	threadhit_gap_buffer=100;
	threadhit_max_gap=100;
ENDBLOCK;

BEGIN TARGETS;
Mammals/Mus_musculus/complete/goldenpath_mm10
ENDBLOCK;
