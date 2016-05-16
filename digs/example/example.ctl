Begin SCREENDB;
    db_name=[your screening database name];
    mysql_server=[localhost];
	mysql_username=[your username];
	mysql_password=[your password];
ENDBLOCK;

BEGIN SCREENSETS;
	query_aa_fasta=./digs/example/herv-k-hml2_rt_probe.fa;
	reference_aa_fasta=./digs/example/erv_rt_references.fa
	output_path=./digs/tmp/;
	redundancy_mode=2;
	seq_length_minimum=30;
	bit_score_min_tblastn=100;
	threadhit_probe_buffer=100;
	threadhit_gap_buffer=100;
	threadhit_max_gap=100;
ENDBLOCK;

BEGIN TARGETS;
Homo_sapiens/complete/goldenpath_hg19/chrY.fa
ENDBLOCK;
