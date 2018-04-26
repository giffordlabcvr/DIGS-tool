Begin SCREENDB;
	db_name=example_db;
	mysql_server=localhost;
ENDBLOCK;

BEGIN SCREENSETS;
	output_path=./example/tmp/;
	query_aa_fasta=./example/libraries/erv_rt_probes.faa;
	reference_aa_fasta=./example/libraries/erv_rt_references.faa;
	seq_length_minimum=60;
	bitscore_min_tblastn=40;
	defragment_range=100;
ENDBLOCK;

BEGIN TARGETS;
	Mammalia/Homo_sapiens/complete/goldenpath_hg19/chrY.fa
ENDBLOCK;
