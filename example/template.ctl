Begin SCREENDB;
	db_name=example_db;
	mysql_server=localhost;
ENDBLOCK;

BEGIN SCREENSETS;
	output_path=./example/tmp/;
	query_aa_fasta=/example/libraries/example-probe-library.faa;
	reference_aa_fasta=/example/libraries/example-reference-library.faa;
	seq_length_minimum=60;
	bitscore_min_tblastn=40;
	defragment_range=100;
ENDBLOCK;

BEGIN TARGETS;
	./example/targets/;
ENDBLOCK;
