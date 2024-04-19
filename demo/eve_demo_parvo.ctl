Begin SCREENDB;
	db_name=rg_eve_demo_parvo;
	mysql_server=localhost;
ENDBLOCK;

BEGIN SCREENSETS;
	query_aa_fasta=/home2/rg128p/DIGS/DIGS-tool/demo/Parvoviridae.DIGS.faa;
	reference_aa_fasta=/home2/rg128p/DIGS/DIGS-tool/demo/NCBI_viruses_none-missing.faa;
	output_path=./tmp/;
	bitscore_min_tblastn=60;
	seq_length_minimum=40;
	defragment_range=1000;
ENDBLOCK;

BEGIN TARGETS;
	Mammalia/
ENDBLOCK;

