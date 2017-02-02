Begin SCREENDB;
    db_name=digs_test_screen;
    mysql_server=localhost;
	mysql_username=root;
	mysql_password=blenat2;
ENDBLOCK;

BEGIN SCREENSETS;
    query_aa_fasta=./test/korv_test2.faa;
    reference_aa_fasta=./test/korv_test2.faa;
    bitscore_min_tblastn=100;
    output_path=./../tmp/;
    seq_length_minimum=50;
    defragment_range=100;
    redundancy_mode=2;
    blast_threads=8;
ENDBLOCK;

BEGIN TARGETS;
	test/fake_species/fake_datatype/fake_version/artificial_test1_korv.fa
ENDBLOCK;
