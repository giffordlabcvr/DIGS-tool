Begin SCREENDB;
    db_name=digs;
    mysql_server=localhost;
    mysql_username=root;    
    mysql_password=digs12345; 
ENDBLOCK;

BEGIN SCREENSETS;
    query_aa_fasta=./test/probes/Hepadna-Pol-probe.faa;
    reference_aa_fasta=./test/probes/Hepadna-Pol-probe.faa;
    output_path=./screens/tmp/;
    bit_score_min_tblastn=60;
    seq_length_minimum=50;
    defragment_range=1000;
    redundancy_mode=2;
    blast_threads=8;
ENDBLOCK;

BEGIN TARGETS;
	Mammalia/Mus_musculus/complete/goldenpath_mm10_sex/
ENDBLOCK;

BEGIN SKIPINDEX;
	test/Locusta_migratoria/low_coverage/LocustGenomeV1/
ENDBLOCK;
