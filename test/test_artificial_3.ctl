# Test - screen with a short probe, then a long one

Begin SCREENDB;
    db_name=digs_artifical;
    mysql_server=localhost;
    mysql_username=digs;    
    mysql_password=digs12345; 
ENDBLOCK;

BEGIN SCREENSETS;
    #query_aa_fasta=./test/probes/Parvo-NS-probe.faa;
    query_aa_fasta=./test/probes/Chinchilla_NS-edited-probe.faa
    reference_aa_fasta=./test/probes/Parvo-NS-probe.faa;
    output_path=./tmp/;
    bit_score_min_tblastn=40;
    seq_length_minimum=50;
    defragment_range=100;
    redundancy_mode=2;
    blast_threads=8;
    #extract_buffer=1000;
ENDBLOCK;

BEGIN TARGETS;
	test/fake/fake/fake/fake4.fa.fasta
ENDBLOCK;
