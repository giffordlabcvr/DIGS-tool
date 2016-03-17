Begin SCREENDB;

	db_name=EVEs;
	mysql_server=localhost;
	mysql_username=root;
	mysql_password=blenat2;

ENDBLOCK;

BEGIN SCREENSETS;

	query_aa_fasta=projects/EVEs/All_EVE_probes.fas;
	reference_aa_fasta=projects/EVEs/All_EVE_references.fas;
	output_path=./tmp/;

	bit_score_min_tblastn=60;
	seq_length_minimum=50;
	threadhit_probe_buffer=100;
	threadhit_gap_buffer=100;
	threadhit_max_gap=100;
	redundancy_mode=2;
	blast_threads=8;

ENDBLOCK;

BEGIN TARGETS;

	Mammalia/Homo_sapiens/complete/goldenpath_hg38
	Mammalia/Pan_troglodytes/complete/goldenpath_panTro4
	Mammalia/Gorilla_gorilla/complete/goldenpath_gorGor3
	Mammalia/Pongo_abelii/complete/goldenpath_ponAbe2
	Mammalia/Nomascus_leucogenys/complete/ensembl_79_Nleu1.0/
	Mammalia/Macaca_mulatta/low_coverage/goldenpath_rheMac3
	Mammalia/Papio_anubis/low_coverage/goldenpath_papAnu2
	Mammalia/Chlorocebus_sabaeus/complete/ensembl_79_ChlSab1.1
	#Mammalia/Cercocebus_atys/low_coverage/ncbi_100_2015
	#Mammalia/Colobus_angolensis_palliatus/low_coverage/ncbi_100_2015
	Mammalia/Nasalis_larvatus/low_coverage/goldepath_nasLar1_2015-12
	Mammalia/Rhinopithecus_roxellana/low_coverage/ncbi_100_2014-11-18
	Mammalia/Callithrix_jacchus/low_coverage/goldenpath_calJac3
	Mammalia/Saimiri_boliviensis/complete/goldenpath_saiBol1
	#Mammalia/Aotus_nancymaae/low_coverage/ncbi_1.0_April_2015
	Mammalia/Tarsius_syrichta/complete/goldenpath_tarSyr1
	Mammalia/Otolemur_garnettii/complete/goldenpath_otoGar3
	Mammalia/Microcebus_murinus/low_coverage/goldenpath_micMur1
	Mammalia/Daubentonia_madagascariensis/low_coverage/broad_1.0_oct_2012
	Mammalia/Galeopterus_variegatus/low_coverage/G-variegatus_3.0.2_apr_2015
	Mammalia/Tupaia_belangeri/low_coverage/goldenpath_tupBel1
	Mammalia/Chinchilla_lanigera/low_coverage/broad_v0_2012
	Mammalia/Heterocephalus_glaber/low_coverage/goldenpath_hetGla2
	Mammalia/Octodon_degus/low_coverage/broad_1.0_july_2012
	Mammalia/Cavia_porcellus/low_coverage/goldenpath_cavPor3
	Mammalia/Rattus_norvegicus/complete/goldenpath_rn6
	Mammalia/Mus_musculus/complete/goldenpath_mm10
	Mammalia/Cricetulus_griseus/complete/goldenpath_criGri1
	Mammalia/Jaculus_jaculus/low_coverage/broad_JacJac1.0
	Mammalia/Mesocricetus_auratus/low_coverage/broad_1.1_april_2013
	Mammalia/Microtus_ochrogaster/complete/broad_1.0_april_2013
	Mammalia/Elephantulus_edwardii/low_coverage/broad_1.0_april_2013
	#Mammalia/Fukomys_damarensis/low_coverage/ncbi_1.0_2014
	Mammalia/Oryctolagus_cuniculus/complete/goldenpath_oryCun2
	Mammalia/Ochotona_princeps/complete/goldenpath_ochPri3
	Mammalia/Tursiops_truncatus/low_coverage/goldenpath_turTru2
	Mammalia/Lipotes_vexillifer/low_coverage/ncbi_100_2014-04-27
	Mammalia/Balaenoptera_acutorostrata/low_coverage/goldenpath_balAcu1
	Mammalia/Orcinus_orca/low_coverage/ncbi_ANOL02_jan_2013
	Mammalia/Capra_hircus/complete/iggc_1.0_april_2013
	Mammalia/Ovis_aires/low_coverage/goldenpath_oviAri3
	Mammalia/Bos_taurus/low_coverage/goldenpath_bosTau7
	Mammalia/Sus_scrofa/complete/goldenpath_susScr3
	Mammalia/Bison_bison/low_coverage/ncbi_UMD1_2015-05-14
	Mammalia/Bos_grunniens/low_coverage/bbu_UMD_2015-05-14
	Mammalia/Bubalus_bubalis/low_coverage/ncbi_100_2015-05-14
	Mammalia/Camelus_ferus/low_coverage/bacsac_1.0_april_2013
	Mammalia/Vicugna_pacos/low_coverage/goldenpath_vicPac2
	Mammalia/Equus_caballus/complete/goldenpath_equCab2_Dec13
	Mammalia/Equus_africanus_asinus/low_coverage/willy_denovo
	Mammalia/Ceratotherium_simum/low_coverage/goldenpath_cerSim1
	Mammalia/Ursus_maritimus/low_coverage/GAJD01_1.0_june_2013
	Mammalia/Ailuropoda_melanoleuca/low_coverage/goldenpath_ailMel1
	Mammalia/Canis_familiaris/low_coverage/goldenpath_canFam3
	Mammalia/Felis_catus/complete/goldenpath_felCat5
	Mammalia/Panthera_tigris_altaica/low_coverage/ncbi_100_2014-03
	Mammalia/Mustela_furo/low_coverage/goldenpath_musFur1
	Mammalia/Odobenus_rosmarus/low_coverage/marinemamm_1.0_april_2013
	Mammalia/Leptonychotes_weddellii/low_coverage/broad_1.1_april_2013
	Mammalia/Megaderma_lyra/low_coverage/ncbi_AWHB00000000.1_2015-05
	Mammalia/Rhinolophus_ferrumequinum/low_coverage/ncbi_AWHA00000000.1_2015-05
	Mammalia/Pteropus_alecto/low_coverage/broad_1.0_april_2012
	Mammalia/Pteropus_vampyrus/low_coverage/goldenpath_pteVam1
	Mammalia/Eidolon_helvum/low_coverage/GCA_000465285.1
	Mammalia/Pteronotus_parnellii/low_coverage/ncbi_AWGZ00000000.1_2015-05
	Mammalia/Eptesicus_fuscus/low_coverage/broad_1.0_april_2013
	Mammalia/Myotis_lucifugus/low_coverage/goldenpath_myoLuc2
	Mammalia/Myotis_davidii/low_coverage/bgi_1.1_april_2013
	Mammalia/Myotis_brandtii/low_coverage/ncbi_100_2013
	Mammalia/Orycteropus_afer/low_coverage/broad_v0_2012
	Mammalia/Echinops_telfairi/complete/goldenpath_echTel2
	Mammalia/Chrysochloris_asiatica/low_coverage/broad_ChrAsi1.0_2012
	Mammalia/Procavia_capensis/low_coverage/goldenpath_proCap1
	Mammalia/Trichechus_manatus/complete/goldenpath_triMan1
	Mammalia/Loxodonta_africana/complete/goldenpath_loxAfr3
	Mammalia/Choloepus_hoffmanni/complete/goldenpath_choHof1
	Mammalia/Dasypus_novemcinctus/complete/goldenpath_dasNov3
	Mammalia/Monodelphis_domestica/complete/goldenpath_monDom5
	Mammalia/Macropus_eugenii/low_coverage/goldenpath_macEug2
	Mammalia/Sarcophilus_harrisii/complete/goldenpath_sarHar1
	Mammalia/Ornithorhynchus_anatinus/low_coverage/goldenpath_ornAna1_2007-03

ENDBLOCK;
