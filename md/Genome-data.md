**Organizing genome data for screening**

The DIGS tool is designed for use with locally stored, FASTA formatted sequence data. Genome data should be stored in a directory tree with a specific structure. 

[[images/TargetFileDirectory.png]]

The 'target genomes directory' should have five levels, as illustrated above. The top directory level should be the path specified under the DIGS environment variable ($DIGS_GENOMES). The top two directory levels can be named anything. 

The bottom three levels should contain directories with names pertaining to the data type (i.e. low coverage, assembly), version, and organism name (Latin binomial with underscore – e.g. Homo_sapiens), of the sequence data files they contain*.

In the above example, the file 'chrX.fa' would be under the path:
```
$DIGS_GENOMES/Mammals/Homo_sapiens/Complete/ncbi_37.3_june_11/ChrX.fa
```

**Indexing genome data for BLAST**

BLAST requires that FASTA files are indexed for similarity searches using the 'makeblastdb' program that is distributed with the BLAST+ package. Because this can be time-consuming when screening many separate files, the digs_tool.pl.pl script includes an option to automatically format target FASTA files. Note that for this to work properly, FASTA files MUST be labeled with the appropriate file extensions (.fa, .fas, or .fasta).

To run the DIGS tool's genome formatting utility, execute the digs_tool.pl.pl script as follows:

```
./digs_tool.pl –m=8
```

**Some common FTP sources for genome sequence data**

| Name | URL |
|---|---|
|NCBI genomes |	ftp://bio-mirror.net/biomirror/ncbigenomes/ |
| Broad Institute genomes | http://www.broadinstitute.org/ftp/pub/assemblies/ |
| Ensembl genomes	| ftp://ftp.ensembl.org/pub/release-67/fasta/ |
| UCSC Goldenpath	| http://hgdownload.soe.ucsc.edu/goldenPath/ |





