## Features

- **Sequence Similarity Search**: Integrate tools like BLAST to find related nucleotide sequences.
- **Relational Database Integration**: Store and manage large sets of sequence data efficiently.
- **Data Analysis**: Utilize a range of tools to analyze and visualize genomic data.
- **Comparative Genomics**: Perform in-depth comparative studies across different genomes.


**Using sequence similarity searches and strategically chosen reference datasets to explore genomic 'dark matter'**

DNA sequencing has advanced dramatically over recent years, leading to exponential increases in the volume of sequence data available for analysis. There are consequently unprecedented opportunities to generate new discoveries by by using 'sequence similarity search' tools in combination with strategically chosen reference datasets. The figure below illustrates the general principles of these approaches, drawing a comparison with the approaches used to characterize other kinds of biological diversity. 

[[images/ExploringGenomicLandscape.png]]

In this analogy, poorly characterized genomic diversity is represented by the flora and fauna of a newly encountered land. The approach to characterizing this diversity is comparative. How similar is the new diversity to the diversity we have previously described?

**Database-Integrated Genome Screening (DIGS)**

Similarity searches enable researchers to survey the immense landscape of sequence diversity and selectively recover similar – thus potentially related – sequences. In database-integrated genome-screening (DIGS), the output of similarity search-based genome 'screens' is captured in a relational database. This facilitates the implementation of automated screens that can be performed on a large scale, and allows for the interrogation and manipulation of output data using structured query language (SQL).

The DIGS tool uses the *basic local alignment search tool* ([BLAST](http://blast.ncbi.nlm.nih.gov/Blast.cgi)) to perform sequence similarity searches. BLAST takes a nucleotide or polypeptide sequence as query input and efficiently searches target databases for similar sequences. Results are reported in the form of a ranked list of ‘hits’, each associated with a statistical significance. 

**Use of 'reciprocal BLAST searches' in DIGS** 

[[images/DIGS-screen.png]]

The DIGS tool utilizes a screening approach based on two rounds of BLAST, a strategy we refer to as 'paired BLAST'. In the first round, query sequences selected from the reference sequence library (the 'probes') are used to search target databases. In the second, sequence ‘hits’ identified by screening are extracted from genomes and assigned a genotype by BLAST comparison to the reference sequence set.

The second BLAST step is included because probe sequences can often cross-match to a wide range of homologous sequences in the initial BLAST screen. For example, consider a gene that has two paralogs, ‘X’ and ‘Y’.  Screening with a probe of type X may yield hits to both X and Y. Comparing hits to a library of representative reference sequences in the second BLAST step provides an efficient means for users to discriminate which hits are more like X’s and which are more like Y’s.


[[images/DIGS-screen-generic22.png]]


**DIGS workflows**
	
The DIGS tool is designed to support exploratory, comparative analyses of genomes. These projects have a fundamentally heuristic nature, and consequently will commonly require iteration. An iterative DIGS workflow can be conceptually divided into two stages. In the first stage, paired BLAST screening is performed. In the second, information recovered via screening is examined using the screening database (i.e. using structured query language (SQL)) or external applications (e.g. phylogenetic packages, statistical analysis software). Following this characterization, new information can be incorporated into reference libraries and/or screening databases, and used to perform a more powerful or precisely-targeted screen in the next iteration. 



**Initial screening run**

Before running a screen for the first time, create the screening database (a one-off procedure), as follows:

```
./pipeline –m=1 –i=[path to control file]
```
Once the screening database is successfully created, a screen can be executed as follows:

```
./pipeline –m=2 –i=[path to control file]
```
Progress is written to the terminal, and can also be monitored by issuing SQL queries against the relevant screening database. A screen can be stopped at any time. The next time the tool is restarted, it will initiate screening at the point it left off.

**Running a screen in the background**

Since similarity search-based screens can take hours or days to complete, you may want to run your screen 'in the background', particularly if you are running DIGS on a networked server that you connect to remotely. This can be done as follows. First start the pipeline as follows:

```
./pipeline –m=2 –i=[path to control file] > [path to a log file, e.g. screen.log]
```

Once this command has been executed, press Ctrl+Z and then type:

```
bg
```
You should then see this output:

```
./pipeline.pl -i ctl/example_1_EVE.ctl -m=2 > screen.log &
```

To check progress on the backgrounded screen, view the log file. For example using:
```
tail screen.log 
```



