# Running a screen 


## Input Data Components

1. **Target Database (TDb):** A collection of whole genome sequence or transcriptome assemblies serving as the target for similarity searches.
2. **Query Sequences (Probes):** Input sequences for similarity searches of the Target Database.
3. **Reference Sequence Library (RSL):** Represents the genetic diversity associated with the genome feature(s) under investigation.

## Screening Process

To initiate screening using the DIGS tool, researchers provide a project-specific command file that serves as the blueprint for the screening process. This command file specifies parameters, including the user-defined name of the screening database, and file paths to the TDb, RSL, and probe sequences.

When a screen is initiated, a project-specific database is created. The core schema of this database can be extended to include any relevant "side data," such as taxonomic information related to the species and sequences included in the screen. This extension increases the power of SQL queries to reveal informative patterns.

Systematic screening proceeds automatically until all searches have been completed. If the process is interrupted at any point or if novel probe/target sequences are incorporated into the project, screening will proceed in a non-redundant way on restarting. Thus, screening projects can readily be expanded to incorporate new TDb files (e.g., recently published WGS assemblies) or novel probe/reference sequences as they become available.

The DIGS tool console allows the reclassification of sequences held in the results table, for example, following an RSL update. To increase efficiency, this process can be tailored to specific subsets of database sequences by supplying SQL constraints via the DIGS tool console.


## Initial screening run

Before running a screen for the first time, you will need to index the target database (TDb) for BLAST searching:

```
./digs_tool.pl –m=1 –i=[path to control file]
```
Once the screening database is successfully created, a screen can be executed as follows:

```
./digs_tool.pl –m=2 –i=[path to control file]
```
Progress is written to the terminal, and can also be monitored by issuing SQL queries against the relevant screening database. A screen can be stopped at any time. The next time the tool is restarted, it will initiate screening at the point it left off.

## Running a screen in the background

Since similarity search-based screens can take hours or days to complete, you may want to run your screen 'in the background', particularly if you are running DIGS on a networked server that you connect to remotely. This can be done as follows. First start the digs_tool.pl as follows:

```
./digs_tool.pl –m=2 –i=[path to control file] > [path to a log file, e.g. screen.log]
```

Once this command has been executed, press Ctrl+Z and then type:

```
bg
```

You should then see this output:


```
./digs_tool.pl.pl -i ctl/example_1_EVE.ctl -m=2 > screen.log &
```

To check progress on the backgrounded screen, view the log file. For example using:
```
tail screen.log 
```


