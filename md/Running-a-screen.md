**Initial screening run**

Before running a screen for the first time, you will need to index the target database (TDb) for BLAST searching:

```
./digs_tool.pl –m=1 
```

Once the TDb has been prepared, a screen can be executed as follows:

```
./digs_tool.pl –m=2 –i=[path to control file]
```
Progress is written to the terminal, and can also be monitored by issuing SQL queries against the relevant screening database. A screen can be stopped at any time. The next time the tool is restarted, it will initiate screening at the point it left off.

**Running a screen in the background**

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
./digs_tool.pl -i ctl/example_1_EVE.ctl -m=2 > screen.log &
```

To check progress on the backgrounded screen, view the log file. For example using:
```
tail screen.log 
```


**Iterating and refining DIGS**

Exploratory screening is a dynamic discovery process, and it is likely that early screening databases will be superseded by later ones that incorporate more refined probe sets and reference sequence libraries.

For example, an initial screen may lead to the creation of a more refined reference sequence library that incorporates a finer-grained classification of reference sequences, or includes novel references identified through the first round of screening.

To update the contents of the Extracted table to reflect changes in the reference library, execute the digs_tool.pl.pl script as follows.

```
./digs_tool.pl –m=3 –i=[path to control file]
```

**Refreshing and deleting DIGS databases**

The digs_tool.pl script can also be used to 'flush' (empty/reset) and delete databases. Execute digs_tool.pl as follows to flush data a screening database.

```
./digs_tool.pl –m=5 –i=[path to control file]
```

To delete a screening database, execute digs_tool.pl.pl as follows:

```
./digs_tool.pl –m=6 –i=[path to control file]
```


This will drop (delete) the screening database specified in the control file.

