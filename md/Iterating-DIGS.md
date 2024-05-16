**Iterating and refining DIGS**

Exploratory screening is a dynamic discovery process, and it is likely that early screening databases will be superseded by later ones that incorporate more refined probe sets and reference sequence libraries.

For example, an initial screen may lead to the creation of a more refined reference sequence library that incorporates a finer-grained classification of reference sequences, or includes novel references identified through the first round of screening.

To update the contents of the Extracted table to reflect changes in the reference library, execute the pipeline.pl script as follows.

```
./pipeline –m=3 –i=[path to control file]
```

**Refreshing and deleting DIGS databases**

The pipeline.pl script can also be used to 'flush' (empty/reset) and delete databases. Execute pipeline.pl as follows to flush data a screening database.

```
./pipeline –m=5 –i=[path to control file]
```

To delete a screening database, execute pipeline.pl as follows:

```
./pipeline –m=6 –i=[path to control file]
```


This will drop (delete) the screening database specified in the control file.
