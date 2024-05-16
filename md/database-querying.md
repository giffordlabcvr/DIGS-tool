DIGS expedites the interrogation of data derived from BLAST-based screening. Data can be interrogated within the framework of the DIGS tool (i.e. by using structured query language (SQL)), or exported to external applications and examined therein.

**Using SQL to interrogate DIGS output**

The example below shows one way in which the output of DIGS can be summarized using SQL.

```
SELECT DISTINCT Assigned_name, COUNT(*) AS Number 
FROM Extracted 
GROUP BY Assigned_name
ORDER BY Number DESC;
```
This command, when entered into a MySQL client, should produce output similar to that shown below:
```
+---------------+--------+
| Assigned_name | Number |
+---------------+--------+
| BDV           |    156 |
| MARV          |     80 |
| Caprine_AAV   |     42 |
| Avian_AAV     |     10 |
| Bovine_AAV    |      8 |
| AAV2          |      6 |
| AMDV          |      4 |
+---------------+--------+
7 rows in set (0.00 sec)
```

**Analyzing DIGS output in external applications**

The extracted sequences can be obtained from the Extracted table using an SQL query.
```
SELECT Record_ID, Sequence FROM Extracted
```
To convert from a two column text tab-delimited file (ID and sequence) to a fasta file, you can use the following command:

```
awk -vOFS='' '{print ">",$1,"\n",$2,"\n";}' two_column_sample_tab.txt > sample1.fa
```