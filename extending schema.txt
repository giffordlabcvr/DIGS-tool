The addition of new tables to the screening database allows additional data to be leveraged in SQL queries. Examples of the kinds of tables that might be added include include extra information about reference sequences, genomes, and chromosomes. The DIGS tool includes a utility function for adding new tables (which can also be done via a MySQL client). In the example below, we use this utility upload a table that includes detailed taxonomy of sequences included in the reference library.

```
rob$ ./digs_tool.pl -i=ctl/eves.ctl -d=1

	######################################################################
	#                                                                    #
	#                              DIGS 1.13                             #
	#                Database-Integrated Genome Screening                #
	#                         Robert J. Gifford                          #
	#                   <robert.gifford@glasgow.ac.uk>                   #
	#                                                                    #
	######################################################################

	  Connecting to DB:  EVEs

	 #### WARNING: This function expects a tab-delimited data table with column headers!

	 Please enter the path to the file with the table data and column headings

	 : projects/EVES/All_EVE_references.txt


	 The following column headers (i.e. table fields) were obtained

		 Column 1: 'Superfamily'
		 Column 2: 'State'
		 Column 3: 'Family'
		 Column 4: 'Genus'
		 Column 5: 'Name'
		 Column 6: 'Gene'

	 Is this correct? (y/n): y

		 1. Create new ancillary table
		 2. Append data to existing ancillary table
		 3. Flush existing ancillary table and upload fresh data
		 4. Drop an ancillary table

	 Choose an option (1/2/3/4): 1

	 What is the name of the new table? : eve_data

	 Creating ancillary table 'eve_data' in EVEs screening database
	 Row count 1: uploading value 'Mononegavirales' to field 'Superfamily'
	 Row count 1: uploading value 'Exogenous' to field 'State'
	 Row count 1: uploading value 'Bornavirus' to field 'Family'
	 Row count 1: uploading value 'Bornavirus' to field 'Genus'
	 Row count 2: uploading value 'ABV' to field 'Name'
	 Row count 2: uploading value 'M' to field 'Gene'
         ...
```
This allows us to execute an SQL query that summarises hits based on data in the new table, for example:

```
SELECT Genus, COUNT(*) AS number 
FROM digs_results, eve_data
WHERE digs_results.assigned_name = eve_data.name
GROUP BY Genus
```
...produces output as follows:

```
+-----------------+--------+
| Genus           | Number |
+-----------------+--------+
| Amdovirus       |      8 |
| Bornavirus      |    406 |
| Dependovirus    |    682 |
| Ebolavirus      |     88 |
| Marburgvirus    |     30 |
| Protoparvovirus |     44 |
+-----------------+--------+
6 rows in set (0.01 sec)
```
