# Extending the Database Schema of a DIGS Project

Adding new tables to the screening database allows you to leverage additional data in SQL queries. Examples of such tables might include extra information about reference sequences, genomes, and chromosomes. The DIGS tool includes a utility function for adding new tables. This can also be done via a MySQL client, but the console-based functions provided by the DIGS tool may be useful for quick management of DIGS screening databases from the terminal.

In the example below, we use this utility to add a table that includes detailed taxonomy of sequences included in the reference library.

## Steps to Upload a New Table

1. **Run the DIGS Tool**: Open your terminal and execute the following command:
    ```bash
    ./digs_tool.pl -i=path/to/eve_screen_file.ctl -d=1
    ```

2. **Console Output**: You will see output similar to the following:
    ```
    ######################################################################
    #                              DIGS 1.13                             #
    #                Database-Integrated Genome Screening                #
    #                         Robert J. Gifford                          #
    #                   <robert.gifford@glasgow.ac.uk>                   #
    ######################################################################

      Connecting to DB:  EVEs

    #### WARNING: This function expects a tab-delimited data table with column headers!

    Please enter the path to the file with the table data and column headings:
    ```

3. **Input the File Path**: Enter the path to your data file, for example:
    ```plaintext
    projects/EVES/All_EVE_references.txt
    ```

4. **Verify Column Headers**: The tool will display the detected column headers:
    ```
    The following column headers (i.e. table fields) were obtained:
        Column 1: 'Superfamily'
        Column 2: 'State'
        Column 3: 'Family'
        Column 4: 'Genus'
        Column 5: 'Name'
        Column 6: 'Gene'

    Is this correct? (y/n):
    ```

5. **Confirm Headers**: Type `y` to confirm the headers.

6. **Choose an Option**: You will then be prompted to choose an action:
    ```
    1. Create new ancillary table
    2. Append data to existing ancillary table
    3. Flush existing ancillary table and upload fresh data
    4. Drop an ancillary table

    Choose an option (1/2/3/4):
    ```

7. **Create New Table**: Select option `1` and enter the name of the new table:
    ```plaintext
    What is the name of the new table? : eve_data
    ```

8. **Table Creation and Data Upload**: The tool will create the table and upload data:
    ```
    Creating ancillary table 'eve_data' in EVEs screening database
    Row count 1: uploading value 'Mononegavirales' to field 'Superfamily'
    Row count 1: uploading value 'Exogenous' to field 'State'
    Row count 1: uploading value 'Bornavirus' to field 'Family'
    Row count 1: uploading value 'Bornavirus' to field 'Genus'
    Row count 2: uploading value 'ABV' to field 'Name'
    Row count 2: uploading value 'M' to field 'Gene'
    ...
    ```

## Example SQL Query

Once the table is uploaded, you can execute SQL queries that leverage the new data. For example, the following query summarizes hits based on data in the new table:

```
SELECT Genus, COUNT(*) AS number 
FROM digs_results, eve_data
WHERE digs_results.assigned_name = eve_data.name
GROUP BY Genus;
```

## Example Output

The query produces output similar to the following:

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

This process allows you to extend the database schema and utilize additional data in your analyses.

