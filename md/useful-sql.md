# Database querying examples

**How many rows in the digs_results table?**

```
SELECT COUNT(*) AS number_of_hits FROM digs_results;
```

**Which organisms have been screened?**

```
SELECT DISTINCT organism
FROM searches_performed
ORDER BY organism;
```

**Summarise DIGS results for each organism**

```
SELECT DISTINCT organism, assigned_name, assigned_gene, COUNT(*) AS Number 
FROM digs_results 
GROUP BY organism, assigned_name, assigned_gene
ORDER BY assigned_name, number DESC;
```

**Create GTF format from digs_results table**

```
SELECT scaffold, assigned_gene, extract_start, extract_end,  
'.' as placeholder1, 
orientation,  
'.' as placeholder2,
CONCAT('gene_id "',assigned_name,'";') AS gene_id, 
CONCAT('gene_name "',assigned_name,'";') AS gene_name,
CONCAT('transcript_id "',assigned_name,'";') AS transcript_id,
CONCAT('tss_id "',assigned_name,'";') AS tss_id
FROM digs_results
ORDER BY scaffold, extract_start;
```

**Create basic BED format from digs_results table**

```
SELECT scaffold, extract_start, extract_end, assigned_gene, '.' as placeholder1, orientation
FROM results
ORDER BY scaffold, extract_start
```

