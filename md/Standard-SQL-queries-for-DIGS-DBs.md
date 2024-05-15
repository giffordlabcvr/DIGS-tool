**How many hits in the extracted table?**

```
# How many hits in the extracted table?
SELECT COUNT(*) AS Number_of_hits FROM Extracted;
```

**Which organisms have been screened?**

```
# Which organisms have been screened?
SELECT DISTINCT Organism
FROM Status
ORDER BY Organism;
```

**Summarise DIGS results for each organism**

```
# Summarise DIGS results for each organism
SELECT Organism, Assigned_name, Assigned_gene, COUNT(*) AS Number 
FROM   Extracted 
GROUP BY Organism, Assigned_name, Assigned_gene
ORDER BY Organism, Number DESC;
```

**Create GTF format from Extracted table**

```
# Create GTF format from Extracted table
SELECT 
Scaffold, 
Assigned_gene, 
Extract_start, 
Extract_end,  
'.' as placeholder1, 
Orientation,  
'.' as placeholder2,
CONCAT('gene_id "',Assigned_name,'";') AS gene_id, 
CONCAT('gene_name "',Assigned_name,'";') AS gene_name,
CONCAT('transcript_id "',Assigned_name,'";') AS transcript_id,
CONCAT('tss_id "',Assigned_name,'";') AS tss_id
FROM Extracted
ORDER BY Scaffold, Extract_start;
```

