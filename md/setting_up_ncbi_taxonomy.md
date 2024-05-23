# Creating a local copy of the NCBI taxonomy database

This document describes how to set up a copy of the NCBI taxonomy database in MySQL.

## Getting the data

Download the taxonomy DB from the NCBI FTP server.


```
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
```

Unpack it:

```
gunzip taxdump.tar.gz
tar -xvf taxdump.tar
```

## Creating the database tables in MySQL

Load the schema into a MySQL database.

```
CREATE DATABASE NCBI_Taxonomy;

USE NCBI_Taxonomy;

CREATE TABLE `ncbi_names` (
  `tax_id` mediumint(11) unsigned NOT NULL default '0',
  `name_txt` varchar(255) NOT NULL default '',
  `unique_name` varchar(255) default NULL,
  `name_class` varchar(32) NOT NULL default '',
  KEY `tax_id` (`tax_id`),
  KEY `name_class` (`name_class`),
  KEY `name_txt` (`name_txt`)
);

CREATE TABLE `ncbi_nodes` (
  `tax_id` mediumint(11) unsigned NOT NULL default '0',
  `parent_tax_id` mediumint(8) unsigned NOT NULL default '0',
  `rank` varchar(32) default NULL,
  `embl_code` varchar(16) default NULL,
  `division_id` smallint(6) NOT NULL default '0',
  `inherited_div_flag` tinyint(4) NOT NULL default '0',
  `genetic_code_id` smallint(6) NOT NULL default '0',
  `inherited_GC_flag` tinyint(4) NOT NULL default '0',
  `mitochondrial_genetic_code_id` smallint(4) NOT NULL default '0',
  `inherited_MGC_flag` tinyint(4) NOT NULL default '0',
  `GenBank_hidden_flag` smallint(4) NOT NULL default '0',
  `hidden_subtree_root_flag` tinyint(4) NOT NULL default '0',
  `comments` varchar(255) default NULL,
  PRIMARY KEY  (`tax_id`),
  KEY `parent_tax_id` (`parent_tax_id`)
);
```

## Populating the database

We now load the node and name dumps into our database. The trick is to do this from the command line, and make sure MySQL can access the directory the *.dmp files are in.

```
LOAD DATA INFILE '/Users/rob/names.dmp' 
INTO TABLE ncbi_names 
FIELDS TERMINATED BY '\t|\t' 
LINES TERMINATED BY '\t|\n' 
(tax_id, name_txt, unique_name, name_class);
```

Allow a little time for the data to upload.

Next load the nodes:

```
LOAD LOCAL DATA INFILE 'nodes.dmp' 
INTO TABLE ncbi_nodes 
FIELDS TERMINATED BY '\t|\t' 
LINES TERMINATED BY '\t|\n' 
(tax_id, parent_tax_id,rank,embl_code,division_id,
 inherited_div_flag,genetic_code_id,inherited_GC_flag,
 mitochondrial_genetic_code_id,inherited_MGC_flag,
 GenBank_hidden_flag,hidden_subtree_root_flag,comments);
```
