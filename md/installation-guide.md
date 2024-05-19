# DIGS Tool Installation and Usage Guide

## Background

The DIGS tool, developed by Rob Gifford and colleagues at the MRC-University of Glasgow Centre For Virus Research, is a screening tool that combines local BLAST search tools, Perl functions, and the MySQL database management system. It enables users to set up highly customized screens with results queryable in multiple ways.

Developers have published an installation guide for local drive use. This guide addresses potential snags and provides detailed instructions for setting up the directory structure, path and environment variables, and downloading genomic data. It covers both local installations and setups on a server.

## How DIGS Tool Works

As an example, consider screening human chromosome 2 for ERVs. The DIGS program requires information on target, probe, and reference sequences, MySQL server credentials, the database name for results, and blast search parameters in a control file. The command to run DIGS uses this control file and goes through formatting the target, screening, and storing results in MySQL tables.

## Installation on a Mac Operating System

### Preparation:

1. **Folders:** Open the terminal and set up DIGS and Genomes directories:
    ```sh
    cd ~
    mkdir DIGS
    mkdir Genomes
    ```

2. **Developer Tools:** Check if Mac developer tools are installed:
    ```sh
    xcode-select -p
    ```
    If not installed, run:
    ```sh
    xcode-select --install
    ```

3. **Shell Environment:** Create or edit the .profile file:
    ```sh
    nano .profile
    ```
    Add the following lines:
    ```sh
    export PATH="/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/mysql/bin:$PATH"
    export PATH="${PATH}:~/DIGS/"
    export DIGS_HOME=~/DIGS/
    export PATH="${PATH}:~/Genomes/"
    export DIGS_GENOMES=~/Genomes/
    ```
    Save and run:
    ```sh
    source .profile
    ```

4. **Homebrew:** Install Homebrew:
    ```sh
    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    brew doctor
    ```

### Installation of BLAST and MySQL

Use Homebrew to install BLAST and MySQL:

```sh
brew install blast
brew install mysql
```

Open .profile and remove the '#' from the MySQL path line. Then run:

`source .profile`

Start the MySQL server:

`mysql.server start`


Create a test database and exit.

### Installation of DBI and DBD::mysql

Install DBI and DBD::mysql

`sudo -H cpan DBI::DBD DBD::mysql`

Navigate to the DBD-mysql directory and run:


```perl Makefile.PL
make
make test
```

Add MySQL user and password variables to .profile and source it.


### Downloading and Testing DIGS Tool

1. Clone or download the DIGS-tool repository and move it to the DIGS directory.
2. Run the DIGS tool help command to verify installation:



## Example control file

```
***********************************************************************
Begin SCREENDB;
        db_name=digs_chr1LTR_db;
        mysql_server=localhost;
        mysql_username=digs;		
        mysql_password=digs-tool;
ENDBLOCK;

BEGIN SCREENSETS;
        output_path=./tmp/;
        query_na_fasta=./libraries/LTR_probes.faa;
        reference_na_fasta=./libraries/LTR_references.faa;
        seq_length_minimum=60;
        bitscore_min_blastn=40;
        defragment_range=100;
ENDBLOCK;

BEGIN TARGETS;
        Mammalia/H_s/complete/goldenpath_hg38/chr1.fa
ENDBLOCK;
***********************************************************************
```

## Resources

The DIGS-tool developers’ site on installation and use:

https://giffordlabcvr.github.io/DIGS-tool/

UNIX basic command line tools:

http://www.westwind.com/reference/OS-X/commandline/navigation.html

UNIX and Perl primer for biologists – excellent guide through UNIX commands and basic perl programming:

http://korflab.ucdavis.edu/Unix_and_Perl

Practical Computing for Biologists, by Haddock and Dunn – book and accompanying online material that covers UNIX, python, mysql, and more:

http://practicalcomputing.org/downloads

Mysql commands cheatsheet:

https://gist.github.com/hofmannsven/9164408





