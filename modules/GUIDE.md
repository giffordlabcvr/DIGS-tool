**Organisation of the DIGS code: PERL module library**
------------------------------------------------------------------------------------

### **Modules in the DIGS directory**

#### Main modules

```
- DIGS.pm             # Database-integrated genome screening - main functions
- TargetDB.pm         # Managing the screening directory
- Utility.pm          # Utility functions
- Nomenclature.pm     # Locus ID creation and management 
- Test.pm             # Tests
```

#### Modules for setting up

```
- Initialise.pm       # General set-up (loading the database etc)
- ScreenBuilder.pm    # Setting up a screen
```

#### Modules used for running screens and merging loci

```
- Classify.pm         # Classify sequences using BLAST
- CrossMatch.pm       # Capture information about cross-matching during DIGS
- Defragment.pm       # Functions for clustering, defragmenting, consolidating loci
- Extract.pm          # Functions for extracting sequences from FASTA files 
```

### **Modules in the Interface directory**

```
- BLAST.pm            # A Perl interface to the BLAST executables
- MySQLtable.pm       # A Perl interface to a MySQL table
```

### **Modules in the Base directory**

```
- Console.pm          # Basic console functions
- FileIO.pm           # Basic file IO etc 
- DevTools.pm         # Tools used for debugging
```
