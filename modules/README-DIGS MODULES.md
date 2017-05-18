**Organisation of the DIGS code**
------------------------------------------------------------------------------------

### **Modules in the DIGS directory**

### **Main modules**

```
- TargetDB.pm         # Main module for managing the screening directory
- DIGS.pm	            # Main module for running screening-related (-m) functions 2-5 
- Utility.pm          # Main module for utility functions
- Nomenclature.pm     # Main module for locus nomenclature functions 
- Test.pm            
```

### **Modules for setting up**

```
- Initialise.pm       # General set-up (loading the database etc)
- ScreenBuilder.pm    # Setting up a screen
```

### **Modules used for running screens and merging loci**

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
