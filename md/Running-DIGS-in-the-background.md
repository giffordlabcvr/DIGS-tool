**Running DIGS in the Background**

Genome screening can take hours, days or even weeks to complete. It is often convenient to run large-scale screens on a remote server, that has a lot of compute power for running BLAST, and plenty of disk space for storing target genome sequences. So that we don't have to depend on a live connection to the server to keep the screen running, DIGS can be made to run 'in the background'. Below is a very brief guide to how to do this for those who - like me - are apt to forget these things.

Normal procedure is as follows:

* First, initiate the screen in the normal way to check that it gets underway as expected:

```
./pipeline.pl -i ../projects/eve/eve_vertebrates.ctl -m=2

	######################################################################
	#                                                                    #
	#                              DIGS 1.1                              #
	#                Database-Integrated Genome Screening                #
	#                         Robert J. Gifford                          #
	#                   <robert.gifford@glasgow.ac.uk>                   #
	#                                                                    #
	######################################################################

	  Connecting to DB:  EVE_vertebrates
	  Reference library: 39411 amino acid sequences
	  Probes:            37 amino acid sequences
	  Targets:           127 target files
	  Searches to run    4427

	  Starting database-integrated genome screening

          Screen: 1 (%0.02 done): 'Homo_sapiens' file 'chrUn.fa' with probe PCV_Cap
		 # 0 newly identified hits
```

* That all seems fine, so press Control-Z to stop the screen, and execute the same command, but this time send the screen output to a file as follows:

```
./pipeline.pl -i ../projects/eve/eve_vertebrates.ctl -m=2 > vertebrate_screen.log

```

* Now press Control-Z to stop the screen again, and execute the background (bg) command.

```
bg
```

You will see the backgrounded command printed to the screen, looking something like this:

```
[5]+ ./pipeline.pl -i ../projects/eve/eve_vertebrates.ctl -m=2 > vertebrate_screen.log &
```

* To check if the screen is running as expected, inspect the log file, e.g using 'tail':

```
tail vertebrate_screen.log 
```

* If the log file shows the output of DIGS screening as expected, you can log out of the server and the screen will continue running.

Note: when you do log out, you will get a warning about stopped jobs, which you can ignore.

