**Initial screening run**

Before running a screen for the first time, create the screening database (a one-off procedure), as follows:

```
./pipeline –m=1 –i=[path to control file]
```
Once the screening database is successfully created, a screen can be executed as follows:

```
./pipeline –m=2 –i=[path to control file]
```
Progress is written to the terminal, and can also be monitored by issuing SQL queries against the relevant screening database. A screen can be stopped at any time. The next time the tool is restarted, it will initiate screening at the point it left off.

**Running a screen in the background**

Since similarity search-based screens can take hours or days to complete, you may want to run your screen 'in the background', particularly if you are running DIGS on a networked server that you connect to remotely. This can be done as follows. First start the pipeline as follows:

```
./pipeline –m=2 –i=[path to control file] > [path to a log file, e.g. screen.log]
```

Once this command has been executed, press Ctrl+Z and then type:

```
bg
```
You should then see this output:

```
./pipeline.pl -i ctl/example_1_EVE.ctl -m=2 > screen.log &
```

To check progress on the backgrounded screen, view the log file. For example using:
```
tail screen.log 
```


