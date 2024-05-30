## Tristan's guide to installing the DIGS tool on Mac

Installing the DIGS tool on your Mac PC can be a convenient way to perform screens without the need for a network connection. Do bear in mind that you'll need to keep local copies of all the genomes you plan to screen, which may eat up a fair amount of disk space.  

[Homebrew](http://brew.sh/) is really useful for installing packages - it does it automatically so you don’t have to fiddle around. It’s got quite a lot of useful stuff for science too.

So the prerequisites are:

1. BLAST, for the searches. Install this using homebrew if you’re inclined. If not, be sure that the directory containing your BLAST executables is in your path in your ['.profile'](http://www.theunixschool.com/2011/07/what-is-profile-file.html) file (note this may be called .bash_profile, or something else, depending on which shell you are using).

2. MySQL plus a client. I use Sequel Pro (http://www.sequelpro.com/), but that’s up to you. SQL is pretty easy to learn, and the manual (attached PDF) has some basic commands. Codeacademy provides quite a useful tutorial for SQL general use - it’s free and very user friendly: https://www.codecademy.com/learn/learn-sql. This is a database management system where your data will be stored for you to retrieve and manipulate.

3. Perl (inc. the DBI and DBD database interface packages). Perl should be installed on your mac by default - check by typing perl -v in terminal.

4. Mac developer tools. These weren’t installed on my computer, but you may have them. These are needed for the 'make' commands in particular. If you don’t have these installed, type this in to terminal: 

```
xcode-select —install
```
...this will install a basic developer toolkit that we need to get the software working. Should be around 130MB and takes a few mins to download and install. Annoyingly, upon installation, your computer may set XCode as the default application for opening 'developer' files (e.g with .pl, .sh, .py .xml etc etc). You'll need to reconfigure according to your preference.

Once you’ve got these sorted, install MySQL. For the most part I followed the instructions on this webpage: http://bixsolutions.net/forum/thread-8.html

So, download MySQL community server: https://dev.mysql.com/downloads/mysql/ (the .pkg), install as you would any other application in the GUI. Open preferences, at the bottom there should be a MySQL icon. Turn on MySQL. I check the option to start it up when my mac boots but up to you.

Then, add the following line to your .profile file: 

```
export PATH="/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/mysql/bin:$PATH"
```
This sets a default path to the mysql directories so you can use mysql more easily.  **NOTE** - after making changes to your '.profile' file, you need to either restart the terminal session, or execute the file using the 'source' command:

```
source .profile
```

...for the changes to take place!

Now, to create your username and password. 

```
mysqladmin -u root password NEW_ROOT_PASSWORD. 
```

Now log in to mysql: 

```
mysql -u root -p
```

...to see if it’s worked. The writer on the website prefers to use a different account to root to admin databases, don’t know why but if you agree, you can too. Log in to mysql and type.

```
mysql> GRANT ALL PRIVILEGES ON *.* to ‘username'@'localhost' IDENTIFIED BY ‘password';
```
Now you need to install the perl module that interfaces the DIGS scripts with mysql. Now we'll use CPAN -network of perl modules and content built into perl. I've always found it a bit of a nightmare to navigate, but I'm assured by my elders and betters that it's actually very useful for those using perl on a regular basis.

So, fire up CPAN:

```
perl -MCPAN -e ‘shell’
```

. It may ask if you’re ready for manual installation. You can configure if you want but I always hit ‘no’, so it will do it automatically.


In the CPAN> prompt, type the following:

```
cpan> get DBI
cpan> get DBD::mysql
cpan> exit
```

This will retrieve the packages. Now ctl-C out of CPAN and install your packages:

```
sudo perl -MCPAN -e 'install DBI'
```

Sudo (the superuser) affords admin privileges, and you’ll need to type in the password associated with the admin user (you, I’m guessing). It should then install fine, and you’ll get an error message at the end if not.

This is where the faff begins :) We’ve had mixed results with this step, but always got it work in the end. You need to build the DBD module yourself so it configures to your database settings. cd into the cpan directories: cd ~/.cpan/build and ls to find the dbd package. Mine looks like this: DBD-mysql-4.033-9Z3F8k

Type:

```
perl Makefile.PL --testuser=‘your_mysql_username' --testpassword='YOUR_PASSWORD'
```

Now it’s time to use the developer tools you downloaded. 

```
make
make test
sudo make install
```

The test may say it fails. Ignore this for now - if it becomes an issue later we will come back to it.


################

So! You now have BLAST, MySQL and the correct PERL modules installed. Now, point your browser to this webpage and download the DIGS stuff. I prefer to keep it in my home directory under DIGS, but that’s a matter of preference.

Next step is to get your genome. Again, up to you where you put it (I also keep them in my home). Read the manual for info on how to create the correct directory structure, pretty straightforward.

You now need to set environment variables for DIGS home and the genome directories- for more info try: http://www.codecoffee.com/tipsforlinux/articles/030.html - should look something like this in your .profile, replacing my file paths with yours. You can keep DIGS and your genomes directory anywhere you like (I keep them in Home), but you need to set the environment variables to reflect this.

```
export PATH="${PATH}:/Path_to_/DIGS/"
export DIGS_HOME=/Path_to_/DIGS/          

export PATH="${PATH}:/Path_to/Genomes/"
export DIGS_GENOMES=/Path_to/Genomes/
```

Restart terminal and try to find them by typing $ and hitting tab until you get a list of available environment variables:

```
$DIGS_GENOMES 
$DIGS_HOME 
```

...should be there.

Try navigating to them by cd’ing into them from the command line - you should be able to tab complete and enter the directories.

Next step is to set up your screen. You will need a files of reference and probe sequences (in FASTA), and a control file. Read the relevant parts of this wiki for more information.

Once your files are in order, navigate to DIGS directory and check to see if digs_tool.pl (the main DIGS script) works. Executing:

```
./digs_tool.pl -h
```

Should give you something like this:

```
pc231-222:~ tristandennis$ digs_tool.pl -h

```
Now you need to create your screening database, so first run as follows:

```
./digs_tool –m=1 –i=[path to control file]
```

Check to see if it’s worked by opening your SQL client and logging in with your MySQL credentials:

Familiarise yourself with the interface, and check the top-left drop down menu for the database. If it’s worked, it should be there.

Now read the rest of the manual, set up a genomes directory, format the genome (option 7), plus the rest of your screen. Of any errors you get, most likely you’ll get an error message about the format of your genome directory, or something about the DBI. If it’s DBI related, navigate back to the mysql-dbd directory and enter this command:

```
$ sudo ln -s /usr/local/mysql/lib/libmysqlclient.18.dylib /usr/lib/libmysqlclient.18.dylib
```
...although bear in mind that your 'libmysqlclient' file may have a different version number (you can find out using tab command completion).

If it's genome related, make sure your genome FASTAs are in a BLAST friendly format. (The format genome directory option included in the DIGS tool will handily format all your genome files for you).

