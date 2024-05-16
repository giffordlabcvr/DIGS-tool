Screening database dump files can be generated using the 'mysqldump program (part of MySQL distribution). For example to backup a single database:

```
mysqldump -u root -p[root_password] [database_name] > dumpfilename.sql
```

Databases can be restored as follows:

```
mysql -u root -p[root_password] [database_name] < dumpfilename.sql
```

**Note**: you must create the database before you perform the above command.