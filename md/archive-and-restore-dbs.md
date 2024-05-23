# Archiving and Restoring Databases

Database dump files can be generated using the `mysqldump` program (part of the MySQL distribution). For example:

```bash
mysqldump -u root -p [database_name] > dumpfilename.sql
```

To restore a database, use the following command:

```
mysql -u root -p [database_name] < dumpfilename.sql
```

**Note**:  Make sure you create the database before running the restore command. You can create a database with:

```
mysql -u root -p -e "CREATE DATABASE [database_name];"
```

Replace [database_name] with the name of your database and dumpfilename.sql with the name of your dump file.
