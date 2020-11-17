
drop database if exists phenDB_devel_LL;
drop user if exists 'phenDB_devel_LL'@'localhost';
create database phenDB_devel_LL;
create user 'phenDB_devel_LL'@'localhost' identified by 'phenDB_devel_LL';
grant all privileges on phenDB_devel_LL.* to 'phenDB_devel_LL'@'localhost';
