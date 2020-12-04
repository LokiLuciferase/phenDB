
drop database SED_REPLACE_ME;
drop user 'SED_REPLACE_ME'@'localhost';
create database SED_REPLACE_ME;
create user 'SED_REPLACE_ME'@'localhost' identified by 'SED_REPLACE_ME';
grant all privileges on SED_REPLACE_ME.* to 'SED_REPLACE_ME'@'localhost';
