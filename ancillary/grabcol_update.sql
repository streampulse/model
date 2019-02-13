use sp;
-- describe grabcols;
alter table grabcols add column method varchar(40);
alter table grabcols add column write_in varchar(40);
alter table grabcols add column addtl varchar(40);

alter table grabdata add column method varchar(40) after value;
alter table grabdata add column write_in varchar(40) after method;
alter table grabdata add column addtl varchar(40) after write_in;

-- alter table grabdata drop column addtl;
