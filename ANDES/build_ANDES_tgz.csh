#!/bin/csh

echo Run me from within ANDES directory
echo Results will be in up one directory.

set date_tag = `date +"%Y%m%d"`
set tar_name = `echo "ANDES_"$date_tag".tgz"`

cd ..
find ANDES | grep "\.svn" > _tar_exclude 

tar -czvf $tar_name -X _tar_exclude ANDES
