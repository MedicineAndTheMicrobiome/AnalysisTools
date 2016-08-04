#!/bin/csh

set date_tag = `date +"%Y%m%d"`
set tar_name = SIMPER_${date_tag}.tgz

cd ..

echo Building list of .svn directories/files to include in tar





find Simper |egrep 'build_SIMPER_tgz.csh|tgz' >/tmp/exclude_files
echo $tar_name >> /tmp/exclude_files

cd Simper

tar -czvf $tar_name -X /tmp/exclude_files ../Simper

echo done.