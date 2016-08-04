#!/bin/csh

set date_tag = `date +"%Y%m%d"`
set tar_name = Corbata_${date_tag}.tgz

cd ..

echo Building list of .svn directories/files to include in tar
find Presence_CDF|egrep 'Assign_Colors.r|Combine_CDFs.r|Compare_CDFs.r|Compute_CDFs.r|Compute_Number_of_Core.r|Plot_CDFs_OtherCore.r|Plot_CDFs.r|Compare_Microbiomes.r|Compare_wAWKS.r' |grep -v '.swp'>/tmp/corbata_files
find SummaryTable_Utilities | egrep 'Summarize_SummaryTable.r|Join_Summary_Table.pl|Generate_Bootstrap_Replicates_on_SummaryTable.r'| grep -v '.swp'>>/tmp/corbata_files
find TaxonomicVariation |egrep 'TaxonomicVariation.r'| grep -v '.swp'>>/tmp/corbata_files

cp `cat /tmp/corbata_files`  Corbata/lib

find Corbata |egrep 'build_Corbata_tgz.csh|tgz' >/tmp/exclude_files
echo $tar_name >> /tmp/exclude_files
cd Corbata

tar -czvf $tar_name -X /tmp/exclude_files ../Corbata

echo done.
