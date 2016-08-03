#!/bin/csh

set PROJECT_DIR = /local/devel/DAS/users/kli/SVN/DAS/16sDataAnalysis/trunk/site_comparison/hierarchical_processing/example

cd $PROJECT_DIR
echo "Create directories in " $PROJECT_DIR "? (Y)"
set result = $<

if ($result == "y" || $result == "Y") then
	echo "Ok... commencing.\n"	
else
	echo "Aborting."
	exit
endif

if (-e $PROJECT_DIR) then
	echo $PROJECT_DIR " found."
else
	echo $PROJECT_DIR " not found.  Making..."
	mkdir $PROJECT_DIR
endif

set root_dir = (\
	Sequencing \
	Annotation \
	Mothur \
	MetaData \
	Sequence_Based_Analyses \
	Count_Based_Analyses \
	Count_Based_Analyses/Summary_Tables \
	Count_Based_Analyses/Summary_Tables/All_Samples \
	Count_Based_Analyses/Analyses \
	Count_Based_Analyses/Analyses/Permanova \
	Count_Based_Analyses/Analyses/MALR \
	Count_Based_Analyses/Analyses/DiversityRegression \
	Count_Based_Analyses/Analyses/Corbata \
	Count_Based_Analyses/Analyses \
)

foreach dir ($root_dir)

	echo $PROJECT_DIR/$dir
	mkdir $PROJECT_DIR/$dir

end
