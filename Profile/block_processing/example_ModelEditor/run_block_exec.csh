#!/bin/csh

set model_dir=example.dir
set st=summary_table.tsv
set factors=factor_file.tsv
set output_dir=outdir


set covariates=$model_dir/covariates

set groups=`ls $model_dir/groups`

foreach grp ($groups)

	echo "Summary Table: " $st
	echo " Factors File: " $factors
	echo "   Covariates: " $covariates
	echo "        Group: " $model_dir/groups/$grp
	echo "   Output Dir: " $output_dir/$grp

	if (0) then
		../block_execute_models.pl \
			-s $st \
			-f $factors \
			-c $covariates  \
			-g $model_dir/groups/$grp \
			-o $output_dir/$grp
	endif 

	echo "Done..."
	echo

end


