#!/bin/bash

# Weighted Diameter
#diameter=(.35 .4 .45 .5 .55 .6 .65 .7 .75 .80 .85 .90)
diameter=(.725 .775 .825 .875 .925 .95)
#diameter=(.975 1)

#qsub -P salabs wd_fuzzyGraph_tasks_prelude.qsub -N wd_fuzzyGraph_tasks_prelude

for h in ${diameter[@]}; do
	#qsub -P salabs -hold_jid wd_fuzzyGraph_tasks_prelude -N wd_fuzzyGraph_task_${h}  wd_fuzzyGraph_tasks.qsub $h
	qsub -P salabs -hold_jid wd_fuzzyGraph_task_${h} -N wd_fuzzyGraph_val_${h} wd_fuzzyGraph_validation.qsub $h
done

