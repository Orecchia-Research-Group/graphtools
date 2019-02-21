#! /bin/bash
#USAGE pbsrun.sh name num_nodes type(cpu2333 || cpu3000) memory(#g) time(:::) command

sed -e 's+%1+'$1'+' -e 's+%2+'$2'+' -e 's+%3+'$3'+' -e 's+%4+'$4'+' -e 's+%5+'"$5"'+' -e 's+%6+'"$6"'+' mytemplate.pbs
