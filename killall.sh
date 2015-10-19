#!/bin/bash

. ./${1}

echo $num_machines

for (( machine=0; machine<$num_machines; machine++ ))
do
    if [ ${load[$machine]} -ne 0 ]
    then
		echo "ssh ${user_name[$machine]}@${machine_index[$machine]}" 
        ssh "${user_name[$machine]}@${machine_index[$machine]}" "killall -9 nonint_sfft"
    fi
done


