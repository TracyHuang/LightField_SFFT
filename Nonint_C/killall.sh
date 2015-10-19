#!/bin/bash

. ./${1}

for (( machine=0; machine<$num_machines; machine++ ))
do
    if [ ${load[$machine]} -ne 0 ]
    then
        ssh "lightfield${machine_index[$machine]}" "killall -9 nonint_sfft"
    fi
done

