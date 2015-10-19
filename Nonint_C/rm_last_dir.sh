#!/bin/bash

. ./${1}

for (( machine=0; machine<$num_machines; machine++ ))
do
    if [ ${load[$machine]} -ne 0 ]
    then
        echo "lightfield${machine_index[$machine]}: rm -rI ${foldername}"
        ssh "lightfield${machine_index[$machine]}" "rm -rI ${foldername}"
    fi
done

