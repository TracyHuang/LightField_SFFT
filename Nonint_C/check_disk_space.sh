#!/bin/bash

. ./${1}

for (( machine=0; machine<$num_machines; machine++ ))
do
    echo "---------------------lightfield${machine} usage:---------------------------"
    ssh "lightfield${machine}" "df -h | grep /dev/sda2"
done



