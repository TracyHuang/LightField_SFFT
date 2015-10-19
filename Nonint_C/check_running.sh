#!/bin/bash

. ./${1}

num_tasks=0
for (( machine=0; machine<$num_machines; machine++ ))
do
    num_tasks=$[${num_tasks}+${load[$machine]}]
done

# specifying the machine
load_index=0
cur_load_task=0
for (( task=0; task<$num_tasks; task++ ))
do
    if [ $cur_load_task -lt ${load[$load_index]} ]
    then
        cur_load_task=$[$cur_load_task+1]
    else
        load_index=$[$load_index+1]
        while [ ${load[$load_index]} -eq 0 ]
        do 
            load_index=$[$load_index+1]
        done
        cur_load_task=1
    fi

    machines[$task]=$load_index
done

for (( machine=0; machine<$num_machines; machine++ ))
do
    if [ ${load[$machine]} -ne 0 ]
    then
        num=$(ssh lightfield${machine_index[$machine]} "ps auwxr"  | grep nonint_sfft | wc -l)
        echo "lightfield ${machine_index[$machine]}: load = ${load[$machine]}, process = ${num}"
    fi
done


for (( task=0; task<$num_tasks; task++ ))
do

    printf -v x_file "${foldername}/z_%04d.dat" ${task}
    printf -v y_file "${foldername}/y_sfft_%04d.dat" ${task}

    if [ ${machine_index[${machines[$task]}]} -eq 102 ] 
    then
        existed=$(ssh lightfield${machine_index[${machines[$task]}]} "test -e TO_SORT/${x_file}" && echo 1 || echo 0)
        if [ $existed == "1" ]
        then

            x_file_size=$(ssh lightfield${machine_index[${machines[${task}]}]} "stat -c%s TO_SORT/${x_file}")
            y_file_size=$(ssh lightfield${machine_index[${machines[${task}]}]} "stat -c%s TO_SORT/${y_file}")

            ratio=$[$[${y_file_size}*100]/${x_file_size}]
            echo "task ${task}: machine lightfield${machine_index[${machines[$task]}]}, complete ${y_file_size}/${x_file_size} = ${ratio}%"
        else
            echo "task ${task}: machine lightfield${machine_index[${machines[$task]}]}, cannot find it"
        fi

    else
        existed=$(ssh lightfield${machine_index[${machines[$task]}]} "test -e ${x_file}" && echo 1 || echo 0)

        if [ $existed == "1" ]
        then

            x_file_size=$(ssh lightfield${machine_index[${machines[${task}]}]} "stat -c%s ${x_file}")
            y_file_size=$(ssh lightfield${machine_index[${machines[${task}]}]} "stat -c%s ${y_file}")

            #if [ $x_file_size -eq "" ] || [ $y_file_size -eq "" ]; then
                #echo "error"
            #fi

            ratio=$[$[${y_file_size}*100]/${x_file_size}]
            echo "task ${task}: machine lightfield${machine_index[${machines[$task]}]}, complete ${y_file_size}/${x_file_size} = ${ratio}%"
        else
            echo "task ${task}: machine lightfield${machine_index[${machines[$task]}]}, cannot find it"
        fi
    fi
done



