#!/bin/bash

. ./${1}


num_tasks=0
for (( machine=0; machine<$num_machines; machine++ ))
do
    echo "lightfield ${machine_index[$machine]}: load = ${load[$machine]}"
    num_tasks=$[${num_tasks}+${load[$machine]}]
done
echo "num_tasks=${num_tasks}"

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
    echo "task ${task}: machine lightfield${machine_index[load_index]}"
done

# the name of the folder
foldername="LightField_"`eval date +%y%m%d_%H%M`
echo "The result folder name: ${foldername}"
# write it to the configuration file
echo "foldername=${foldername}" >> "${1}"


for (( machine=0; machine<$num_machines; machine++ ))
do
    if [ ${load[$machine]} -ne 0 ]
    then
        
        if [ ${machine_index[$machine]} -eq 102 ] 
        then
        
            ssh "lightfield${machine_index[$machine]}" "test -d TO_SORT/${foldername} || mkdir -p TO_SORT/${foldername}"
            ssh "lightfield${machine_index[$machine]}" "cp ~/TO_SORT/sfft_temp/* ~/TO_SORT/${foldername}/"
            scp ${config_file} lightfield${machine_index[$machine]}:~/TO_SORT/${foldername}/
        else 
            ssh "lightfield${machine_index[$machine]}" "test -d ${foldername} || mkdir -p ${foldername}"
            if [ ${machine_index[$machine]} -eq 21 ] || [ ${machine_index[$machine]} -eq 105 ]
            then
                ssh "lightfield${machine_index[$machine]}" "cp ~/sfft_temp/* ~/${foldername}/"
            else
                if [ ${machine_index[$machine]} -ge 106 ] && [ ${machine_index[$machine]} -le 115 ]
                then
                    ssh "lightfield${machine_index[$machine]}" "cp ~/sfft_temp/* ~/${foldername}/"
                else
                    scp build/nonint_sfft lightfield${machine_index[$machine]}:~/${foldername}/nonint_sfft
                    scp /usr/lib/libfftw3.so.3 lightfield${machine_index[$machine]}:~/${foldername}/
                    scp /usr/lib/libconfig++.so.8 lightfield${machine_index[$machine]}:~/${foldername}/
                fi
            fi
            scp ${config_file} lightfield${machine_index[$machine]}:~/${foldername}/
        fi
    fi
done

# actually running the stuff
for (( task=0; task<$num_tasks; task++ ))
do
    printf -v x_file "z_%04d.dat" ${task}
    printf -v out_file "log_%04d.out" ${task}
    printf -v err_file "log_%04d.err" ${task}


    echo "task ${task}: Run on machine lightfield${machine_index[${machines[${task}]}]} in folder ${foldername}"

    if [ ${machine_index[${machines[$task]}]} -ge 106 ] && [ ${machine_index[${machines[$task]}]} -le 115 ]
    then
        scp data/${x_file} lightfield${machine_index[${machines[${task}]}]}:~/${foldername}/${x_file}
        gnome-terminal -x ssh lightfield${machine_index[${machines[${task}]}]} \
            "cd ${foldername}; "\
            "chmod +x nonint_sfft; "\
            "export LD_LIBRARY_PATH=~/${foldername}; "\
            "./nonint_sfft ${config_file} ${task} ${num_tasks} > ${out_file} 2> ${err_file} < /dev/null"
    else
        if [ ${machine_index[${machines[$task]}]} -eq 102 ] 
        then
            scp data/${x_file} lightfield${machine_index[${machines[${task}]}]}:~/TO_SORT/${foldername}/${x_file}
            ssh lightfield${machine_index[${machines[${task}]}]} \
                "cd TO_SORT/${foldername}; "\
                "chmod +x nonint_sfft; "\
                "export LD_LIBRARY_PATH=~/TO_SORT/${foldername}; "\
                "nohup ./nonint_sfft ${config_file} ${task} ${num_tasks} > ${out_file} 2> ${err_file} < /dev/null &"
        else
            scp data/${x_file} lightfield${machine_index[${machines[${task}]}]}:~/${foldername}/${x_file}
            ssh lightfield${machine_index[${machines[${task}]}]} \
                "cd ${foldername}; "\
                "chmod +x nonint_sfft; "\
                "export LD_LIBRARY_PATH=~/${foldername}; "\
                "nohup ./nonint_sfft ${config_file} ${task} ${num_tasks} > ${out_file} 2> ${err_file} < /dev/null &"
        fi
    fi

done

