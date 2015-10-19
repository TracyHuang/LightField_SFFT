#!/bin/bash

. ./${1}

#num_tasks=0
#for (( machine=0; machine<$num_machines; machine++ ))
#do
    #echo "lightfield ${machine}: load = ${load[$machine]}"
    #num_tasks=$[${num_tasks}+${load[$machine]}]
#done
#echo "num_tasks=${num_tasks}"

# specifying the machine
#load_index=0
#cur_load_task=0
#for (( task=0; task<$num_tasks; task++ ))
#do
    #if [ $cur_load_task -lt ${load[$load_index]} ]
    #then
        #cur_load_task=$[$cur_load_task+1]
    #else
        #load_index=$[$load_index+1]
        #while [ ${load[$load_index]} -eq 0 ]
        #do 
            #load_index=$[$load_index+1]
        #done
        #cur_load_task=1
    #fi

    #machines[$task]=$load_index
    #echo "task ${task}: machine lightfield${load_index}"
#done

# the name of the folder
# foldername="LightField_${1}"

# collecting the stuff
# only get the date
date_e=`echo ${foldername} | rev | cut -d '_' -f '1 2' | rev`
date_p=`echo ${date_e} | cut -d '.' -f 1`
result_folder=result_${date_p}/result
output_folder=${result_folder}/output
test -d ${result_folder} || mkdir -p ${result_folder}
test -d ${output_folder} || mkdir -p ${output_folder}

echo "Result is in ${result_folder}"

y_file="y_sfft.dat"
debug_file="result.debug"
x_file="x.dat"
log_out_file="log.out"
log_err_file="log.err"
#for (( task=0; task<$num_tasks; task++ ))
for (( machine=0; machine<$num_machines; machine++ ))
do
    if [ ${load[$machine]} -ne 0 ]
    then 
        #printf -v result_y "y_sfft_%04d.dat" ${task}
        #printf -v result_debug "result_%04d.debug" ${task}
        #printf -v data_x "x_%04d.dat" ${task}
        #printf -v log_out "log_%04d.out" ${task}
        #printf -v log_err "log_%04d.err" ${task}

        if [ ${machine_index[$machine]} -eq 102 ] 
        then
            existed=$(ssh lightfield${machine_index[${machine}]} "test -e TO_SORT/${foldername}" && echo 1 || echo 0)

            if [ $existed == "1" ]
            then
                scp lightfield${machine_index[$machine]}:~/TO_SORT/${foldername}/*.dat ${result_folder}/
                scp lightfield${machine_index[$machine]}:~/TO_SORT/${foldername}/*.debug ${result_folder}/
                scp lightfield${machine_index[$machine]}:~/TO_SORT/${foldername}/*.out ${output_folder}
                scp lightfield${machine_index[$machine]}:~/TO_SORT/${foldername}/*.err ${output_folder}
            else
                echo "task ${task}: machine lightfield${machine_index[${machine}]}, cannot find it"
            fi

        else
            existed=$(ssh lightfield${machine_index[${machine}]} "test -e ${foldername}" && echo 1 || echo 0)

            if [ $existed == "1" ]
            then
                scp lightfield${machine_index[$machine]}:~/${foldername}/*.dat ${result_folder}/
                scp lightfield${machine_index[$machine]}:~/${foldername}/*.debug ${result_folder}/
                scp lightfield${machine_index[$machine]}:~/${foldername}/*.out ${output_folder}
                scp lightfield${machine_index[$machine]}:~/${foldername}/*.err ${output_folder}
                #scp lightfield${machine}:~/${foldername}/${log_err_file} ${output_folder}/${log_err}
            else
                echo "task ${task}: machine lightfield${machine_index[${machine}]}, cannot find it"
            fi

        fi
        
    fi
done

