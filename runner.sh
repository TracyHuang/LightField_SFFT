#!/bin/bash

. ./${1}

#get start time
res1=$(date +%s.%N)

# first check whether all need programs and config files exist
if [ ! -f $input_program_path ];
then
	echo "1" 
	echo " $file_program_path does not exist"
	exit
fi

if [ ! -f $config_file_path ];
then 
	echo "2"
	echo " $config_file_path does not exist"
	exit
fi

if [ ! -f $main_program_path ];
then 
	echo "3"
	echo " $main_progeam_path does not exist"
	exit
fi

if [ ! -f $output_file_path ];
then 
	echo "4"
	echo " $output_file_path does not exist"
	exit
fi


# first we ask the user to decide the starting point of the program
echo "Please select the starting point of the script."
echo "1--->start at the input generation."
echo "2--->start at the main program."
echo "3--->start at the checking loop."
echo "4--->start at the result collection."
echo "Please check documentation if you are not sure where to start"
read START_POINT

#echo "Do you want to run program with full RGB or gray scale? (1 for RGB/ 0 for gray scale"
#read FULL_DATA

#if [ $FULL_DATA -eq 1 ] 
#then 
#max_channel=2
#else
#max_channel=0
#fi

if [ $START_POINT -eq 1 ] && [ $n3 -eq $n4 ] && [ $n2 -eq $n1 ]
then 
echo "Do you want to use copying to speed up calculation? (1 for YES/0 for NO)"
read if_copy
else
if_copy=0
fi



for channel in {0..2}   
do

		

	# first we will generate input and copy them to data_file_folder
	if [ $START_POINT -eq 1 ]
	then
		echo "start input generation---"
		echo "$input_program_path $input_images_path $n1 $n2 $n3 $n4 $num_splits $if_copy $channel"
		$input_program_path $input_images_path $n1 $n2 $n3 $n4 $num_splits $if_copy $channel
		echo "---------input generation finished"
		rm -r $data_file_folder
		mkdir $data_file_folder
		mv ./*.dat $data_file_folder
	fi

	# before run the main program 
	# we add Control master to enable us using same ssh channel for several scp/ssh
	# echo "Host *" >> ~/.ssh/config
	# echo "ControlMaster auto" >> ~/.ssh/config
	# echo "ControlPath ~/.ssh/sockets/ssh-socket-%r-%h-%p" >> ~/.ssh/config




	# run the main program
	if [ $START_POINT -le 2 ]
	then
		num_tasks=0
		for (( machine=0; machine<$num_machines; machine++ ))
		do
			echo "${machine_index[$machine]}: load = ${load[$machine]}"
			num_tasks=$[${num_tasks} + ${load[$machine]}]
		done
		echo "num_tasks=${num_tasks}"

		# specifying which machine gets which tasks
		load_index=0
		cur_load_task=0
		for (( task=0; task<$num_tasks; task++))
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
			echo "task ${task}: machine ${machine_index[load_index]}"
		done


		foldername="LightField_"`eval date +%y%m%d_%H%M`
		echo "The result folder name: ${foldername}"
		#write it to the configuration file
		echo "foldername=${foldername}" >> "${1}"

		task=0
		for (( machine=0; machine<$num_machines; machine++ ))
		do
			if [ ${load[$machine]} -ne 0 ]
			then
				ssh "${user_name[$machine]}@${machine_index[$machine]}" "test -d ${foldername} || mkdir -p ${foldername}"
				scp $main_program_path "${user_name[$machine]}@${machine_index[$machine]}":~/${foldername}/nonint_sfft
				scp ./lib/libconfig++.so.9 "${user_name[$machine]}@${machine_index[$machine]}":~/${foldername}/
				scp ./lib/libfftw3.so.3 "${user_name[$machine]}@${machine_index[$machine]}":~/${foldername}/
				scp ./lib/libconfig++.so.8 "${user_name[$machine]}@${machine_index[$machine]}":~/${foldername}/
				scp ${config_file_path} "${user_name[$machine]}@${machine_index[$machine]}":~/${foldername}		

				index=0
				for (( ; index<${load[$machine]}; index++ ))
				do	
					printf -v x_file "x_%04d.dat" ${task}
					printf -v out_file "log_%o4d.out" ${task}
					printf -v err_file "log_%o4d.err" ${task}

					echo "task ${task}: Run on machine ${user_name[$machine]}@${machine_index[${machines[${task}]}]} in folder ${foldername}"
						
					scp $data_file_folder/${x_file} "${user_name[$machine]}@${machine_index[$machine]}":~/${foldername}/${x_file}
					ssh "${user_name[$machine]}@${machine_index[$machine]}" \
						"cd ${foldername}; " \
						"chmod +x nonint_sfft; " \
						"export LD_LIBRARY_PATH=~/${foldername}; " \
						"nohup ./nonint_sfft ${config_file} ${task} ${num_tasks} > ${out_file} 2> ${err_file} < /dev/null &"

						task=$[$task+1]
				done
			fi
		done
		echo "all tasks start running"
	fi


	# check every ten minutes whether the program finishes running or not
	if [ $START_POINT -le 3 ]
	then	
		total_num=0
		while true
		do
			for (( machine=0; machine<$num_machines; machine++ ))
			do
				if [ ${load[$machine]} -ne 0 ]
				then
					num=$(ssh ${user_name[$machine]}@${machine_index[$machine]} "ps auwxr"  | grep nonint_sfft | wc -l)
					total_num=$[$total_num+$num]
					echo "machine ${machine_index[$machine]}: load = ${load[$machine]}, process = ${num}"
				fi
			done


			for (( task=0; task<$num_tasks; task++ ))
			do

				printf -v x_file "~/${foldername}/x_%04d.dat" ${task}
				printf -v y_file "~/${foldername}/y_sfft_%04d.dat" ${task}
				#echo $task
				#echo ${machines[${task}]}
				#echo ${user_name[${machines[${task}]}]}
				existed=$(ssh ${user_name[${machines[${task}]}]}@${machine_index[${machines[$task]}]} "test -e ${x_file}" && echo 1 || echo 0)

				if [ $existed == "1" ]
				then

					x_file_size=$(ssh ${user_name[${machines[${task}]}]}@${machine_index[${machines[${task}]}]} "stat -c%s ${x_file}")
					y_file_size=$(ssh ${user_name[${machines[${task}]}]}@${machine_index[${machines[${task}]}]} "stat -c%s ${y_file}")

					#if [ $x_file_size -eq "" ] || [ $y_file_size -eq "" ]; then
						#echo "error"
					#fi

					ratio=$[$[${y_file_size}*100]/${x_file_size}]
					echo "task ${task}: machine ${machine_index[${machines[$task]}]}, complete ${y_file_size}/${x_file_size} = ${ratio}%"
				else
					echo "task ${task}: machine ${machine_index[${machines[$task]}]}, cannot find it"
				fi
			done

			if [ $total_num -eq 0 ] 
			then
				echo "Main program finished execution"
				break
			else 
				echo "There are still $total_num of processes running"
				total_num=0		
				sleep 5m
			fi

		done
	fi

	# collect result and do ifft to get reconstructed photos

	if [ $START_POINT -le 4 ]
	then
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

		for (( machine=0; machine<$num_machines; machine++ ))
		do
			if [ ${load[$machine]} -ne 0 ]
			then
				existed=$(ssh ${user_name[$machine]}@${machine_index[${machine}]} "test -e ${foldername}" && echo 1 || echo 0)
				if [ $existed == "1" ]
				then
					scp ${user_name[$machine]}@${machine_index[$machine]}:~/${foldername}/*.dat ${result_folder}/
					scp ${user_name[$machine]}@${machine_index[$machine]}:~/${foldername}/*.debug ${result_folder}/
					scp ${user_name[$machine]}@${machine_index[$machine]}:~/${foldername}/*.out ${output_folder}
					scp ${user_name[$machine]}@${machine_index[$machine]}:~/${foldername}/*.err ${output_folder}
					#scp lightfield${machine}:~/${foldername}/${log_err_file} ${output_folder}/${log_err}
				else
					echo "task ${task}: machine ${machine_index[${machine}]}, cannot find it"
				fi
			fi

			
		done

		# run output handle program
		echo "start to generate image--"
		test -d data_result && rm -rf data_result
		mkdir -p data_result
		cp ${result_folder}/y_sfft_*.dat ./data_result
		$output_file_path ./data_result $n1 $n2 $n3 $n4 $channel

		if [ $channel -eq 0 ]
		then 
		cur_folder=y_result
		elif [ $channel -eq 1 ]
		then
		cur_folder=u_result
		elif [ $channel -eq 2 ]
		then
		cur_folder=v_result
		fi

		test -d $cur_folder && rm -rf $cur_folder
		mkdir -p $cur_folder
		mv ./*.dat $cur_folder
	fi
done

$merge_file_path ./y_result/ ./u_result/ ./v_result/ $n1 $n2 $n3 $n4 

echo "finished"

# get finished time and print out execution time information
res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)

printf "Total runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds


 
