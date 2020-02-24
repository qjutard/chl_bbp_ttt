
# check the git status
if [ ! -z "$(git status --porcelain)" ]; then
	read -p "The git working copy is not clean so these tests cannot be associated to a commit. Continue anyway (y/n) ?" choice
	case "$choice" in
		y|Y ) echo "continuing..." && sha="";;
		n|N ) exit 0;;
		*) echo "invalid choice" && exit 1;;
	esac
else
	sha=$(git rev-parse HEAD)
fi


# go to test environment
cd ~/Documents/cornec_chla_qc/test_env/work_dir


DMMC() {
	~/Documents/cornec_chla_qc/chl_bbp_ttt/DMMC.sh "$@"
}

clean() {
	rm test_n_errors.t
	rm ../netcdf/coriolis/6901524/profiles/DMMC/DMMC_profiles/*
	rm ../netcdf/coriolis/6901524/profiles/DMMC/DMMC_filecheck/*
}

test_finished() {
	if [ -f test_n_errors.t ]; then
		echo "success"
	else
		echo "failure"
	fi
}

test_no_errors() {
	num_errors=$(sed '4q;d' test_n_errors.t)
	if [ $num_errors == 0 ]; then
		echo "success"
	else
		echo "failure"
	fi
}

test_filechecker() {
	# Run the filechecker
	java -jar ~/Documents/filechecker/ValidateSubmit.jar -text-result coriolis ~/Documents/filechecker/spec/ ~/Documents/cornec_chla_qc/test_env/netcdf/coriolis/6901524/profiles/DMMC/DMMC_filecheck/ ~/Documents/cornec_chla_qc/test_env/netcdf/coriolis/6901524/profiles/DMMC/DMMC_profiles/
	rm ./ValidateSubmit_LOG

	# Verify each .filechecker
	all_fck_ok=true
	for filename in $(ls ~/Documents/cornec_chla_qc/test_env/netcdf/coriolis/6901524/profiles/DMMC/DMMC_filecheck/*)
	do
		line=$(sed '3q;d' $filename)
		if [ "$line" != "STATUS: FILE-ACCEPTED" ]; then
			all_fck_ok=false
		fi
	done
	
	if $all_fck_ok; then
		echo "success"
	else
		echo "failure"
	fi
}


# initialize json
report_name=../test_logs/test_report_$(date +%Y-%m-%d_%H-%M).json
report=$(jq -n \
			--arg sha "$sha" \
			'{sha:$sha, tests:[]}')


#clean
#rm DEEP_EST.t

# loop on tests
while read -r line; do
	
	# Clean outputs
	#clean
	
	# Run the test
	echo $line
	#$line
	
	# Verify success	

	finished=$(test_finished)	
	
	if [ $finished == "success" ]; then
		no_errors=$(test_no_errors)
	else
		no_errors="pending"
	fi

	if [ $no_errors == "success" ]; then
		filechecker=$(test_filechecker)
	else
		filechecker="pending"
	fi
	
		
	# Update report
	report=$(echo $report | jq \
				--arg cmd "$line" \
				--arg fin "$finished" \
				--arg err "$no_errors" \
				--arg fck "$filechecker" \
				'.tests += [{command:$cmd,finished:$fin,no_errors:$err,filechecker:$fck}]')

done < ~/Documents/cornec_chla_qc/chl_bbp_ttt/Tests/test_runs.t 


# Write global commentary
n_tests=$(echo $report | jq '.tests[].finished' | wc -l)
n_success_finished=$(echo $report | jq '.tests[].finished' | grep -o "success" | wc -l)
n_success_no_errors=$(echo $report | jq '.tests[].no_errors' | grep -o "success" | wc -l)
n_success_filechecker=$(echo $report | jq '.tests[].filechecker' | grep -o "success" | wc -l)

report=$(echo $report | jq \
			--arg fin "$n_success_finished/$n_tests successfully finished" \
			--arg err "$n_success_no_errors/$n_tests had no errors"\
			--arg fck "$n_success_filechecker/$n_tests passed the filechecker" \
			'. += {finished_description:$fin,no_errors_description:$err,filechecker_description:$fck}')


# Write global status
if [ $n_success_finished == $n_tests ]; then
	all_success_finished="success"
	if [ $n_success_no_errors == $n_tests ]; then
		all_success_no_errors="success"
		if [ $n_success_filechecker == $n_tests ]; then
			all_success_filechecker="success"
		else
			all_success_filechecker="failure"
		fi
	else
		all_success_no_errors="failure"
		all_success_filechecker="pending"
	fi
else
	all_success_finished="failure"
	all_success_no_errors="pending"
	all_success_filechecker="pending"
fi

report=$(echo $report | jq \
			--arg fin "$all_success_finished" \
			--arg err "$all_success_no_errors"\
			--arg fck "$all_success_filechecker" \
			'. += {finished:$fin,no_errors:$err,filechecker:$fck}')

echo $report | jq '.' > $report_name
