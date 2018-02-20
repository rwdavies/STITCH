extra_list=("--gridWindowSize=1000 " "--gridWindowSize=10000 " "--gridWindowSize=100000 ")
version_list=(1.4.0 1.4.0 1.4.0)
name_list=(1.4.0_1000 1.4.0_10000 1.4.0_100000)
options_list=("gridWindowSize=1000" "gridWindowSize=10000" "gridWindowSize=100000")
interface_list=("cli" "cli" "cli")

versions=(1.4.0 1.3.7 1.3.6 1.3.5 1.3.4 1.3.3 1.2.5 1.1.1)
n_e=${#extra_list[@]}
for i_v in $(seq 0 $((${#versions[@]} - 1)))
do
    version=${versions[${i_v}]}
    extra_list[${i_v} + ${n_e}]=" "
    options_list[${i_v} + ${n_e}]="NA"
    version_list[${i_v} + ${n_e}]=${version}
    name_list[${i_v} + ${n_e}]=${versions[${i_v}]}
    if [ "$version" == "1.2.5" ] || [ "$version" == "1.1.1" ]
    then
	interface_list[${i_v} + ${n_e}]="R"
    else
	interface_list[${i_v} + ${n_e}]="cli"	
    fi
done
