extra_list=("--output_format='bgen' " "--gridWindowSize=1000 " "--gridWindowSize=10000 " "--gridWindowSize=100000 ")
version_list=(1.5.2 1.5.2 1.5.2 1.5.2)
name_list=(1.5.2_bgen 1.5.2_1000 1.5.2_10000 1.5.2_100000)
options_list=("output_format=\'bgen\'" "gridWindowSize=1000" "gridWindowSize=10000" "gridWindowSize=100000")
interface_list=("cli" "cli" "cli" "cli")
extension_list=("bgen" "vcf.gz" "vcf.gz" "vcf.gz")
    
versions=(1.5.2 1.5.1 1.5.0 1.4.2 1.4.1 1.4.0 1.3.7 1.3.6 1.3.5 1.3.4 1.3.3 1.2.5 1.1.1)
n_e=${#extra_list[@]}
for i_v in $(seq 0 $((${#versions[@]} - 1)))
do
    version=${versions[${i_v}]}
    extra_list[${i_v} + ${n_e}]=" "
    options_list[${i_v} + ${n_e}]="NA"
    version_list[${i_v} + ${n_e}]=${version}
    extension_list[${i_v} + ${n_e}]="vcf.gz"
    name_list[${i_v} + ${n_e}]=${versions[${i_v}]}
    if [ "$version" == "1.2.5" ] || [ "$version" == "1.1.1" ]
    then
	interface_list[${i_v} + ${n_e}]="R"
    else
	interface_list[${i_v} + ${n_e}]="cli"	
    fi
done
