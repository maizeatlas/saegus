#!/bin/bash


echo "Run ID: $1, Number of Replicates $2"
run_id=$1
number_of_replicates=$2
final_rep_index="$((number_of_replicates - 1))"

echo "Beginning TASSEL analysis of Run ID: $run_id"
echo "Number of Replicates: $number_of_replicates"
echo "First configuration file: small_0_gwas_pipeline.xml"onca
echo "Final configuration file: small_"$final_rep_index"_gwas_pipeline.xml"

for i in `seq 0 $final_rep_index`
do
    config_file_name=$run_id$i"_gwas_pipeline.xml"
    echo "$config_file_name"
    ./run_pipeline.pl -Xmx6g -configFile $config_file_name
done
