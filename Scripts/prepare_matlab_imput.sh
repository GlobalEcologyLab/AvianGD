mkdir ./Alignments/Matlab_fasta

num_loops=`grep -cve '^\s*$' ./Data/dataset_for_matlab.csv`
#mkdir bird_cytb
for ((loops=1;loops<=${num_loops};loops++))
do
    sp_name=`cat ./Data/dataset_for_matlab.csv | head -n ${loops} | tail -n 1 | cut -d ',' -f 2 | sed 's/ /_/g'`
    accnum=`cat ./Data/dataset_for_matlab.csv | head -n ${loops} | tail -n 1 | cut -d ',' -f 1 | sed 's/ /_/g'`
    #geo_info=`cat grid_ids.csv | head -n ${loops} | tail -n 1 | cut -f 5-`
    if [ -n "$(grep $accnum ./Alignments/Checked_Alignments.fasta)" ]; then
        cat ./Alignments/Checked_Alignments.fasta | grep -w -A1 $accnum | head -n 2 | grep -v '>' | sed 's/A/1 /g' | sed 's/C/2 /g' | sed 's/G/3 /g' | sed 's/T/4 /g' | sed 's/-/N/g' | sed 's/[A-Z]/0 /g' >> ./Alignments/Matlab_fasta/${sp_name}.fasta
        #cat all_cytb_fasta.fasta | grep -A1 $accnum | head -n 2 | grep '>' >> ${sp_name}.title #sed 's/A/1 /g' | sed 's/C/2 /g' | sed 's/G/3 /g' | sed 's/T/4 /g' | sed 's/-/N/g' | sed 's/[A-Z]/0 /g' >> ${sp_name}.fasta
        #echo ${geo_info} | sed 's/_/\t/g' >> co1_birds_matlab/${sp_name}.coords
    fi



    echo " $loops of $num_loops "





done
