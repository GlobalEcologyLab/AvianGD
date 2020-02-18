mkdir matlab_fasta_02dec2019

num_loops=`grep -cve '^\s*$' matlab_input_02Dec2019.csv`
#mkdir bird_cytb
for ((loops=1;loops<=${num_loops};loops++))
do
    sp_name=`cat matlab_input_02Dec2019.csv | head -n ${loops} | tail -n 1 | cut -d ',' -f 2 | sed 's/ /_/g'`
    accnum=`cat matlab_input_02Dec2019.csv | head -n ${loops} | tail -n 1 | cut -d ',' -f 1 | sed 's/ /_/g'`
    #geo_info=`cat grid_ids.csv | head -n ${loops} | tail -n 1 | cut -f 5-`
    if [ -n "$(grep $accnum big_fasta_02Dec2019.fasta)" ]; then
        cat big_fasta_02Dec2019.fasta | grep -w -A1 $accnum | head -n 2 | grep -v '>' | sed 's/A/1 /g' | sed 's/C/2 /g' | sed 's/G/3 /g' | sed 's/T/4 /g' | sed 's/-/N/g' | sed 's/[A-Z]/0 /g' >> matlab_fasta_02dec2019/${sp_name}.fasta
        #cat all_cytb_fasta.fasta | grep -A1 $accnum | head -n 2 | grep '>' >> ${sp_name}.title #sed 's/A/1 /g' | sed 's/C/2 /g' | sed 's/G/3 /g' | sed 's/T/4 /g' | sed 's/-/N/g' | sed 's/[A-Z]/0 /g' >> ${sp_name}.fasta
        #echo ${geo_info} | sed 's/_/\t/g' >> co1_birds_matlab/${sp_name}.coords
    fi
    
    
    
    echo " $loops of $num_loops "
    
    
    
    
    
done
