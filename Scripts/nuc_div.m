%% put seqs in grid cell
clear;
%fid=fopen('species_list_coi','r');
fid=fopen('list_species_for_matlab','r');
c=textscan(fid,'%s');
fclose(fid);
species_list=c{1,:};
num_files=size(species_list,1);
%coords_real=[];
mammals_biome={};
for iter_file=1:num_files        
    file_name_seq=strcat('Matlab_fasta/',cell2mat(species_list(iter_file)),'.fasta');
    %file_name_coords=strcat('cytb_birds_matlab/',cell2mat(species_list(iter_file)),'.coords');
    species_seqs=load(file_name_seq);
    
    length_seq=size(species_seqs,2);
    num_seqs2=size(species_seqs,1);
    %species_coords=load(file_name_coords);
%    coords_biome=[ceil((species_coords(:,1)+2)/4)*4-4,ceil(species_coords(:,2)/4)*4-2]; % 4*4 biome cell
    %species_biome=species_coords(:,3);
    %uniq_biome=unique(species_biome,'rows');
    %coords_biome=2*(ceil(species_coords/2)-0.5); % 2*2 biome cell
    %uniq_biome=unique(coords_biome,'rows');
    %for iter_biome=1:size(uniq_biome,1)
        %seqs_biome=species_seqs(find(ismember(species_biome,uniq_biome(iter_biome,:),'rows')==1),:);
        %num_seqs_biome=size(seqs_biome,1);       
        %seg_sites_biome=FunSegSites(seqs_biome);
     if num_seqs2>1
        [avg_pair,pair_muts, dropped_pairs]=FunPairSegSites(species_seqs);
     else
         pair_muts=0;
         avg_pair=0;
            %pair_muts2=0;
            %avg_pair2=0;
            %pair_muts3=0;
            %avg_pair3=0;
     end
     %num_seg_sites_biome=size(seg_sites_biome,2);
     %total_valid_sites=length(find(reshape(seqs_biome,1,[]))>0);
     mammals_biome(end+1,:)=[species_list(iter_file) avg_pair pair_muts dropped_pairs num_seqs2] ;
                                         %uniq_biome(iter_biome,1) ... 
                                         %num_seqs_biome ... 
                                         %avg_pair pair_muts ...
                                         %total_valid_sites];
    %end
    iter_file
end
