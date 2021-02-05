%% function for pairwised avg. number of segregating sites per nuc. site
function [avg_segs_bp,tot_num_segs, dropped_pairs] = FunPairSegSites (species_seqs)

num_seqs=size(species_seqs,1);
num_per_bp=[];
num_seg_sites=[];
dropped_pairs=0;
% num_per_bp2=[];
% num_seg_sites2=[];
%num_per_bp3=[];
%num_seg_sites3=[];
for iter=1:num_seqs
    main_seq=species_seqs(iter,:);    
    temp_seqs=species_seqs;
    temp_seqs(iter,:)=[];
    other_seqs=temp_seqs;
    for iter_sub=1:(num_seqs-1)
        pair_seqs=[main_seq;other_seqs(iter_sub,:)];
        common_sites=pair_seqs(:,find((main_seq.*other_seqs(iter_sub,:))>0));
        if size(common_sites,2)<0.5*max(length(find(pair_seqs(1,:)>0)),length(find(pair_seqs(2,:)>0)))
            dropped_pairs=dropped_pairs+1;
            continue;
            %num_per_bp(end+1)=0;
        else
           num_segs=length(find(abs(common_sites(1,:)-common_sites(2,:))>0));
           num_seg_sites(end+1)=num_segs;
           num_per_bp(end+1)=num_segs/size(common_sites,2);
        end
        
%         if size(common_sites,2)<0.5*max(length(find(pair_seqs(1,:)>0)),length(find(pair_seqs(2,:)>0)))
%             %continue;
%             %num_per_bp(end+1)=0;
%         else
%             num_segs2=length(find(abs(common_sites(1,:)-common_sites(2,:))>0));
%             num_seg_sites2(end+1)=num_segs2;
%             num_per_bp2(end+1)=num_segs2/size(common_sites,2);
%         end
        
        %if isempty(common_sites)
            %continue;
        %    num_per_bp3(end+1)=0;
        %    num_seg_sites3(end+1)=0;
        %else
        %    num_segs3=length(find(abs(common_sites(1,:)-common_sites(2,:))>0));
        %    num_seg_sites3(end+1)=num_segs3;
        %    num_per_bp3(end+1)=num_segs3/size(common_sites,2);
        %end
        
    end
    
end

avg_segs_bp=sum(num_per_bp)/length(num_per_bp);
tot_num_segs=sum(num_seg_sites)/2;
dropped_pairs;
% avg_segs_bp2=sum(num_per_bp2)/length(num_per_bp2);
% tot_num_segs2=sum(num_seg_sites2)/2;
%avg_segs_bp3=sum(num_per_bp3)/length(num_per_bp3);
%tot_num_segs3=sum(num_seg_sites3)/2;