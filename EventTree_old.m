function [eventtree] = EventTree(spiketimes,alpha,mmax)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%Output
%   eventtree{m} is the numneurons^m array set of m-event trees.
%   eventtree{m}(j1,j2,...jm) is the histogram of event chains
%   j1 -> j2 -> ... -> jm
%
%TO DO
%   -Loop
%   -add nmax
%   -clean up: order spikesmat by time, only look at entries before current
%   spike and ~1000 after?
%   -make sparse matrix
%
%%
%Load the UP state spiking data
%Load the UP state spiking data
% load('Database/BWRat19_032513_SSubtypes.mat')
% 
% timelim = 100; %only look at first hour of recording
% 
% numcells = length(Se);
% for c = 1:numcells
%     CellSpikes{c} = Range(Se{c},'s');
%     CellSpikes{c} = CellSpikes{c}(CellSpikes{c}<timelim);
% end
% 
% spiketimes = CellSpikes;
% mmax = 5;
% alpha = 0.015;

%%

numspikes = length(vertcat(spiketimes{:}));
numcells = length(spiketimes);

%Make numspikes x 2 vector, (:,1) is spiketimes (:,2) is cell numbers
spikesmat = zeros(numspikes,2);
indtrack = 1;
for c = 1:numcells
    cellspikes = length(spiketimes{c});
    spikesmat(indtrack:(indtrack+cellspikes-1),1) = spiketimes{c};
    spikesmat(indtrack:(indtrack+cellspikes-1),2) = c.*ones(size(spiketimes{c}));
    indtrack = indtrack + cellspikes;
end

%sort spikesmat
[~,sort_stime] = sort(spikesmat(:,1));
spikesmat = spikesmat(sort_stime,:);
%%
% eventtree{1} = sparse([],[],[],numcells,1);
% eventtree{2} = repmat(eventtree{1},[1,numcells]);
% eventtree{3} = repmat(eventtree{2},[1,1,numcells]);
% 
% %Error using repmat
% %N-D sparse output arrays are not supported.
% %bummer... use a large diagonalized sparse matrix instead
% eventtree{2} = sparse([],[],[],numcells^1,numcells^1);
% eventtree{3} = sparse([],[],[],numcells^2,numcells^2);
% eventtree{4} = sparse([],[],[],numcells^3,numcells^3);
% 
% %idea from david
% %eventtree{3}(j1,j2,j3) =>
% %eventtree{3}((j1+(j3-1)*numcells),(j2+(j3-1)*numcells))
% %etc...

%%
eventtree{1} = zeros(numcells,1);
eventtree{2} = zeros(numcells,numcells);
eventtree{3} = zeros(numcells,numcells,numcells);
eventtree{4} = zeros(numcells,numcells,numcells,numcells);
eventtree{5} = zeros(numcells,numcells,numcells,numcells,numcells);


j = zeros(mmax,1);
for sp = 1:numspikes
    if mod(sp,1000)==0
        display(['Spike ',num2str(sp),' of ',num2str(numspikes)])
    end
    
    sptime = spikesmat(sp,1);
    %1 event chain
    j(1) = spikesmat(sp,2);
    eventtree{1}(j(1)) = eventtree{1}(j(1))+1;
    
    %2 event chains
    chainspikes2 = find((spikesmat(sp:min(end,sp+1000),1)-sptime)<alpha & (spikesmat(sp:min(end,sp+1000),1)-sptime)>0)+sp-1;
    for sp2 = chainspikes2'
        
        sptim2 = spikesmat(sp2,1);
        j(2) = spikesmat(sp2,2);
        eventtree{2}(j(1),j(2)) = eventtree{2}(j(1),j(2))+1;
        
        %3 event chains
        chainspikes3 = find((spikesmat(sp2:min(end,sp2+1000),1)-sptim2)<alpha & (spikesmat(sp2:min(end,sp2+1000),1)-sptim2)>0)+sp2-1;
        for sp3 = chainspikes3'
            sptim3 = spikesmat(sp3,1);
            j(3) = spikesmat(sp3,2);
            eventtree{3}(j(1),j(2),j(3)) = eventtree{3}(j(1),j(2),j(3))+1;
            
            %4 event chains
            chainspikes4 = find((spikesmat(sp3:min(end,sp3+1000),1)-sptim3)<alpha & (spikesmat(sp3:min(end,sp3+1000),1)-sptim3)>0)+sp3-1;
            for sp4 = chainspikes4'
                sptim4 = spikesmat(sp4,1);
                j(4) = spikesmat(sp4,2);
                eventtree{4}(j(1),j(2),j(3),j(4)) = eventtree{4}(j(1),j(2),j(3),j(4))+1;
                
                %5 event chains
                chainspikes5 = find((spikesmat(sp4:min(end,sp4+1000),1)-sptim4)<alpha & (spikesmat(sp4:min(end,sp4+1000),1)-sptim4)>0)+sp4-1;
                for sp5 = chainspikes5'
                    sptim5 = spikesmat(sp5,1);
                    j(5) = spikesmat(sp5,2);
                    eventtree{5}(j(1),j(2),j(3),j(4),j(5)) = eventtree{5}(j(1),j(2),j(3),j(4),j(5))+1;

                end
            end
        end
    end
end
        


end

