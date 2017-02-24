function [tree,counts,siz] = EventTree(varargin)
%[tree,counts,siz] = EventTree(varargin)
%
%INPUTS
%   spiketimes  {ncells} cell array of [nspikes] vector of spike times for each cell
%   alpha       number of seconds to look ahead after each spike
%   mmax        longest tree
%   chainlimit  recursion limit - number of spikes to look ahead before
%               giving up
%   treetype    'all', 'noburst', 'norepeat'   %BUG%
%
%OUTPUTS
%   tree        (mmax*numcells.^(mmax)) x 1 sparse matrix. To access
%               coordinates of a specific sequence, use sparse2mat function to pull out
%               the coordinates(cell #s) of a sequence
%   counts      D x 1 vector that is the number of total sequences observed
%               at that depth (Depth is the number of spikes in a sequence)
%               Note: counts includes repeats
%
%
%D.Tingley and D.Levenstein 2016
%% EXTRAS/TODO
%   -Add timestamps as an output variable (first or last spike in sequence)
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

% to get seq ID:
% cellseq = sparse2mat(repmat(numcells,mmax,1),seqID)

%add argument that is a shank list, so sequences must jump shanks?
%
%%
if nargin < 2
    error('not enough inputs');
end
spiketimes = varargin{1};
alpha = varargin{2};

if nargin < 3
    mmax=5;
else
    mmax = varargin{3};
end
if nargin < 4
    chainlimit = 3;
else
    chainlimit = varargin{4};
end

if nargin < 5
    treetype = 'norepeat';
else
    treetype = varargin{5};
end

if size(spiketimes,1) > size(spiketimes,2)
    numspikes = length(vertcat(spiketimes{:}));
else
    numspikes = length(vertcat(spiketimes{:}));
end
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

% remove duplicate spikes here


% set up siz variable which mat2sparse.m uses to access the correct sparse
% coordinates

siz = repmat(numcells,mmax,1);

level = mmax; % set the inital level to that of the sequence depth
% set sequence counts to zeros

for i = 1:mmax
    m = numcells.^(i);
    tree{i} = sparse(m,1);
    counts(i) = 0;
end
%    tree = sparse(0);
%    counts = [];

tic



for sp = 1:numspikes
%     [tree, counts] = onespike(sp,numspikes,spikesmat,level,mmax,siz,numcells,alpha,chainlimit);

    %Plotting and Time Counter
    if mod(sp,10000)==0              
        timesofar = toc;
        percdone = sp./numspikes;
        totaltimeestimate = timesofar./percdone;
        timeleft = totaltimeestimate-timesofar;
        display(['Spike ',num2str(sp),' of ',num2str(numspikes), ...
            '.  Est. Total Time: ',num2str(round(totaltimeestimate./60,1)),...
            'min.  ETA: ',num2str(round(timeleft./60,1)),...
            'min.'])
        
        subplot(2,2,1)
            plot(toc,numspikes-sp,'.k')
            axis([0 toc 0 numspikes])
            hold on
            xlabel('Time (s)');ylabel('Num. Spks. Left')
        subplot(2,2,2)
            plot(toc,(numspikes-sp)./(sp./toc),'.r')
            hold on
            xlabel('Time');
        subplot(2,2,3)
            hold off
            for i=1:mmax
                f = find((tree{i}>0));
                [a(i,:) b{i}]=polyfit(log(1:length(f)),log((sort(tree{i}(f),'descend')')),1);
                lens(i)=length(f);
                plot(log10(1:length(f)),log10(sort(tree{i}(f),'descend')))
                hold on
            end
            LogScale('xy',10)
            xlabel('Sequence Rank Order');ylabel('# Occurances')
        subplot(2,2,4)
            plot(a(:,1))
            title(num2str(lens))
            xlabel('Tree Length');ylabel('Slope')
        drawnow
%         save('temp_tree.mat','tree','counts','-v7.3')
    end
    
    %Get spike time/cell number of the first spike
    sptime = spikesmat(sp,1);    
    cellNum = spikesmat(sp,2);
    
    %Add the 1spike chain to the count
    coords = [cellNum];
    [z] = mat2sparse(siz,coords);
    tree{mmax-level+1}(z) = tree{mmax-level+1}(z)+1;
    counts(mmax-level+1)=counts(mmax-level+1)+1;
    
    %Go to the next spike in the tree
    if level > 1
        [tree, level_blah,counts] = onedeeper(counts,spikesmat,sp,...
            sptime,alpha,tree,numcells,mmax,coords,level,chainlimit,treetype);
    end
        
end

display('DONE!')
end

%Recursive function that takes a spike index and looks at the next level of
%spikes that fit alpha and chainlimit criteria
function [eventtree, level,counts] = onedeeper(counts,spikesmat,sp,sptime,alpha,...
    eventtree,numcells,mmax,coords,level,chainlimit,treetype)

    %Find the indices of spikes within alpha seconds (or chainlimit spikes) in the future
    timelag = spikesmat(sp+1:min(end,sp+chainlimit),1)-sptime;
    chainspikes = find(timelag<alpha & timelag>0)+sp;
    
    %Remove repeats of the same cell (i.e. bursty cells), if applicable
    switch treetype
        case 'noburst'
            repeats = ismember(spikesmat(chainspikes,2),coords(end));
            chainspikes(repeats) = [];
        case 'norepeat'
            repeats = ismember(spikesmat(chainspikes,2),coords);
%             if isequal([36,11,5,1],coords) %BUG
%                 keyboard
%             end
            chainspikes(repeats) = [];
    end
    
    %If there's no cells next, end the tree.
    if isempty(chainspikes); return; end
    
    %Manage which level you're in
    level = level-1;  % take us down a level...
    treedepth = mmax-level+1; %ex: treedepth=2 for the second spike in a sequence
    siz=repmat(numcells,treedepth,1);
    
    %Loop through all the next spikes in the chain
    for sp5 = chainspikes'
        if isempty(chainspikes); keyboard; end
        %Get spike time/cell number of the next spike
        sptime5 = spikesmat(sp5,1);
        cellNum = spikesmat(sp5,2);
        
        %Add the new spike chain to the count
        coords_new = [coords cellNum];
        [z] = mat2sparse(siz,coords_new);
        eventtree{treedepth}(z) = eventtree{treedepth}(z)+1;
        counts(treedepth)=counts(treedepth)+1;
        
        %If you're not yet at mmax, go on to the next spikes in the chain
        if level > 1
            [eventtree, level_blah,counts] = onedeeper(counts,spikesmat,sp5,sptime5,alpha,...
                eventtree,numcells,mmax,coords_new,level,chainlimit,treetype);
        else
            
            % comment these lines out above and uncomment here if you only
            % want to search for sequences of a specific depth
            
            %             [z] = mat2sparse(siz,coords_new);
            %             eventtree{mmax-level+1}(z) = eventtree{mmax-level+1}(z)+1;
            %             counts(mmax-level+1)=counts(mmax-level+1)+1;
        end

    end

end



%Is this function still needed?
% function [eventtree counts] = onespike(sp,numspikes,spikesmat,level,mmax,siz,numcells,alpha,chainlimit)  
%     for i = 1:mmax
%         m = numcells.^(i);
%         eventtree{i} = sparse(m,1);
%     end
%     counts=zeros(1+mmax,1);
%     
% %     if mod(sp,1000)==0
% %         display(['Spike ',num2str(sp),' of ',num2str(numspikes) ', '  ])
% % %         toc
% %     end
%     sptime = spikesmat(sp,1);
%     cellNum = spikesmat(sp,2);
%     coords = [cellNum];
%     [z] = mat2sparse(siz,coords);
%     eventtree{mmax-level+1}(z) = eventtree{mmax-level+1}(z)+1;
%     chainspikes = find((spikesmat(sp:min(end,sp+1000),1)-sptime)<alpha...
%         & (spikesmat(sp:min(end,sp+1000),1)-sptime)>0)+sp-1;
%     
%     counts(mmax-level+1)=counts(mmax-level+1)+1;
%     if level > 1
%         [eventtree, level_blah,counts] = onedeeper(counts,spikesmat,sp,...
%             sptime,alpha,eventtree,chainspikes,numcells,mmax,siz,coords,level,chainlimit);
%     end
% end


