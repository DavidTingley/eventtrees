function [ts,numoccurance,...
    ts_nonoverlap,numoccurance_nonoverlap] = findchain(spiketimes,cellseq,alpha,chainlimit)

%INPUT
%     spiketimes    cell array where each element is the time stamps of a single
%                   neuron
%     cellseq       [nsequences x mmax] matrix of the sequences we want to
%                   hunt for. each row is a sequence of cell numbers.
%     alpha         time window to constrain each transition to
%     chainlimit    number of spikes ahead to look - saves time, but be careful!
%                   
%
%OUTPUT
%   ts              {nsequences x 1} cell array with time stamps for event 
%                   timesof each sequence, each entry is a vector constining  
%                   spike time of the [FIRST LAST] spike in the sequence.
%   numoccurance    [nsequences] vector of the counts of how many times 
%                   each sequence was observed
%
%
%D.Levenstein and D.Tingley 2016/17
%%
if ~exist('overlap','var')
    overlap = true;
end

%Make numspikes x 2 vector, (:,1) is spiketimes (:,2) is cell numbers
if iscell(spiketimes)
    
    if size(spiketimes,1) > size(spiketimes,2)
        numspikes = length(vertcat(spiketimes{:}));
    else
        numspikes = length(vertcat(spiketimes{:}));
    end
    numcells = length(spiketimes);
   
    spikesmat = zeros(numspikes,2);
    indtrack = 1;
    for c = 1:numcells
        cellspikes = length(spiketimes{c});
        spikesmat(indtrack:(indtrack+cellspikes-1),1) = spiketimes{c};
        spikesmat(indtrack:(indtrack+cellspikes-1),2) = c.*ones(size(spiketimes{c}));
        indtrack = indtrack + cellspikes;
    end
elseif isnumeric(spiketimes) && size(spiketimes,2)==2
    spikesmat = spiketimes;
    numspikes = size(spikesmat,1);
    numcells = length(unique(spikesmat(:,2)));
else
    display('Something is wrong with you spiketimes input')
end

%sort spikesmat
[~,sort_stime] = sort(spikesmat(:,1));
spikesmat = spikesmat(sort_stime,:);

cellnums = unique(spikesmat(:,2));
maxcell = max(cellnums);

% set up siz variable which mat2sparse.m uses to access the correct sparse
% coordinates
mmax = size(cellseq,2);
%     eventtree = sparse(m,1); clear m
siz=repmat(maxcell,mmax,1);  %This doesn't get used...

level = mmax; % set the inital level to that of the sequence depth
% set sequence counts to zeros
counts=zeros(1+mmax,1);
tic

numseqs = size(cellseq,1);
ts = cell(numseqs,1);
numoccurance = zeros(numseqs,1);


%% This is the top level of the recursive loop
%Start with each spike that matches the first cell in the sequence(s)
firstspikepossibilities = unique(cellseq(:,1));
firstspikeindices = find(ismember(spikesmat(:,2),firstspikepossibilities));
for sp = firstspikeindices'
    
    if mod(sp,10000)==0              
        timesofar = toc; percdone = sp./numspikes;
        totaltimeestimate = timesofar./percdone;
        timeleft = totaltimeestimate-timesofar;
        display(['Spike ',num2str(sp),' of ',num2str(numspikes), ...
            '.  Est. Total Time: ',num2str(round(totaltimeestimate./60,1)),...
            'min.  ETA: ',num2str(round(timeleft./60,1)),'min.'])
    end
    
    %Get spike time/cell number of the first spike
    sptime = spikesmat(sp,1);
    firstsptime = sptime;
    cellNum = spikesmat(sp,2);
    
    %Pass on the sequences that start with this cell
    seqnum = find(cellseq(:,1)==cellNum);
    keepseqs = cellseq(seqnum,:);
    
    %Go to the next spike in the tree
    [ts, level_blah] = onedeeper(spikesmat,sp,...
        sptime,alpha,ts,maxcell,mmax,siz,level,...
        keepseqs,firstsptime,chainlimit,seqnum);
end

%% This is filtering to remove overlaps

%This could be faster if combined into one cellfun...
[~,uniquestartends] = cellfun(@(X) unique(X(:,1)),ts,'UniformOutput',false);
[~,uniquestartends] = cellfun(@(X,Y) unique(X(Y,2)),ts,uniquestartends,'UniformOutput',false);
ts_nonoverlap = cellfun(@(X,Y) X(Y,:),ts,uniquestartends,'UniformOutput',false);

numoccurance_nonoverlap = zeros(size(numoccurance));
%Count the occurances (make this a cellfun), don't need to count 0's as
%they're already 0's
for ss = 1:numseqs
    if isempty(ts{ss})
        numoccurance(ss) = 0;
        numoccurance_nonoverlap(ss) = 0;
    else
        numoccurance(ss) = length(ts{ss}(:,1));
        numoccurance_nonoverlap(ss) = length(ts_nonoverlap{ss}(:,1));
    end
end


end


%%
%Recursive function that takes a spike index and looks at the next level of
%spikes that fit alpha and chainlimit criteria
function [ts, level] = onedeeper(spikesmat,sp,sptime,alpha,...
    ts,maxcell,mmax,siz,level,cellseq,firstsptime,chainlimit,seqnum)

    %Find the indices of spikes within alpha seconds (or chainlimit spikes) in the future
    timelag = spikesmat(sp+1:min(end,sp+chainlimit),1)-sptime;
    chainspikes = find(timelag<alpha & timelag>0)+sp;
    
    %Manage which level you're in
    level = level-1;  % take us down a level...
    treedepth = mmax-level+1; %ex: treedepth=2 for the second spike in a sequence

    %Keep only next spikes that match the next spike in the desired sequence
    nextspikepossibilities = unique(cellseq(:,treedepth));
    nextspikeindices = chainspikes(ismember(spikesmat(chainspikes,2),nextspikepossibilities));
    %If there's no cells next, end the tree.
    if isempty(nextspikeindices); return; end
    for sp5 = nextspikeindices'
        
        %Get spike time/cell number of the next spike
        sptime5 = spikesmat(sp5,1);
        cellNum = spikesmat(sp5,2);
        
        %Pass on the sequences that continue the chain with this cell
        keepseqs5 = cellseq(cellseq(:,treedepth)==cellNum,:);
        seqnum5 = seqnum(cellseq(:,treedepth)==cellNum);
                
        if level > 1
            %If you're not at the end of the chain, go to the next spike in the tree
            [ts, level_blah] = onedeeper(spikesmat,sp5,...
                sptime5,alpha,ts,maxcell,mmax,siz,level,...
                keepseqs5,firstsptime,chainlimit,seqnum5);
        elseif level == 1
            %If you make it to the last spike in the chain, add the first
            %and last spike time in the sequence to the output vector for
            %the correct cellseq
            ts{seqnum5} = [ts{seqnum5}; firstsptime, sptime5];          
        else display('Welcome to Bugland'); keyboard;
        end
    end
    
end

