function [eventtree,counts] = EventTree_test(varargin)
%% INPUT
%
%
%
%
%% Output
%   eventtree   (mmax*numcells.^(mmax)) x 1 sparse matrix. To access
%               coordinates of a specific sequence, use sparse2mat function to pull out
%               the coordinates(cell #s) of a sequence
%
%   counts      D x 1 vector that is the number of total sequences observed
%               at that depth (Depth is the number of spikes in a sequence)
%
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
%

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

% set up siz variable which mat2sparse.m uses to access the correct sparse
% coordinates

for i = 1:mmax
    m = numcells.^(i);
    eventtree{i} = sparse(m,1);
end
%     eventtree = sparse(m,1); clear m
siz=repmat(numcells,mmax,1);

level = mmax; % set the inital level to that of the sequence depth
% set sequence counts to zeros
counts=zeros(1+mmax,1);
tic

for sp = 1:numspikes
    if mod(sp,1000)==0
        display(['Spike ',num2str(sp),' of ',num2str(numspikes)])
        toc
    end
    
    
    sptime = spikesmat(sp,1);
    
    cellNum = spikesmat(sp,2);
    coords = [cellNum];
    [z] = mat2sparse(siz,coords);
    eventtree{mmax-level+1}(z) = eventtree{mmax-level+1}(z)+1;
    
    chainspikes = find((spikesmat(sp:min(end,sp+1000),1)-sptime)<alpha...
        & (spikesmat(sp:min(end,sp+1000),1)-sptime)>0)+sp-1;
    
    
    counts(mmax-level+1)=counts(mmax-level+1)+1;
    
    if level > 1
        [eventtree, level_blah,counts] = onedeeper(counts,spikesmat,sp,...
            sptime,alpha,eventtree,chainspikes,numcells,mmax,siz,coords,level,chainlimit);
    end
end
end

function [eventtree, level,counts] = onedeeper(counts,spikesmat,sp,sptime,alpha,...
    eventtree,chainspikes,numcells,mmax,siz,coords,level,chainlimit)
%     disp(['running level: ' num2str(level)])  % make this accessible on
%     debug?
chainspikes5 = find((spikesmat(sp:min(end,sp+1000),1)-sptime)<alpha & (spikesmat(sp:min(end,sp+1000),1)-sptime)>0)+sp-1;
level = level-1;  % take us down a level...
siz=repmat(numcells,mmax-level+1,1);
% check here to remove repeats of the same cell (i.e. bursty cells)
repeats = find(spikesmat(chainspikes5,2)==coords(end));
chainspikes5(repeats) = [];

if size(chainspikes5,1) > chainlimit  % limit number of spikes we search through
    chainspikes5(chainlimit+1:end) = [];
end

if ~isempty(chainspikes5)
    for sp5 = chainspikes5'
        sptime5 = spikesmat(sp5,1);
        cellNum = spikesmat(sp5,2);
        coords_new = [coords cellNum];
        
        [z] = mat2sparse(siz,coords_new);
        eventtree{mmax-level+1}(z) = eventtree{mmax-level+1}(z)+1;
        counts(mmax-level+1)=counts(mmax-level+1)+1;
        
      
        
        if level > 1
            [eventtree, level_blah,counts] = onedeeper(counts,spikesmat,sp5,sptime5,alpha,...
                eventtree,chainspikes5,numcells,mmax,siz,coords_new,level,chainlimit);
        else
            
            % comment these lines out above and uncomment here if you only
            % want to search for sequences of a specific depth
            
            %             [z] = mat2sparse(siz,coords_new);
            %             eventtree{mmax-level+1}(z) = eventtree{mmax-level+1}(z)+1;
            %             counts(mmax-level+1)=counts(mmax-level+1)+1;
        end
        
    
        
    end
end
end
