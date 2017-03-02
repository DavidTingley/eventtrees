function [tree,siz,cumulativeUIDs] = hash_eventtree(varargin)
%% INPUT
%
%   spikes -  1 x N cell array of FMAToolbox formatted spike trains
%   alpha_min  -  double value that is the minimum time over which to search 
%   alpha_max  -  double value that is the maximum time over which to search 
%   seqLength  -  length of sequences we're looking for
%   chainLimit -  the # of spikes we can skip before giving up search (this
%                 indirectly limits recursion depth)
%
%
%% Output
%   tree       -  Map container with UID and observation count. To access
%                 coordinates of a specific sequence, use tree.key(UID)
%                 where UID is returned by the mat2sparse.m function
%                 (formerly z)
%   sparseDims -  repmatted size matrix that is the number of cells recorded,
%                 repeated by the sequence length  essentially the dimmensions
%                 a sparse matrix would be (used for mat2sparse.m
%                 conversion)
%   cumulativeUIDs - N X 1 vector where N is the number of spikes in the
%                    recording and each value is the cumulative number of 
%                    unique sequences up to that spike in the recording
%
%
%% EXTRAS/TODO
%   - Add timestamps as an output variable (first or last spike in sequence)
%   - add argument that is a shank list, so sequences must jump shanks?
%   - other recommendations?
%
% This version of the eventtree's analysis utilizes a matlab Container.map object
% to store sequence observations.  The advantage, being able to look for sequences 
% of arbitrarily long lengths (8-30 realistically) 
%
% D.Levenstein and D.Tingley 2016

warning('this function is experimental and maybe not be functioning (ha) correctly...')
if nargin < 2
    error('not enough inputs');
end

%% rework inputs to be flexible/stable // add debug variable
spikes = varargin{1};
alpha_min = varargin{2};
alpha_max = varargin{3};
file = varargin{6};
if nargin < 4
    seqLength=5;
else
    seqLength = varargin{4};
end
if nargin < 5
    chainlimit = 3;
else
    chainlimit = varargin{5};
end

if size(spikes,1) > size(spikes,2)
    numspikes = length(vertcat(spikes{:}));
else
    numspikes = length(vertcat(spikes{:}));
end  % dafuq does this do?

numcells = length(spikes);

spikesmat = zeros(numspikes,2);
indtrack = 1;
for c = 1:numcells  % throw all spikes into a single two column matrix..
    cellspikes = length(spikes{c});
    spikesmat(indtrack:(indtrack+cellspikes-1),1) = spikes{c};
    spikesmat(indtrack:(indtrack+cellspikes-1),2) = c.*ones(size(spikes{c}));
    indtrack = indtrack + cellspikes;
end

% sort spikesmat by time stamps
[~,sort_stime] = sort(spikesmat(:,1));
spikesmat = spikesmat(sort_stime,:);

% remove duplicate spikes here? check for any?

% set up sparseDimsvariable which mat2sparse.m uses to access the correct UID
sparseDims= repmat(numcells,seqLength,1);

level = seqLength; % set the initial level to that of the sequence length
keySet{1}=0;  % make dummy key/value pair to initialize the tree dictionary
valueSet{1}=0;
tree = containers.Map(keySet,valueSet);

tic

for sp = 1:numspikes-21 % don't search the last ~20 spikes for sequences
    if mod(sp,1000)==0              
%         display(['Finished spike ',num2str(sp),' of ',num2str(numspikes)])
%         display(['Spike ',num2str(sp),' of ',num2str(numspikes), ', ~' num2str((numspikes-sp)./(sp./toc)) ' seconds left..'])
%         
%         subplot(2,2,1)
%         plot(toc,numspikes-sp,'.k')
%         ylabel('# of spikes left')
%         axis([0 toc 0 numspikes])
%         hold on 
%         subplot(2,2,2)
%         plot(toc,(numspikes-sp)./(sp./toc),'.r')
%         hold on
%         ylabel('estimated time left')
%         subplot(2,2,3)
%         plot(sp,tree.Count,'.b')
%         hold on
%         title('numbers of unique sequences')
%         subplot(2,2,4)
        vals = cell2mat(tree.values);
        f = find(vals>0);
        if ~isempty(f)
            [a b]=polyfit(log(1:length(f)),log(sort(vals(f),'descend')),1);
            display([num2str(seqLength) ' ' num2str(a(:,1))])
%             plot(toc,a(:,1),'.')
%             hold on
        end
%         drawnow
        
        save(['eventtree_' file '_' num2str(seqLength) '.mat'],'-v7.3')
    end
    
    sptime = spikesmat(sp,1);
    coords = spikesmat(sp,2);
    if level > 1
        [tree, level_blah] = onedeeper(spikesmat(sp+1:sp+21,:),...
            sptime,alpha_min,alpha_max,tree,numcells,seqLength,sparseDims,coords,level,chainlimit);
    end
    cumulativeUIDs(sp) = tree.Count;
end
end

function [eventtree, level] = onedeeper(spikesmat,sptime,alpha_min,alpha_max,...
    eventtree,numcells,seqLength,sparseDims,coords,level,chainlimit)

chainspikes5 = find((spikesmat(:,1)-sptime)<alpha_max ...
                    & (spikesmat(:,1)-sptime)>alpha_min);
level = level-1;  % take us down a level...
% check here to remove repeats of the same cell (i.e. bursty cells)
chainspikes5(spikesmat(chainspikes5,2)==coords(end)) = [];

if size(chainspikes5,1) > chainlimit  % limit number of spikes we search through
    chainspikes5(chainlimit+1:end) = [];
end

if ~isempty(chainspikes5)
    for sp5 = chainspikes5'
        sptime5 = spikesmat(sp5,1);
        cellNum = spikesmat(sp5,2);
        if ~ismember(cellNum,coords)            
            %need to find a way to remove sequences with the same
            %start/stop timesstampes
        coords_new = [coords cellNum];
%         [z] = mat2sparse(siz,coords_new);
%         eventtree{seqLength-level+1}(z) = eventtree{seqLength-level+1}(z)+1;
        if level > 1
            [eventtree, level_blah] = onedeeper(spikesmat,sptime5,alpha_min,alpha_max,...
                eventtree,numcells,seqLength,sparseDims,coords_new,level,chainlimit);
        else
            % comment these lines out above and uncomment here if you only
            % want to search for sequences of a specific depth
            [z] = mat2sparse(sparseDims,coords_new);
            
            if eventtree.isKey(z)
                eventtree(z) = eventtree(z)+1;
            else
                newMap = containers.Map({z},1);
                eventtree = [eventtree;newMap];
            end
        end
    end
    end
end
end

