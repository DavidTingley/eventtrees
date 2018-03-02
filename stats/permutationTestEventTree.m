function [] = permutationTestEventTree(tree,siz,spiketimes)
  %% this function is designed to take the output from EventTree.m (or the 
  %% hash table version) and return sequences that pass a series of 
  %% criteria that are:
  %% 1) permutation/actual ratio greater than XXX
  %% 2) observation rate higher than XXX
  %% 3) 
level = 3;

f = find(tree{level}>50);
[a b] = (sort(tree{level}(f),'descend'));
for i=1:length(f)
    cellseq(i,:) = sparse2mat(siz(1:level),f(b(i)));
    idx(i) = f(b(i));
end
disp(['starting with ' num2str(length(f)) ' potential sequences..'])


obs  = zeros(length(idx),length(perms(cellseq(1,:))));
actual  = zeros(length(idx),1);
% te_summed = zeros(length(idx),length(perms(cellseq(1,:))));
% te_actual = zeros(length(idx),1);

for s = 1:length(idx)
   seq = cellseq(s,:);
   perm = perms(seq);
   f = find(actual>100);
   ff = find((actual./sum([obs, actual],2))>.3);
   subplot(2,2,3)
   plot(s,((actual(s)./sum([obs(s,:), actual(s)]))),'.k')
   actual(s) = tree{level}(idx(s));

   for p=1:length(perm)
       ndx = mat2sparse(siz(1:level),perm(p,:));
%        if ndx ~= idx(s)
          obs(s,p) = tree{level}(ndx);
%        end
% t=zeros(24,1);
% for i = 1:24
%     for j=1:length(seq)-1
%         te_summed(s,p) = te_summed(s,p) + te_result(perm(p,j),perm(p,j+1));
%     end
% end
   end
   for p=1:length(perm)
     if ((actual(s)./sum([obs(s,:), actual(s)])))>.2
        [ts,numoccurance,...
        ts_nonoverlap,numoccurance_nonoverlap] = findchain(spiketimes,perm(p,:),.22,chainlimit);

        subplot(2,2,1);

        [times groups] = spikes2sorted({ts{1}(:,2),spikes.times{44}});
        [ccg t] = CCG(times,groups,'binSize',.001,'duration',.2);
        plot(t,ccg(:,1,2)./obs(s,p),'r')
        hold on
        subplot(2,2,2);
        [times groups] = spikes2sorted({ts{1}(:,2),spikes.times{43}});
        [ccg t] = CCG(times,groups,'binSize',.001,'duration',.2);
        plot(t,ccg(:,1,2)./obs(s,p),'r')
        hold on
     end
   end
    if ((actual(s)./sum([obs(s,:), actual(s)])))>.2
    [ts,numoccurance,...
    ts_nonoverlap,numoccurance_nonoverlap] = findchain(spiketimes,cellseq(s,:),.22,chainlimit);
   
    subplot(2,2,1);
    hold on
    [times groups] = spikes2sorted({ts{1}(:,2),spikes.times{44}});
    [ccg t] = CCG(times,groups,'binSize',.001,'duration',.2);
    plot(t,ccg(:,1,2)./actual(s),'k')
    hold off

    subplot(2,2,2);
    hold on
    [times groups] = spikes2sorted({ts{1}(:,2),spikes.times{43}});
    [ccg t] = CCG(times,groups,'binSize',.001,'duration',.2);
    plot(t,ccg(:,1,2)./actual(s),'k')
    hold off
    
    %    for j=1:length(seq)-1
%         te_actual(s) = te_actual(s) + te_result(seq(j),seq(j+1));
%    end

    subplot(2,2,3)
    plot(s,((actual(s)./sum([obs(s,:), actual(s)]))),'.k')
    pause
    clf
    hold on
    
    
    end
   


%% shuffling here


%% prediction here


% [aic cellseq logl dev] = sequence_peer_prediction(ts,spiketimes,cellseq,seqLength)

end
return