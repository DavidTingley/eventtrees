function [] = permutationTestEventTree(tree,siz(1:level))
  %% this function is designed to take the output from EventTree.m (or the 
  %% hash table version) and return sequences that pass a series of 
  %% criteria that are:
  %% 1) permutation/actual ratio greater than XXX
  %% 2) observation rate higher than XXX
  %% 3) 
level = 5;

f = find(tree{level}>100);
for i=1:length(f)
    cellseq(i,:) = sparse2mat(siz(1:level),f(i));
    idx(i) = f(i);
end
disp(['starting with ' num2str(length(f)) ' potential sequences..'])


obs  = zeros(length(idx),length(perms(cellseq(1,:))));
% te_summed = zeros(length(idx),length(perms(cellseq(1,:))));
% te_actual = zeros(length(idx),1);

for s = 1:length(idx)
   seq = cellseq(s,:);
   perm = perms(seq);
   for p=1:length(perm)
       ndx = mat2sparse(siz(1:level),perm(p,:));
%        if ndx ~= idx(s)
%        [ts] = findchain(spikes,perm(p,:),alpha,chainlimit);
          obs(s,p) = tree{level}(ndx);
%        end
% t=zeros(24,1);
% for i = 1:24
%     for j=1:length(seq)-1
%         te_summed(s,p) = te_summed(s,p) + te_result(perm(p,j),perm(p,j+1));
%     end
% end

   end
   actual(s) = tree{level}(idx(s));
%    for j=1:length(seq)-1
%         te_actual(s) = te_actual(s) + te_result(seq(j),seq(j+1));
%    end
end

f = find(actual>100);
ff = find((actual./sum([obs'; actual]))>.3);

%% shuffling here


%% prediction here


[aic cellseq logl dev] = sequence_peer_prediction(ts,spikes,cellseq,seqLength)

    


return