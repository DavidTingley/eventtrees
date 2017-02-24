function cellseq = sparse2mat(siz,ndx)
%Input
%   siz: repmat(numcells,nmax,1)
%   ndx: linear index

nout = length(siz);
siz = double(siz);
    k = [1 cumprod(siz(1:end-1))'];
    for i = nout:-1:1,
        vi = rem(ndx-1, k(i)) + 1;
        vj = (ndx - vi)/k(i) + 1;
        cellseq{i} = vj;
        ndx = vi;
    end
    
    cellseq = cell2mat(cellseq);

