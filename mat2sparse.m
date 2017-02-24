function ndx = mat2sparse(siz,coords)

siz = double(siz);
lensiz = length(siz);
if lensiz < 2
    error(message('MATLAB:sub2ind:InvalidSize'));
end

numOfIndInput = length(coords);


    
    %Compute linear indices
    k = [1 cumprod(siz(1:end-1))'];
    ndx = 1;
    for i = 1:numOfIndInput
        v = coords(i);
        if (any(v(:) < 1)) || (any(v(:) > siz(i)))
            %Verify subscripts are within range
            error(message('MATLAB:sub2ind:IndexOutOfRange'));
        end
        ndx = ndx + (v-1)*k(i);
    end

