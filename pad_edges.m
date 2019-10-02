%pads outside of matrix mat, by size of matrix h divided by l
function [out_mat,npad_rows,npad_cols] = pad_edges(mat,h,l)
t = ceil(size(h)/l);
npad_rows=t(1); npad_cols=t(2);
out_size = size(mat) + [2*npad_rows 2*npad_cols];
out_mat = nan(out_size);
out_mat(npad_rows+1:npad_rows+size(mat,1),npad_cols+1:npad_cols+size(mat,2)) = mat;
end

%{
%original code
function [out_mat,npad_i,npad_j] = pad_edges(mat,h,l)
npad_ij = ceil(size(h)/l);
npad_i = npad_ij(1);
npad_j = npad_ij(2);
in_size = size(mat);
out_size = in_size + [2*npad_i 2*npad_j];
out_mat = nan(out_size);
out_mat(npad_i+1:npad_i+in_size(1),npad_j+1:npad_j+in_size(2)) = mat;
disp('')
end
%}