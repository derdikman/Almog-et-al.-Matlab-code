function out_mat=Cross_Correlation(mat1,mat2)

[ma,na] = size(mat1);
[mb,nb] = size(mat2);
out_mat = nan(ma+mb-1,na+nb-1);

i_size = size(mat2,1); j_size = size(mat2,2);
[work_mat,npad_i,npad_j] = pad_edges(mat1,mat2,1);

for i = 1:size(out_mat,1)
    for j = 1:size(out_mat,2)
        
        % for each i and j, choose the correct sub-mat (size of mat 2) to
        % multiply with mat2
        
        sub_mat = work_mat(npad_i+i-floor(i_size):npad_i+i-1, ...
            npad_j+j-floor(j_size):npad_j+j-1  );
        nan_sub_mat=sub_mat .* mat2;
        notnan_inds = find(~isnan(nan_sub_mat));  %normalized to the number of nontnan components (average)
        
        n=length(notnan_inds);
        
        if n < 20
            out_mat(i,j) = NaN;
            continue;
        end
        
        sigma_x_y =sum(nan_sub_mat(notnan_inds));
        sigma_x =      sum(sub_mat(notnan_inds));
        sigma_y =      sum(mat2(notnan_inds));
        sigma_x2 = sum(sub_mat(notnan_inds).^2);
        sigma_y2 = sum(mat2(notnan_inds).^2);
        
        out_mat(i,j) = (n*sigma_x_y - sigma_x.*sigma_y) ./ ...
            sqrt(n*sigma_x2-sigma_x.^2) ./ ...
            sqrt(n*sigma_y2-sigma_y.^2);
    end % for j
end % for i
disp('')
end