function xc2mat=xcorr2g(A,B)%(mat1,mat2)
if nargin==1
    B=A;
end
[mA,nA] = size(A);
[mB,nB] = size(B); 
xc2mat = nan(mA+mB-1,nA+nB-1);
%isize = size(B,1); jsize = size(B,2);
[workmat,npadi,npadj] = pad_edges(A,B,1); %returns A padded by size of B with nans

ioff = npadi-floor(mB):npadi-1; 
joff = npadj-floor(nB):npadj-1;
B=B(:);
il= size(xc2mat,1); jl = size(xc2mat,2);
%for each i and j, choose the correct sub-matrix of A (size of B) to multiply with B
for i = 1:il
    for j = 1:jl %get correct submatrix       
        %subA = workmat(npadi+i-floor(mB):npadi+i-1 , npadj+j-floor(nB):npadj+j-1 );
        SA = workmat(ioff+i,joff+j); 
        SA=SA(:);
        SAtimesB = SA .* B;
        notnanids = ~isnan(SAtimesB);
        n=sum(notnanids);
        
        if n<20; xc2mat(i,j) = NaN; continue; end
        
        SAn = SA(notnanids);
        Bn  =  B(notnanids);
        
        SAdotB =    sum( SAtimesB(notnanids));
        sumSA =     sum(      SAn           );
        sumB =      sum(       Bn           );
        sumSAsqrd = sum(      SAn.^2        );
        sumBsqrd =  sum(       Bn.^2        );
        
        xc2mat(i,j) = (     n*SAdotB - sumSA.*sumB) ./ ...
                      ( sqrt( n*sumSAsqrd - sumSA.^2 ) .* ...
                        sqrt( n* sumBsqrd - sumB .^2 ) ); 
    end % for j
end % for i

end
%{

a=cellsn(2); a = a.before.rm;

tic; for g=[1:1]; b=Cross_Correlation(a,a); end; toc;
tic; for g=[1:1]; c=xcorr2g(a,a); end; toc;
tic; for g=[1:1]; d=xcorr2go(a,a); end; toc;

figure(1); clf;
imgsc(xcorr2g(a,a));
imgsc(xcorr2go(a,a));

%}


% function out_mat=PointCorr(Rxx,RxxR)
%     nan_mat=Rxx .* RxxR;
%     notnan_inds = find(~isnan(nan_mat));  %normalized to the number of nontnan components (average)
%     n=length(notnan_inds);
%     if n < 2
%         out_mat = NaN;
%     end
%     sigma_x_y =sum(nan_mat(notnan_inds));
%     sigma_x =      sum(Rxx(notnan_inds));
%     sigma_y =      sum(RxxR(notnan_inds));
%     sigma_x2 = sum(Rxx(notnan_inds).^2);
%     sigma_y2 = sum(RxxR(notnan_inds).^2);
% 
%     out_mat = (n*sigma_x_y - sigma_x.*sigma_y) ./ ...
%         sqrt(n*sigma_x2-sigma_x.^2) ./ ...
%         sqrt(n*sigma_y2-sigma_y.^2);
% end

%         xc2mat(i,j) = (     n*sumSAdotB - sumSA.*sumB) ./ ...
%                       sqrt( n*sumSAsqrd - sumSA.^2   ) ./ ...
%                       sqrt( n* sumBsqrd - sumB .^2   )  ;
