%scales nan-containing dot product values by difference in number of dot product
%elelements compared to [max] number of elements in non-nan cross product. 
%A,B may contain nans.
%mB  optional, can be calculated outside function for optimizing, see how below
function xc2nansscaled = xcorr2n(A,B,mB)%(mat1,mat2)
if nargin == 1
    B=A;
end 
Anel = ones(size(A)); 
Bnel = ones(size(B));
%nDotProdsAll = xcorr2(Anel,Bnel); sizeo=size(nDotProdsAll);
sizeo = [ size(A,1)+size(B,1)-1, size(A,2)+size(B,2)-1 ];
if nargin == 3 %could pass this in to speed up;   
    m = mB;
else 
    nDotProdsAll = xcorr2(Anel,Bnel); %simple 2D cross correlation, dot product of sliding images.
    m = max(nDotProdsAll(:)); %length of dot product vector at maximum overlap (center)
end
nDotProdsAll = ones(sizeo)*m; %max number of dot products used in calculating non-nan xcorr (center).
%another option (below) which is closer in result to what xcorr2 does, which treats values
 %beyond image as 0's, is to just use nDotProdsAll calculated abve.
Anel(isnan(A)) = 0;
Bnel(isnan(B)) = 0;
nDotProdsNonNan = xcorr2(Anel,Bnel);% number of dot products used when removing nans in xcorr.
%scale nans, ignore values with a lot of nans to by setting to nan [or 0]
scalenans = nDotProdsAll./nDotProdsNonNan;
%scalenans = abs(scalenans);
%scalenans = scalenans-min(scalenans(:))+0.1;
%scalenans=scalenans.^(-2/1);
%stats(scalenans(:));
l=100;
scalenans(scalenans>l) = nan;%nan or 0? %ignore values where number of nans to non nans is > 100:1 
scalenans(scalenans<1/l) = nan;%nan or 0? %ignore values where number of nans to non nans is > 100:1 
%scalenans=scalenans./max(scalenans(:));
%scalenans(scalenans<1e-5) = 0; 
%xcorr setting nans to 0
A(isnan(A)) = 0;
B(isnan(B)) = 0;
xc0 = xcorr2(A,B);%normxcorr2(A,B);
%nans scaled
xc2nansscaled = xc0 .* scalenans;

%figure(1);imgsc(xc2nansscaled,2);

end

%set set vals with a lot of nans to 0
%scalenans0=scalenans;
%scalenans0(scalenans>100) = 0;
%xc2nansscaled=xc0.*scalenans0;