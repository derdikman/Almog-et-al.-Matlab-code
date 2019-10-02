function xc = normxcorr2ni(A,B)
A = interpnanrm(A);    
if nargin == 1
  B = A; 
else
  B = interpnanrm(B);
end

xc=normxcorr2(A,B);
end

%{
X = rand(25,4)*4-2;
V = X(:,1).^2 + X(:,2).^2 + X(:,3).^2 + X(:,4).^2;
%W = sum(X.^2,2); sum(V~=W)
%[X W]
X(5,1) = NaN; X(10,1) = NaN; X(12,2) = NaN; X(4,3) = NaN; X(1,4) = NaN; V(14) = NaN;
%nn=isnan(X);
%n=any(isnan(X),2);
%This code errors:
%F = scatteredInterpolant(X,V);
%Find the sample point indices of the NaN and then construct the interpolant:
nan_flags = isnan(X(:,1)) | isnan(X(:,2)) | isnan(X(:,3)) | isnan(X(:,4)) | isnan(V);

X(nan_flags,:) = [];
V(nan_flags) = [];
F = scatteredInterpolant(X,V);

drm = delaunay(rm(isnan()));
figure(1);
subplot(211);
imgsc(rm)
subplot(212);
imgsc(drm)
%}