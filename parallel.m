n = 50;
A = 500;

%serial
tic
s = zeros(n);
for i = 1:n
    s(i) = max(abs(eig(rand(A))));
end
toc
%parallel 
tic
ticBytes(gcp);
p = zeros(n);
parfor i = 1:n
    p(i) = max(abs(eig(rand(A))));
end
tocBytes(gcp)
toc