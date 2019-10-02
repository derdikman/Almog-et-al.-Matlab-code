function c = driftccof(a,b)
nbins = 50; %px = px+max(px); py = py+max(py);

px = a(:,1);pxi = discretize(px, min(px):range(px)/nbins:max(px));
py = a(:,2);pyi = discretize(py, min(py):range(py)/nbins:max(py));
rma = accumarray([pyi, pxi], 1, [nbins nbins]);

px = b(:,1);pxi = discretize(px, min(px):range(px)/nbins:max(px));
py = b(:,2);pyi = discretize(py, min(py):range(py)/nbins:max(py));
rmb = accumarray([pyi, pxi], 1, [nbins nbins]);

c = ccof(rma,rmb);

end