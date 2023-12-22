% how to parallize the loop in Mtalab 
%   https://de.mathworks.com/help/parallel-computing/parfor.html
%%
tic
n = 200;
A = 500;
a = zeros(1,n);
parfor i = 1:n
    a(i) = max(abs(eig(rand(A))));
end
toc
%%%
%%

cluster = parcluster;
values = [3 3 3 7 3 3 3];
parfor (i=1:numel(values),cluster)
    out(i) = norm(pinv(rand(values(i)*1e3)));
end



