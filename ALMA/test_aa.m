%% Using https://github.com/kreutz-lab/Rcall
addpath('~/Dropbox/projects/sandbox/Rcall/')
Rinit
load('../test_A.mat')

%%

n = size(A, 2);
L = size(A, 1);

params.K = [2, 3, 5];
params.L = L;
params.M = 3;
params.n = n;

%%
z = zeros(L);
xi = zeros(params.M, n);

[Q1, W1] = AltMin(A, params);
labelhat = kmeans(W1, params.M, "Replicates", 100);

%%
Rpush('ztru', cast(z_tru,'double'), 'zhat', cast(labelhat,'double'))
Rrun('nmi <- nett::compute_mutual_info(ztru, zhat)')
Rpull('nmi')