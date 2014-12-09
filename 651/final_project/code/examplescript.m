load '~/files/afrl/matlab/data/sample.mat';
addpath(genpath('~/files/afrl/matlab/afrl_matlab'));
[simData, obsData] = chem_struct(sampleSims, sampleField, 2, 250);
nburn = 250;
nlev = 20;
nmcmc = 10000;
doPred = 1;
npvec = 1000;
sample_pout = chem_runcode(simData, obsData, 'nburn', nburn, ...
    'nlev', nlev, 'nmcmc', nmcmc, 'doPred', doPred, 'npvec', npvec);
chem_plots(sample_pout, 1, 'topdf', 1);
chem_plots(sample_pout, 2, 'subset', 1:5, 'topdf', 1);
chem_plots(sample_pout, 2, 'subset', 6:9, 'topdf', 0);
chem_plots(sample_pout, 3, 'topdf', 1);
chem_plots(sample_pout, 4, 'topdf', 1);
chem_plots(sample_pout, 5, 'topdf', 1);
chem_plots(sample_pout, 6, 'topdf', 1);
accept = get_accept(sample_pout.pvals);
mean(accept)
[min(mean(accept)) max(mean(accept))]
v = zeros(1001, 27);
for jj = 1:3
    for ii = 1:9
        a = sample_pout.obsData(ii).orig.y(jj);
        b = sample_pout.pred.zeta(jj, :, ii);
        b = log(b./(1-b));
        v(:,ii + 9*(jj-1)) = [a b];
    end
end
fileID = fopen('pred_v.txt','w');
fprintf(fileID,'%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f \n', v');
fclose(fileID);
