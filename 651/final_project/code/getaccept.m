function accepts = get_accept(pvals)
niter = length(pvals);
f = fieldnames(pvals);
keep = [];
for ii = 1:length(f)
    temp = char(f(ii));
    if strcmp(temp(end-2:end), 'Acc')
        keep = [keep ii];
    end
end
f = f(keep);
nf = length(f);
nparams = [];
for ii = 1:length(f)
    nparams(ii) = length(pvals(1).(char(f(ii))));
end
cp = [0 cumsum(nparams)];
accepts = zeros(niter, sum(nparams));
for ii = 1:niter
    ll = 0;
    for jj = 1:nf
        for kk = 1:nparams(jj)
            ll = ll + 1;
            accepts(ii, ll) = pvals(ii).(char(f(jj)))(kk);
        end
    end
end
end
