function pout = chem_runcode(simData, obsData, varargin)
pout = [];
nburn = 250;
nlev = 20;
nmcmc = 20000;
npvec = nmcmc;
doMCMC = 1;
doPred = 1;
parseAssignVarargs({'pout', 'nburn', 'nlev', 'nmcmc', 'npvec', ...
                    'doMCMC', 'doPred'});
if npvec > nmcmc
    npvec = nmcmc;
end
if doMCMC
    if isempty(pout)
        pout = setupModel(obsData, simData);
        % do some burn in get to make the posterior draws better
        if nburn > 0 && nlev > 0
            pout = stepsize(pout, nburn, nlev);
        end
    end
    pout = gpmmcmc(pout, nmcmc, 'step', 1);
    nmcmc = nmcmc+nburn*nlev;
    pout.pvec = floor(linspace(nburn*nlev+1, nmcmc, npvec));
end
if doPred
    % xpred is the vector of (scaled) x locations where predictions
    % are made, i.e. at what value of initial electron concentration
    % should predictions be made at? 
    % By default xpred is a column vector of the scaled x's
    xpred = pout.data.x;
    % for a vector of 50 equally spaced points between 0 and 1.
    % xpred = linspace(0, 1, 50)';
    npvec = length(pout.pvec)
    pout.pred = gPredict(xpred, pout.pvals(pout.pvec), pout.model, ...
        pout.data);
    % Put the predictions back on the original scale (i.e. [0, 1])
    n = size(pout.simData.orig.y, 1);
    if ~isfield(pout.pvals, 'theta')
        % only simulations were used
        zeta = (pout.simData.Ksim * pout.pred.w');
        for ii = 1:n
            zeta(ii, :) = zeta(ii, :) * pout.simData.orig.ysd(ii) + ...
                pout.simData.orig.ymean(ii);
        end
        zeta = 1./(1+exp(-zeta));
    else 
        % both field and simlulation
        % initalize 3-dimensional matrices
        eta = zeros(size(pout.simData.Ksim, 1), npvec, ...
            length(xpred));
        delta = zeros(size(eta));
        zeta = zeros(size(eta));
        for ii = 1:length(xpred)
            % back-transform the predictions
            eta(:, :, ii) = pout.simData.Ksim * pout.pred.u(:, ...
                :, ii)' .* repmat(pout.simData.orig.ysd, 1, npvec);
            delta(:, :, ii) = pout.simData.orig.Dsim * ...
                pout.pred.v(:, :, ii)' .* repmat( ...
                pout.simData.orig.ysd, 1, npvec);
            % sum the modeled simulations and discrepancy
            zeta(:, :, ii) = eta(:, :, ii) + delta(:, :, ii) + ...
                repmat(pout.simData.orig.ymean, 1, npvec);
            % put into [0, 1]
            eta(:, :, ii)= 1./(1+exp(-(eta(:, :, ii) + repmat( ...
                pout.simData.orig.ymean, 1, npvec))));
            zeta(:, :, ii) = 1./(1+exp(-zeta(:, :, ii)));
            % impose sum-to-one constraint
            eta(:, :, ii) = eta(:, :, ii) ./ repmat(sum(eta( ...
                :, :, ii)), n, 1);
            zeta(:, :, ii) = zeta(:, :, ii) ./ repmat(sum(zeta( ...
                :, :, ii)), n, 1);
            delta(:, :, ii) = zeta(:, :, ii) - eta(:, :, ii);
        end
    end
    % summarize the predictions further by taking the median and
    % 2.5th and 97.5th quantiles. these should be used in plotting
    pout.pred.xpred = xpred;
    pout.pred.zeta = zeta;
    pout.pred.zmedian = squeeze(median(zeta, 2));
    pout.pred.zlower = squeeze(prctile(zeta, 2.5, 2));
    pout.pred.zupper = squeeze(prctile(zeta, 97.5, 2));
    pout.pred.dmedian = squeeze(median(delta, 2));
    pout.pred.dlower = squeeze(prctile(delta, 2.5, 2));
    pout.pred.dupper = squeeze(prctile(delta, 97.5, 2));
    pout.pred.emedian = squeeze(median(eta, 2));
    pout.pred.elower = squeeze(prctile(eta, 2.5, 2));
    pout.pred.eupper = squeeze(prctile(eta, 97.5, 2));
end
end
