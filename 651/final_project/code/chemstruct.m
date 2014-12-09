function [simData, obsData] = chem_struct(simulationData, ...
    fieldData, pu, subset)
if ~exist('pu','var')
    pu = 1; % the number of principal components
end
if ~exist('subset', 'var')
    subset = 1:length(simulationData);
else
    if isscalar(subset)
        % subset is the index for the observations in
        % the simulations that will be used
        % subset(1) must equal 1 since the first row
        % in the simulations variable is just a header
        % (but it gets removed)
        subset = abs(subset);
        if subset > 1   % is the number of sims to use
            subset = floor(subset);
            subset = randperm(length(simulationData), subset);
            subset(1) = 1;
            subset = sort(subset);
        else    % is a proportion
            subset = randperm(length(simulationData),  ceil(length( ...
                simulationData)*subset)+1);
            subset(1) = 1;
            subset = sort(subset);
        end
    end
    subset(1) = 1;
end
sims = simulationData(subset, :);
Start = Inf;
for col = 1:size(sims, 2)
    if all(sims(2:end, col) == sims(2:end, end))
        design = sims(2:end, 1:col-1);
        Start = col + 2;
    end
    if col > Start
        if (all(isnan(sims(:, col))) ...
            || all(sims(:, col) == sims(1, col)) ...
            || col == size(sims, 2))
            n_exp = col - Start;
            n_out = floor ((size(sims, 2) - Start) / n_exp);
            break
        end
    end
end
for ii = 1:n_out
    [~, order] = sort(sims(1, Start+n_exp*(ii-1)+ii-1:Start+ ...
        n_exp*ii+ii-2));
    sims(:, Start+n_exp*(ii-1)+ii-1:Start+n_exp*ii+ii-2) = sims(:, ...
        order+Start+n_exp*(ii-1)+ii-2);
end
xsim = repmat(sims(1, Start:(Start+n_exp-1))', size(sims, 1)-1, 1);
design = [xsim design(ceil([1:n_exp*size(design,1)]/n_exp),:)];
m = size(design, 1); 
xmin = min(design);
xrange = range(design);
design = transform(design, 'col', [0 1]);
ysim = zeros(n_out, m);
for yy = 1:(m/n_exp)
    for zz = 1:n_exp
        ysim(:, (yy-1)*n_exp+zz) = sims(yy+1, Start + zz-1:(n_exp+ ...
            1):Start + n_out*n_exp + zz);
    end
end
ysim = log(ysim./(1-ysim));
simData = struct('x', [], 'yStd', [], 'Ksim', [], 'orig', ...
    struct('y', [], 'ymean', [], 'ysd', [], 'Dsim', [], 'xmin', [], ...
    'xrange', []));
ysimmean = mean(ysim, 2); % the mean of each row
ysimStd = ysim - repmat(ysimmean, 1, m); % subtract each row by mean
ysimsd = std(ysimStd, 0, 2); % get standard deviated for each row
ysimStd = ysimStd ./ repmat(ysimsd, 1, m); % dive rows by s.d.
[U, S, ~] = svd(ysimStd, 0);
Ksim = U(:, 1:pu)*S(1:pu, 1:pu)./sqrt(m);
PCAvar = sum(sum(S(1:pu, 1:pu))) / sum(sum(S));
pv = n_out;
Dsim = eye(pv);
if isempty(fieldData)
    obsData = [];
else
    obs = fieldData;
    % field output is transformed to put on (-Inf, +Inf)
    obs(:, 2:end) = log(obs(:, 2:end)./(1-obs(:, 2:end)));
    obsBuild = struct('x', [], 'yStd', [], 'Kobs', [], 'Dobs', [], ...
        'Sigy', [], 'orig', struct('y', [], 'ymean', [], 'sd', []));
    obsData = repmat(obsBuild, 1, n_exp);
    for ii = 1:n_exp
        obsData(ii).x = (obs(ii, 1)-xmin(1))/xrange(1);
        % this assumes known observation errors, probably
        % not the best idea
        obsData(ii).Sigy = diag((0.1^2) * ones(n_out, 1));
        obsData(ii).orig.y = obs(ii, 2:end)';
        obsData(ii).orig.ymean = ysimmean;
        obsData(ii).orig.sd = sqrt(diag(obsData(ii).Sigy));
        obsData(ii).yStd = (obsData(ii).orig.y - ysimmean) ./ ysimsd;
        obsData(ii).Kobs = Ksim;
        obsData(ii).Dobs = eye(n_out);
    end
end
simData.x = design;
simData.yStd = ysimStd;
simData.Ksim = Ksim;
simData.orig.y = ysim;
simData.orig.ymean = ysimmean;
simData.orig.ysd = ysimsd;
simData.orig.Dsim = Dsim;
simData.orig.xmin = xmin;
simData.orig.xrange = xrange;
simData.orig.subset = subset(2:end);
simData.orig.PCAvar = PCAvar;
end
