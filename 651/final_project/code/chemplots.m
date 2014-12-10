function chem_plots(pout, plotnum, varargin)
model = pout.model;
pu = model.pu;
p = model.p;
q = model.q;
data = pout.data;
etaMod = 0;
n_exp = length(pout.obsData);
if ~isfield(pout.pvals, 'theta')  % eta-only model: simulation data; no
    etaMod = 1;                   % observed data used in the analysis
end
if etaMod
    % eta-only model
    subset = 1:p;
    labs = cell(1, p);
    for ii = 2:p
        labs{ii} = int2str(ii - 1);
    end
    labs{1} = 'D';
else
    % observed data included
    subset = 1:q;
    labs = cell(1, q);
    for ii = 1:q
        labs{ii} = int2str(ii);
    end
end
ngrid = 21;
ctr = 0;    % center
topdf = 1;    % save figure as a .pdf
filename = ['Figure' int2str(plotnum) '.pdf'];
pvec = pout.pvec;
perc = [0.5 0.9];
parseAssignVarargs({'ngrid', 'subset', 'ctr', 'topdf', 'filename', ... 
                    'pvec', 'perc'});
if isfield(pout, 'pred')
    pred = pout.pred;
else
    pred = [];
end
pvals = pout.pvals(pvec);
doPlot = zeros(5, 1);
if exist('plotnum', 'var');
    doPlot(plotnum) = 1;
end
lsub = length(subset);
grid = linspace(0, 1, ngrid);
if doPlot(1)
    figure(1);
     clf;
    colormap([0 0 0]);
    bu = [pvals.betaU]';
    ru = exp(-bu/4);
    for ii = 1:pu
        if etaMod
            r = ru(:, (ii-1)*p+1:ii*p);
        else
            r = ru(:, (ii-1)*(p+q)+p+1:ii*(p+q));
        end
        r = r(:, subset);
        sep = 0.3;    % Needs to be changed according to value of pu
        gPackSubplot(pu, 1, ii, 1, sep);
        boxplot(r);
        ylab = ['\rho_{w', num2str(ii), 'k}'];
        ylabel(ylab, 'FontSize', 18);
        set(gca,'XTickLabel',{' '})    % removes the ticks
        if ii == pu
            ylim([0 1.01]);
            xlabel('Simulator Inputs', 'FontSize', 18);
        else
             ylim([0.01 1.00]);
        end
        if ii == 1
            text(1:length(labs(subset)), ones([1 length(labs(subset ...
                ))])+0.08, labs(subset), 'FontSize', 18, ...
                'HorizontalAlignment', 'center');
        end
    end
    figure(1);
end
if doPlot(2)
    if length(subset) > 12
        display('NOTE: With a large number of inputs, plotting all two-way');
        display('      distributions make take a bit of time. It is possible');
        display('      to look at a selection of inputs by specifying the');
        display('      ''subset'' variable in the argument line to a smaller ');
        display('      set of numbers. Doing so will also make the plots');
        display('      easier to see.');
    end
    figure(2);
    clf;
    if length(pvec) > 1000
        pvec2 = pvec(floor(linspace(1, length(pvec), 1000)));
    else
        pvec2 = pvec;
    end
    if etaMod
        t = [pout.pvals(pvec2).betaU]';    % These values aren't in [0,1]
                                        % like theta is.
         % temporary fix (scales betaU to be in [0,1])
        for ii = subset
             t(:, ii) = transform(t(:, ii),'col',[0 1]);
         end
        % end temp fix
    else
        t = [pout.pvals(pvec2).theta]';
    end
    t = t(:, subset);
    gPlotMatrix(t, 'Pcontours', perc, 'shade', 1, 'lstyle', 'imcont', ...
        'ngrid', ngrid, 'ksd', 0.1, 'labels', labs(subset));
    figure(2);
    % PROBLEMS:
    % Getting issues where some contours are not showing up. (Triflate1,
    % inputs 16 and 17, for example).
    % Possibly only on eta-only models, the plot matrix isn't symmetric.
    % Maybe because of the dummy variable?
    % Find out what the arguments to gPlotMatrix do.
    % See if the diagonals have the right range.
end
if doPlot(3)
    if ~isfield(pout, 'pred')
        error(['No predictions found. Use ''chem_runcode''.']);
    else 
        figure(3)
        clf
        sep = 0.50;
        % simulator
        gPackSubplot(1, 3, 1, 1, sep, 0);
        hold on;
        plot(pout.pred.xpred, pout.pred.emedian, 'LineWidth', 3)
        plot(pout.pred.xpred, pout.pred.elower, '--')
        plot(pout.pred.xpred, pout.pred.eupper, '--')
        for ii = 1:n_exp
            plot(pout.obsData(ii).x, 1./(1+exp(-pout.obsData(ii). ...
                orig.y)), 'ko', 'MarkerSize', 4, 'LineWidth', 3, ...
                'MarkerFaceColor', 'k')
        end
        ylabel('Rel. Abund.', 'FontSize', 18);
        title('Simulator', 'FontSize', 18);
        % discrepancy
        gPackSubplot(1, 3, 1, 2, sep, 0);
        hold on;
        plot(pout.pred.xpred, pout.pred.dmedian, 'LineWidth', 3)
        plot(pout.pred.xpred, pout.pred.dlower, '--')
        plot(pout.pred.xpred, pout.pred.dupper, '--')
        plot(pout.pred.xpred, zeros( 1, length(pout.pred.xpred)), 'k-', 'LineWidth', 2)
        title('Discrepancy', 'FontSize', 18);
        xlabel('scaled initial electron', 'FontSize', 18);
        % discprency-adjusted
        gPackSubplot(1, 3, 1, 3, sep, 0);
        hold on;
        plot(pout.pred.xpred, pout.pred.zmedian, 'LineWidth', 3)
        plot(pout.pred.xpred, pout.pred.zlower, '--')
        plot(pout.pred.xpred, pout.pred.zupper, '--')
        for ii = 1:n_exp
            plot(pout.obsData(ii).x, 1./(1+exp(-pout.obsData(ii). ...
                orig.y)), 'ko', 'MarkerSize', 4, 'LineWidth', 3, ...
                'MarkerFaceColor', 'k')
        end
        title('Calibrated Predictions', 'FontSize', 18);
        figure(3)
    end
end
if doPlot(4)
    if ~isfield(pout, 'pred')
        error(['No predictions found. Use ''chem_runcode''.']);
    else 
        figure(4)
        clf
        cols = get(0, 'DefaultAxesColorOrder');
        n = size(pout.simData.yStd, 1);
        sep = 0.75;
        top = 0.0;
        for jj = 1:n;
            if (jj > 1)
                top = 0.5;
            end
            % simulator
            gPackSubplot(n, 3, jj, 1, sep, top);
            hold on;
            plot(pout.pred.xpred, pout.pred.emedian(jj, :), 'Color', cols(jj, :), 'LineWidth', 3)
            plot(pout.pred.xpred, pout.pred.elower(jj, :), '--', 'Color', cols(jj, :))
            plot(pout.pred.xpred, pout.pred.eupper(jj, :), '--', 'Color', cols(jj, :))
            for ii = 1:n_exp
                plot(pout.obsData(ii).x, 1./(1+exp(-pout.obsData(ii). ...
                    orig.y(jj, :))), 'ko', 'MarkerSize', 4, 'LineWidth', 3, ...
                    'MarkerFaceColor', 'k')
            end
            ylabel('Rel. Abund.', 'FontSize', 18);
            if (jj == 1)
                title('Simulator', 'FontSize', 18);
            end
            % discrepancy
            gPackSubplot(n, 3, jj, 2, sep, top);
            hold on;
            plot(pout.pred.xpred, pout.pred.dmedian(jj, :), 'LineWidth', 3, 'Color', cols(jj, :))
            plot(pout.pred.xpred, pout.pred.dlower(jj, :), '--', 'Color', cols(jj, :))
            plot(pout.pred.xpred, pout.pred.dupper(jj, :), '--', 'Color', cols(jj, :))
            plot(pout.pred.xpred, zeros( 1, length(pout.pred.xpred)), 'k-', 'LineWidth', 2)
            if jj == 1
                title('Discrepancy', 'FontSize', 18);
            end
            if jj == n
                xlabel('scaled initial electron', 'FontSize', 18);
            end
            % discprency-adjusted
            gPackSubplot(n, 3, jj, 3, sep, top);
            hold on;
            plot(pout.pred.xpred, pout.pred.zmedian(jj, :), 'LineWidth', 3, 'Color', cols(jj, :))
            plot(pout.pred.xpred, pout.pred.zlower(jj, :), '--', 'Color', cols(jj, :))
            plot(pout.pred.xpred, pout.pred.zupper(jj, :), '--', 'Color', cols(jj, :))
            for ii = 1:n_exp
                plot(pout.obsData(ii).x, 1./(1+exp(-pout.obsData(ii). ...
                    orig.y(jj, :))), 'ko', 'MarkerSize', 4, 'LineWidth', 3, ...
                    'MarkerFaceColor', 'k')
            end
            if jj == 1
                title('Calibrated Predictions', 'FontSize', 18);
            end
        end
        figure(4)
    end
end
if doPlot(5)
    % plotting the data
    n = size(pout.simData.yStd, 1);
    cols = get(0, 'DefaultAxesColorOrder');
    x = zeros(n_exp, 1);
    y = zeros(n_exp, n);
    for ii = 1:n_exp
        x(ii) = pout.obsData(ii).x;
        y(ii, :) = 1./(1+exp(-pout.obsData(ii).orig.y));
    end
    figure(5)
    clf
    hold on;
    n_sims = size(pout.simData.orig.y, 2) / n_exp;
    for ii = 1:n
        for jj = 1:floor(n_sims/2)
            plot(x, 1./(1+exp(-pout.simData.orig.y(ii, (1:n_exp) + n_exp*(jj-1)))), ...
                'Color', min(1, (cols(ii, :)+0.7)/1.0), 'LineWidth', 1)
        end
    end
    for ii = 1:n
        plot(x, y(:, ii), 'Color', cols(ii, :), 'LineWidth', 4)
        plot(x, y(:, ii), 'ko', 'MarkerSize', 4, 'LineWidth', 4)
    end
    title('Field and Simulation Data', 'FontSize', 18);
    xlabel('scaled initial electron', 'FontSize', 18);
    ylabel('Rel. Abund.', 'FontSize', 18);
    figure(5)
end
if doPlot(6)
    showPvals(pout.pvals(pout.pvec));
end
if topdf
    if doPlot(1)
        if pu > 1
            margs = [-0.25-0.015*lsub -1.75-0.018*lsub -0.75 -0.25];    % [Left Right Top Bottom]
        else
            margs = [-0.25-0.015*lsub -0.75-0.018*lsub -0.25 -0.25];
        end
        paper = [6.5+0.1*lsub min(4*pu+0.25, 11)];    % [Width Height]
                                    % Height should changed based on pu
    end
    if doPlot(2)
        margs = [-0.75 -0.75 -0.75 -0.75];    % [Left Right Top Bottom]
        paper = [9 9];    % [Width Height]
    end
    if doPlot(3)
        if etaMod
            margs = [-0.25 -0.75 -0.5 -0.25];    % [Left Right Top Bottom]
            paper = [10 6];    % [Width Height]
        else
            margs = [-1.25 -3.0 0 0.25];    % [Left Right Top Bottom]
            paper = [18 6];    % [Width Height]
        end
    end
    if doPlot(4)
        margs = [-0.25 -1.5 0 0.25];    % [Left Right Top Bottom]
        paper = [12 3.5*n];
    end
    if doPlot(5)
        margs = [-0.5 -0.5 0 0.25];    % [Left Right Top Bottom]
        paper = [9 9];
    end
    if doPlot(6)
        margs = [-0.0 -0.5 -0.5 -0.25];    % [Left Right Top Bottom]
        paper = [9 9];
    end
    fp = fillPage(gcf, 'margins', margs, 'papersize', paper);
    if ~strcmp(filename(length(filename)-3:length(filename)), '.pdf')
        filename = [filename '.pdf'];
    end
    set(gcf, 'InvertHardCopy', 'off');
    print(gcf, '-dpdf', '-r300', filename)
    set(gcf, fp);
end
end
