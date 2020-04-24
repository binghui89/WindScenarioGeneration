function run_WIND_toolkit_9sites(red_flag, season, method_sr)
% red_flag: true or false, true is the full case (45 sites), 
% false is the reduced case (20 sites). Defaults to true.
% season indicates the target date for scenario generation
% 1: Feb. 20, 2012, 2: May 10, 2012, 3: Aug. 10, 2012, 4: Nov. 10, 2012
% Defaults to 2.
% method_sr switch between the sum-based reduction (1) and site-based
% reduction (2). Defauts to 2.
addpath('./packages/copula-matlab');

if nargin==0
    method_sr = 2;  % 1: sum-based reduction; 2: site-based reduction
    season = 2;
    red_flag = true; 
elseif nargin==1
    season = 2;
    method_sr = 2;  % 1: sum-based reduction; 2: site-based reduction
elseif nargin==2
    method_sr = 2;  % 1: sum-based reduction; 2: site-based reduction
end


nS_target = 10; % Number of scenarios we'd like to reduce to.
% nT1: % Number of hours of the training dataset
% nT2: % Number of hours of the probabilistic forecast
if season == 1
    nT1 = 44808 - 24; % Feb. 20, 2012
    nT2 = 24;
elseif season == 2
    nT1 = 46968 - 24; % May 10, 2012
    nT2 = 24;
elseif season == 3
    nT1 = 49176 - 24; % Aug. 10, 2012
    nT2 = 24;
elseif season == 4
    nT1 = 51384-24; % Nov. 10, 2012
    nT2 = 24;
end
nS1 = 8000; % Number of burn-in scenarios
nS2 = 1000; % Number of generated scenarios

data = read_all();

% This is randomly selected 20 sites
if red_flag
    selected_sites = randsample(data.sid, 20);
else
    selected_sites = data.sid;
end

I_selected = zeros(data.nI, 1);
for i = 1: length(selected_sites)
    s = selected_sites(i);
    I_selected(data.sid==s) = 1;
end
I_selected = logical(I_selected);

data_s = data;
data_s.nI = sum(I_selected);
data_s.nT = data.nT;
data_s.cap = data.cap(I_selected);
data_s.va = data.va(:, I_selected);
data_s.xa = data.xa(:, I_selected);
data_s.xf = data.xf(:, I_selected);
data_s.xa_MW = data.xa_MW(:, I_selected);
data_s.xf_MW = data.xf_MW(:, I_selected);
data_s.sid = data.sid(I_selected);

%% First, let's do cluster analysis
Y = 2007: 2012;
Karray = 1: 1: 20; % # of clusters
f_cluster = 'kmeans'; % We can also test kmedoids
s_cluster   = nan(length(Karray), length(Y)); % Silhouette metric
d_cluster   = nan(length(Karray), length(Y)); % Dunn metric
idxC = nan(data_s.nI, length(Karray), size(Y, 1));
for j = 1: length(Y)
    y = Y(j);
    x = data_s.va(data_s.yyyy==y, :); % Let's try wind speed
    idxC(:, :, j) = cluster_analysis_tmp(x,f_cluster,Karray);
    x_t = x';
    for i = 1: length(Karray)
        [s, ~] = silhouette(x_t, idxC(:, i, j), 'correlation');
        s_cluster(i, j) = mean(s);
    end
end

% Add K = 1
if red_flag
    Karray = [1 Karray];
    idxC = cat(2, true(data_s.nI, 1, length(Y)), idxC);
end

%% Start generating scenarios
xnew_all = cell(length(Karray), 1);
time_taken_sc = nan(length(Karray), 1); % Scenario generation time

mse_allK = cell(length(Karray), 1);
mae_allK = cell(length(Karray), 1);
for nK = Karray
    tic;
    idx = idxC(:, Karray==nK, 1); % Just use cluster results from the first year
    mse = nan(nK, 1);
    mae = nan(nK, 1);
    xnewK = cell(nK, 1);
    for k = 1: nK
        i_cluster = (idx==k);
    %     xa = data_s.xa(:, i_cluster);
    %     xf = data_s.xf(:, i_cluster);
    %     sid = data_s.sid(i_cluster);
        xnew = gibbs2(data_s, nT1, nT2, nS1, nS2, i_cluster);
        xnewK{k} = xnew;
        abserr = abs(repmat(data_s.xa(nT1+1:nT1+nT2, i_cluster), 1, 1, nS2) - xnew(:, :, nS1+1: nS1+nS2));
%         abserr = abs(data_s.xa(nT1+1:nT1+nT2, i_cluster) - xnew);
        sqrerr = abserr.^2;
        mae(k) = mean(abserr(:));
        mse(k) = sqrt(mean(sqrerr(:)));
    end
    xnew_all{Karray==nK} = xnewK;
    mae_allK{Karray==nK} = mae;
    mse_allK{Karray==nK} = mse;
    time_taken_sc(Karray==nK) = toc;
end

mse_mean = nan(length(Karray), 1);
mae_mean = nan(length(Karray), 1);
for i = 1: length(Karray)
    mse_mean(i) = mean(mse_allK{i});
    mae_mean(i) = mean(mae_allK{i});
end

% Calculate MW
xnewMW_all = cell(size(xnew_all));
for nK = Karray
    idx = idxC(:, Karray==nK, 1); % Just use cluster results from the first year
    xnewMWK = cell(nK, 1);
    for k = 1:nK
        sid_cluster = data_s.sid((idx==k));
        xnew = xnew_all{Karray==nK}{k};
        xnewMW = nan(nT2, length(sid_cluster), nS2);
        for i = 1: length(sid_cluster)
            s = sid_cluster(i);
            xnewMW(:, i, :) = xnew(:, i, nS1+1: nS1+nS2).*data_s.cap(data_s.sid==s);
        end
        xnewMWK{k} = xnewMW;
    end
    xnewMW_all{Karray==nK} = xnewMWK;
end

%% Reduce scenarios
xnewMWred_all = cell(size(xnew_all));
time_taken_red = zeros(max(Karray), length(Karray));
array_prob_red = nan(nS_target, length(Karray)); % Probabilities
for nK = Karray
    % This is a cell of arrays, each array represents scenarios from a 
    % cluster. This cell includes nK clusters.
    xnewMWK = xnewMW_all{Karray==nK};

    xnewMWredK = cell(size(xnewMWK)); % This cell contains all reduced scenarios.
    for k = 1: nK
        t_start = tic;
        xnewMW = xnewMWK{k};
        if k == 1
            scenario_index = 1:nS2;
            p = 1/nS2*ones(nS2, 1);
            switch method_sr
                case 1 % Sum-based reduction
                    x_in = squeeze(sum(xnewMW, 2));
                case 2 % Site-based reduction
                    tmp_size = size(xnewMW);
                    x_in = reshape(xnewMW, tmp_size(1)*tmp_size(2), tmp_size(3));
                otherwise
                    warning('Unexpected method index.');
            end
            [x_out, p_red, b_selected] = reduction_forward(x_in, p(:), nS_target);
            scenario_index_red = scenario_index(b_selected);

        else
            scenario_index_dn = 1:nS2;
            scenario_index = combvec(scenario_index_red, scenario_index_dn);
            p_dn = 1/nS2*ones(nS2, 1);
            p = prod(combvec(p_red(:)', p_dn(:)'), 1);
            
            switch method_sr
                case 1
                    x_in_dn = squeeze(sum(xnewMW, 2));
                case 2
                    tmp_size = size(xnewMW);
                    x_in_dn = reshape(xnewMW, tmp_size(1)*tmp_size(2), tmp_size(3));
                otherwise
                    warning('Unexpected method index.');
            end
            x_in = combvec(x_out, x_in_dn);
            [x_out, p_red, b_selected] = reduction_forward(x_in, p(:), nS_target);
            scenario_index_red = scenario_index(:, b_selected);
        end
        time_taken_red(k, Karray==nK) = toc(t_start); % Record time for each cluster
    end
    for k = 1: nK
        xnewMW = xnewMWK{k};
        xnewMWredK{k} = xnewMW(:, :, scenario_index_red(k, :));
    end
    xnewMWred_all{Karray==nK} = xnewMWredK;
    array_prob_red(:, Karray==nK) = p_red(:);
end

%% Plot results, all years of cluster results
% figure();
% if length(mse_mean) == length(Karray)
%     plotyy(Karray, s_cluster, Karray, mse_mean);
% elseif length(mse_mean) == length(Karray) + 1
%     plotyy(Karray, s_cluster, [1 Karray], mse_mean);
% end
% figure();
% if length(mse_mean) == length(Karray)
%     plotyy(Karray, s_cluster, Karray, mae_mean);
% elseif length(mse_mean) == length(Karray) + 1
%     plotyy(Karray, s_cluster, [1 Karray], mae_mean);
% end

%% Plot results, only the first year cluster results
fig = figure('Position', [100 100 560 420*0.8]);
% if length(mse_mean) == length(Karray)
%     plotyy(Karray, s_cluster(:, 1), Karray, mse_mean);
% elseif length(mse_mean) == length(Karray) + 1
%     plotyy(Karray, s_cluster(:, 1), [1 Karray], mse_mean);
% end
% figure();
% [AX,H1,H2] = plotyy(Karray, s_cluster(:, 1), [1 Karray], mse_mean);
ax = subplot(1, 1, 1);
yyaxis(ax, 'left'); 
if Karray(1) == 1 
    plot(Karray(2:end), s_cluster(:, 1), 'Marker', 's', 'LineWidth', 2, 'Color', 'k');
else
    plot(Karray(1:end), s_cluster(:, 1), 'Marker', 's', 'LineWidth', 2, 'Color', 'k');
end
ylabel('Silhouette coefficient') % left y-axis
yyaxis(ax, 'right');
plot(Karray, mse_mean, 'LineStyle', '--', 'Marker', 's', 'LineWidth', 2,'Color', 'k');
ylabel('RMSE');
xlabel('Number of clusters (K)');
set(findall(gcf,'-property','FontSize'),'FontSize',18);
legend('Silhouette', 'RMSE'); legend('boxoff');
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%% Plot results, all year mean cluster results
% figure();
% if length(mse_mean) == length(Karray)
%     plotyy(Karray, mean(s_cluster, 2), Karray, mse_mean);
% elseif length(mse_mean) == length(Karray) + 1
%     plotyy(Karray, mean(s_cluster, 2), [1 Karray], mse_mean);
% end
% figure();
% if length(mse_mean) == length(Karray)
%     plotyy(Karray, mean(s_cluster, 2), Karray, mae_mean);
% elseif length(mse_mean) == length(Karray) + 1
%     plotyy(Karray, mean(s_cluster, 2), [1 Karray], mae_mean);
% end

%% Explore relationship between correlation coefficient and MSE/MAE
% mse_per_c  = nan((1+Karray(end))*Karray(end)/2, 1);
% mae_per_c  = nan((1+Karray(end))*Karray(end)/2, 1);
% coef_per_c = nan((1+Karray(end))*Karray(end)/2, 1);
% mse_max = nan(Karray(end), 1);
% mse_min = nan(Karray(end), 1);
% 
% i = 0;
% for nK = 1: Karray(end)
%     mse = nan(nK, 1);
%     mae = nan(nK, 1);
%     for k = 1: nK
%         i = i + 1;
%         if nK == 1
%             i_cluster = true(data_s.nI, 1);
%         else
%             i_cluster = (idxC(:, Karray==nK, 1) == k);
%         end
%         xnew = xnew_all{nK}{k};
%         abserr = abs(repmat(data_s.xa(nT1+1:nT1+nT2, i_cluster), 1, 1, nS2) - xnew_all{nK}{k}(:, :, nS1+1: nS1+nS2));
%         sqrerr = abserr.^2;
%         mae_per_c(i) = mean(abserr(:));
%         mse_per_c(i) = sqrt(mean(sqrerr(:)));
%         mse(k) = mse_per_c(i);
%         mae(k) = mae_per_c(i);
%         coef_matrix = corrcoef(data_s.va(:, i_cluster));
%         n = sum(i_cluster);
%         coef_matrix_upper = triu(coef_matrix);
%         coef_per_c(i) = sum(coef_matrix_upper(:))/(n*(n+1)/2);
%     end
%     mse_max(nK) = max(mse);
%     mse_min(nK) = min(mae);
% end

%% PICP analysis
% xnew_all{1} = {xnew_all{1}}; % Otherwise xnew_all{nK}{K} later will raise an error
p_norminal = 10:10:90;

Kmax = max(Karray);

array_average_p_actual = nan(length(Karray), 9);

for nK = Karray
    idx = idxC(:, Karray==nK, 1);
    array_p_actual = nan(nK, 9);
    array_nsite = nan(1, nK);
%     figure;
%     plot(0:10:100, 0:10:100, 'r');
%     hold on;
    for K = 1: nK
        i_cluster = (idx==K);
        xnew = xnew_all{Karray==nK}{K}; 
        p_actual = picp(xnew(:, :, nS1+1: nS1+nS2) , data_s.xa(nT1+1:nT1+nT2, i_cluster), p_norminal);
%         plot([0 all_p_norminal 100], [0 all_p_actual 100], 'k');
        
        array_p_actual(K, :) = p_actual;
        array_nsite(K) = sum(i_cluster);
        average_p_actual = (array_nsite./sum(array_nsite))*array_p_actual;
        array_average_p_actual(Karray==nK, :) = average_p_actual;
    end
end

cell_legends = cell(length(Karray), 1);
for nK = Karray
    cell_legends{Karray==nK} = int2str(nK);
end

% PICP figure 1, all PICP curves
figure();
plot([0 p_norminal 100]', [zeros(length(Karray), 1) array_average_p_actual 100.*ones(length(Karray), 1)]');
xlabel('Nominal coverage rate');
ylabel('Actual coverage rate');
legend(cell_legends);

% PICP figure 2, average delta PICP vs. nK
figure()
delta_actual_nominal = abs(...
    array_average_p_actual - repmat(10:10:90, size(array_average_p_actual, 1), 1)...
    );
plot(Karray, mean(delta_actual_nominal, 2));
xlabel('Number of clusters (K)');
ylabel('Average delta PICP');

% PICP figure 3, max and min average delta PICP vs. nK
[~, rowmax] = max(mean(delta_actual_nominal, 2));
[~, rowmin] = min(mean(delta_actual_nominal, 2));
figure();
plot([0 p_norminal 100], [0 array_average_p_actual(rowmin, :) 100], 'r');
hold on;
plot([0 p_norminal 100], [0 array_average_p_actual(rowmax, :) 100], 'b');
plot(0:10:100, 0:10:100, 'k');
xlabel('Nominal coverage rate');
ylabel('Actual coverage rate');
cell_legend = {strcat('K=', int2str(Karray(rowmin))), strcat('K=', int2str(Karray(rowmax)))};
legend(cell_legend);
% end

%% Sharpness, based on the Wan paper
p_norminal = 10:10:90;

array_average_interval_scores = nan(length(Karray), 9);

for nK = Karray
    idx = idxC(:, Karray==nK, 1);
    array_interval_scores = nan(nK, 9);
    array_nsite = nan(1, nK);
    for K = 1: nK
        i_cluster = (idx==K);
        xnew = xnew_all{Karray==nK}{K}; 
        interval_scores = interval_score(xnew(:, :, nS1+1: nS1+nS2) , data_s.xa(nT1+1:nT1+nT2, i_cluster), p_norminal);
%         plot([0 all_p_norminal 100], [0 all_p_actual 100], 'k');
        
        array_interval_scores(K, :) = interval_scores;
        array_nsite(K) = sum(i_cluster);
        average_p_actual = (array_nsite./sum(array_nsite))*array_interval_scores;
        array_average_interval_scores(Karray==nK, :) = average_p_actual;
    end
end

cell_legends = cell(length(Karray), 1);
for nK = Karray
    cell_legends{Karray==nK} = int2str(nK);
end

% Sharpness figure 1, all Sharpness curves
figure();
plot(p_norminal', -array_average_interval_scores');
xlabel('Nominal coverage rate');
ylabel('Interval scores');
legend(cell_legends);

% Sharpness figure 2, average sharpness vs. nK
figure()
plot(Karray(:), mean(-array_average_interval_scores, 2));
xlabel('Number of clusters (K)');
ylabel('Average interval scores');

% Sharpness figure 3, highest and lowest sharpness curves vs. nK
[~, rowmax] = max(mean(-array_average_interval_scores, 2));
[~, rowmin] = min(mean(-array_average_interval_scores, 2));
figure();
plot(p_norminal, -array_average_interval_scores(rowmin, :), 'r');
hold on;
plot(p_norminal, -array_average_interval_scores(rowmax, :), 'b');
xlabel('Nominal coverage rate');
ylabel('Interval scores');
cell_legend = {strcat('K=', int2str(Karray(rowmin))), strcat('K=', int2str(Karray(rowmax)))};
legend(cell_legend);

%% Interval sizes
p_norminal = 10:10:90;

array_average_interval_sizes = nan(length(Karray), 9);

for nK = Karray
    idx = idxC(:, Karray==nK, 1);
    array_interval_size = nan(nK, 9);
    array_nsite = nan(1, nK);
    for K = 1: nK
        i_cluster = (idx==K);
        xnew = xnew_all{Karray==nK}{K}; 
        itvl_size = interval_size(xnew(:, :, nS1+1: nS1+nS2) , data_s.xa(nT1+1:nT1+nT2, i_cluster), p_norminal);
%         plot([0 all_p_norminal 100], [0 all_p_actual 100], 'k');
        
        array_interval_size(K, :) = itvl_size;
        array_nsite(K) = sum(i_cluster);
        average_p_actual = (array_nsite./sum(array_nsite))*array_interval_size;
        array_average_interval_sizes(Karray==nK, :) = average_p_actual;
    end
end

cell_legends = cell(length(Karray), 1);
for nK = Karray
    cell_legends{Karray==nK} = int2str(nK);
end

% Interval size figure 1, all Sharpness curves
figure();
plot(p_norminal', array_average_interval_sizes');
xlabel('Nominal coverage rate');
ylabel('Interval size');
legend(cell_legends);

% Interval size figure 2, average sharpness vs. nK
figure()
plot(Karray(:), mean(array_average_interval_sizes, 2));
xlabel('Number of clusters (K)');
ylabel('Average interval size');

% Interval size figure 3, highest and lowest sharpness curves vs. nK
[~, rowmax] = max(mean(array_average_interval_sizes, 2));
[~, rowmin] = min(mean(array_average_interval_sizes, 2));
figure();
plot(p_norminal, array_average_interval_sizes(rowmin, :), 'r');
hold on;
plot(p_norminal, array_average_interval_sizes(rowmax, :), 'b');
xlabel('Nominal coverage rate');
ylabel('Interval size');
cell_legend = {strcat('K=', int2str(Karray(rowmin))), strcat('K=', int2str(Karray(rowmax)))};
legend(cell_legend);

%% Time analysis
figure();
bar(Karray, [time_taken_sc sum(time_taken_red, 1)';], 'stacked');
legend('Scenario generation' ,'Scenario reduction')
xlabel('Number of clusters (K)')
ylabel('Time (s)')

%% Write data to file
csv_sid = [];
csv_xa_MW = [];
csv_xf_MW = [];
csv_xs_MW = [];

k_optim = 10; % The optimal number of clusters, based on visual inspection of PICP
xnewMWredK = xnewMWred_all{Karray==k_optim};
idx = idxC(:, Karray==k_optim, 1); % Just use cluster results from the first year

for k = 1: k_optim
    i_cluster = (idx==k);
    csv_xa_MW = cat(2, csv_xa_MW, data_s.xa_MW(nT1+1:nT1+nT2, i_cluster));
    csv_xf_MW = cat(2, csv_xf_MW, data_s.xf_MW(nT1+1:nT1+nT2, i_cluster));
    csv_xs_MW = cat(2, csv_xs_MW, xnewMWredK{k});
    sid = data_s.sid(i_cluster);
    csv_sid = cat(1, csv_sid, sid(:));
end

header = {...
    'xa', 'xf', ...
    'N1', 'N2', 'N3', 'N4', 'N5', ...
    'N6', 'N7', 'N8', 'N9', 'N10' ...
    };
% dirname = 'season2';

if exist(dirname, 'var')
    write_scenario(...
        csv_sid, ...
        csv_xa_MW, ...
        csv_xf_MW, ...
        csv_xs_MW, ...
        header, ...
        dirname ...
        );
end

end

function idxC = cluster_analysis_tmp(x,f_cluster,K)
% Note this is a temp version of cluster_analysis, only for the 9 sites
% analysis. Should rewrite cluster_analysis.m and call it in the future.

nI = size(x, 2);

s_cluster   = nan(length(K), 1); % Silhouette metric
idxC = nan(nI, length(K));
x_t = x'; % Let's try wind speed
for i = 1: length(K)
    k = K(i);
    [idxC(:, i), ~] = feval(...
        f_cluster, x_t, k,...
        'Replicates',10,...
        'Distance', 'correlation');
end

end

function all_p_actual = picp(samples, actual, all_p_norminals)
% samples should be a nsamples x nsite x nscenarios
% actual should be nsamples x nsite
% p_norminal is the nominal coverage rate, e.g., [10, 20, ..., 90]
% should return actual coverage rate, the same length as all_p_norminals

all_p_actual = nan(size(all_p_norminals));
for i = 1: length(all_p_norminals)
    p_norminal = all_p_norminals(i);
    p_upper = 100 - (100-p_norminal)/2;
    p_lower = (100-p_norminal)/2;
    bounds = prctile(samples, [p_lower p_upper], 3);
    hitted = (actual<=bounds(:, :, 2)) & (actual>=bounds(:, :, 1));
    all_p_actual(i) = sum(hitted(:))/numel(actual)*100;
end

end

function scores = interval_score(samples, actual, all_p_norminals)
% Follow equation (32) in Wan, Can, et al. "Probabilistic forecasting of 
% wind power generation using extreme learning machine." IEEE Transactions 
% on Power Systems 29.3 (2013): 1033-1044.
% samples should be a nsamples x nsite x nscenarios
% actual should be nsamples x nsite
% p_norminal is the nominal coverage rate, e.g., [10, 20, ..., 90]
% should return actual coverage rate, the same length as all_p_norminals
% Note: nominal covergage rate is prediction interval

[nT, nI, nS] = size(samples);
sample_reshape = reshape(samples, nT*nI, nS);
actual_reshape = reshape(actual, nT*nI, 1);

scores = nan(size(all_p_norminals));
for i = 1: length(all_p_norminals)
    p_norminal = all_p_norminals(i);
    p_upper = 100 - (100-p_norminal)/2;
    p_lower = (100-p_norminal)/2;
    bounds = prctile(sample_reshape, [p_lower p_upper], 2);
    lb = bounds(:, 1);
    ub = bounds(:, 2);
    i_lower  = find(actual_reshape < lb);
    i_higher = find(actual_reshape > ub);
    i_hitted = find((actual_reshape > lb) & (actual_reshape < ub));
    
    alpha = 100 - p_norminal;
    s_lower  = -2*alpha.*(ub(i_lower) - lb(i_lower)) - 4.*(lb(i_lower) - actual_reshape(i_lower));
    s_higher = -2*alpha.*(ub(i_higher) - lb(i_higher)) - 4.*(actual_reshape(i_higher) - ub(i_higher));
    s_hitted = -2*alpha.*(ub(i_hitted) - lb(i_hitted));
    s = mean([s_lower; s_higher; s_hitted]);
    scores(i) = s;
end

end

function sizes = interval_size(samples, actual, all_p_norminals)
% Follow equation (32) in Wan, Can, et al. "Probabilistic forecasting of 
% wind power generation using extreme learning machine." IEEE Transactions 
% on Power Systems 29.3 (2013): 1033-1044.
% samples should be a nsamples x nsite x nscenarios
% actual should be nsamples x nsite
% p_norminal is the nominal coverage rate, e.g., [10, 20, ..., 90]
% should return actual coverage rate, the same length as all_p_norminals
% Note: nominal covergage rate is prediction interval

[nT, nI, nS] = size(samples);
sample_reshape = reshape(samples, nT*nI, nS);
actual_reshape = reshape(actual, nT*nI, 1);

sizes = nan(size(all_p_norminals));
for i = 1: length(all_p_norminals)
    p_norminal = all_p_norminals(i);
    p_upper = 100 - (100-p_norminal)/2;
    p_lower = (100-p_norminal)/2;
    bounds = prctile(sample_reshape, [p_lower p_upper], 2);
    lb = bounds(:, 1);
    ub = bounds(:, 2);
    sizes(i) = mean(ub-lb);
end

end