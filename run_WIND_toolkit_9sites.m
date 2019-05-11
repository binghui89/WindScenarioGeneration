function run_WIND_toolkit_9sites()
addpath('./packages/copula-matlab');

data = read_all();

% This is Mucun's 9 sites
% selected_sites = [1342,2061,4816,8979,9572,10069,10526,10527,11038];

% This is randomly selected 20 sites
selected_sites = randsample(data.sid, 9);

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
Karray = 2: 1: sum(I_selected); % # of clusters
f_cluster = 'kmeans';
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

% plot(Karray, s_cluster);


%% Start generating scenarios
xnew_all = cell(length(Karray), 1);

mse_allK = cell(length(Karray), 1);
mae_allK = cell(length(Karray), 1);
for nK = Karray
    idx = idxC(:, Karray==nK, 1); % Just use cluster results from the first year
    mse = nan(nK, 1);
    mae = nan(nK, 1);
    xnewK = cell(nK, 1);
    for k = 1: nK
        i_cluster = (idx==k);
    %     xa = data_s.xa(:, i_cluster);
    %     xf = data_s.xf(:, i_cluster);
    %     sid = data_s.sid(i_cluster);
        nT1 = 46968 - 24; nT2 = 24;
        nS1 = 2900; nS2 = 100;
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
end

mse_mean = nan(length(Karray), 1);
mae_mean = nan(length(Karray), 1);
for i = 1: length(Karray)
    mse_mean(i) = mean(mse_allK{i});
    mae_mean(i) = mean(mae_allK{i});
end

%% And, don't forget when K = 1
i_cluster = true(data_s.nI, 1);
xnew = gibbs2(data_s, nT1, nT2, nS1, nS2, i_cluster);
xnew_all = [{xnew};xnew_all];
abserr = abs(repmat(data_s.xa(nT1+1:nT1+nT2, i_cluster), 1, 1, nS2) - xnew(:, :, nS1+1: nS1+nS2));
% abserr = abs(data_s.xa(nT1+1:nT1+nT2, i_cluster) - xnew);
sqrerr = abserr.^2;

mse_mean = [sqrt(mean(sqrerr(:))); mse_mean];
mae_mean = [mean(abserr(:)); mae_mean];

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
ax = subplot(1, 1, 1);
% if length(mse_mean) == length(Karray)
%     plotyy(Karray, s_cluster(:, 1), Karray, mse_mean);
% elseif length(mse_mean) == length(Karray) + 1
%     plotyy(Karray, s_cluster(:, 1), [1 Karray], mse_mean);
% end
% figure();
if length(mse_mean) == length(Karray)
    [AX,H1,H2] = plotyy(Karray, s_cluster(:, 1), Karray, mse_mean);
elseif length(mse_mean) == length(Karray) + 1
%     [AX,H1,H2] = plotyy(Karray, s_cluster(:, 1), [1 Karray], mse_mean);
    ax = subplot(1, 1, 1);
    yyaxis(ax, 'left'); 
    plot(Karray, s_cluster(:, 1), 'Marker', 's', 'LineWidth', 2, 'Color', 'k');
    ylabel('Silhouette coefficient') % left y-axis
    yyaxis(ax, 'right');
    plot([1 Karray], mse_mean, 'LineStyle', '--', 'Marker', 's', 'LineWidth', 2,'Color', 'k');
    ylabel('RMSE');
end
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