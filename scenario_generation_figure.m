function scenario_generation_figure(season)
% The figures in the journal paper
figure_per_season(season);
% sem_multiple_figures()
end

function figure_per_season(season)
switch season
    case 1
        s_a = load('full_case.S9000_da_45sites_logitn.s1.mat'); % K = 2 to 20
        s_b = load('season1_20191031.mat'); % K = 1
    case 2
        s_a = load('full_case.S9000_da_45sites_logitn.s2.mat');
        s_b = load('season2_20191031.mat');
    case 3
        s_a = load('full_case.S9000_da_45sites_logitn.s3.mat');
        s_b = load('season3_20191031.mat');
    case 4
        s_a = load('full_case.S9000_da_45sites_logitn.s4.mat');
        s_b = load('season4_20191031.mat');
end
Karray = [s_b.Karray(1) s_a.Karray];
idxC = cat(2, s_b.idxC(:, 1, :), s_a.idxC);
xnew_all = cat(1, {s_b.xnew_all{1}}, s_a.xnew_all);
data_s = s_b.data_s;
nS1 = s_b.nS1;
nS2 = s_b.nS2;
nT1 = s_b.nT1;
nT2 = s_b.nT2;
time_taken_sc = cat(1, s_b.time_taken_sc(1), s_a.time_taken_sc);
time_taken_red = cat(2, [s_b.time_taken_red(1, 1); zeros(19, 1)], s_a.time_taken_red);
clear s_a s_b;

%% PICP analysis
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

delta_picp_pinc = abs(...
    array_average_p_actual - repmat(10:10:90, size(array_average_p_actual, 1), 1)...
    );

% PICP figure 1, all PICP curves
figure();
plot([0 p_norminal 100]', [zeros(length(Karray), 1) array_average_p_actual 100.*ones(length(Karray), 1)]', 'color', [0.5, 0.5, 0.5]);
% legend(cell_legends);
[~, rowmax] = max(mean(delta_picp_pinc, 2));
[~, rowmin] = min(mean(delta_picp_pinc, 2));
hold on;
h1 = plot([0 p_norminal 100], [0 array_average_p_actual(rowmin, :) 100], 'r');
h2 = plot([0 p_norminal 100], [0 array_average_p_actual(rowmax, :) 100], 'b');
plot(0:10:100, 0:10:100, '--k');
set(findall(gcf,'-property','FontSize'),'FontSize',22);
cell_legend = {strcat('K=', int2str(Karray(rowmin))), strcat('K=', int2str(Karray(rowmax)))};
legend([h1; h2], cell_legend);
legend boxoff;
xlabel('PINC (%)');
ylabel('PICP (%)');

% PICP figure 2, average delta PICP vs. nK
figure()
plot(Karray, mean(delta_picp_pinc, 2), 'LineWidth', 2, 'Marker', 's', 'Color', 'k');
xlabel('Number of clusters (K)');
ylabel('ACE (%)');
set(findall(gcf,'-property','FontSize'),'FontSize',22);

% PICP figure 3, max and min average delta PICP vs. nK
[~, rowmax] = max(mean(delta_picp_pinc, 2));
[~, rowmin] = min(mean(delta_picp_pinc, 2));
figure();
plot([0 p_norminal 100], [0 array_average_p_actual(rowmin, :) 100], 'r');
hold on;
plot([0 p_norminal 100], [0 array_average_p_actual(rowmax, :) 100], 'b');
plot(0:10:100, 0:10:100, 'k');
xlabel('Nominal coverage rate');
ylabel('Actual coverage rate');
cell_legend = {strcat('K=', int2str(Karray(rowmin))), strcat('K=', int2str(Karray(rowmax)))};
legend(cell_legend);

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
plot(p_norminal', array_average_interval_scores'.*100, 'color', [0.5, 0.5, 0.5]);
[~, rowmax] = max(mean(array_average_interval_scores, 2));
[~, rowmin] = min(mean(array_average_interval_scores, 2));
hold on;
h1 = plot(p_norminal, array_average_interval_scores(rowmin, :).*100, 'b');
h2 = plot(p_norminal, array_average_interval_scores(rowmax, :).*100, 'r');
xlabel('PINC (%)');
ylabel('Interval score (%)');
% legend(cell_legends);
cell_legend = {strcat('K=', int2str(Karray(rowmin))), strcat('K=', int2str(Karray(rowmax)))};
legend([h1; h2], cell_legend);
legend boxoff;
set(findall(gcf,'-property','FontSize'),'FontSize',22);

% Sharpness figure 2, average sharpness vs. nK
figure()
plot(Karray(:), mean(array_average_interval_scores, 2).*100, 'LineWidth', 2, 'Marker', 's', 'Color', 'k');
xlabel('Number of clusters (K)');
ylabel('AIS (%)');
set(findall(gcf,'-property','FontSize'),'FontSize',22);

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
plot([0 p_norminal 100]', [zeros(length(Karray), 1) array_average_interval_sizes.*100 100.*ones(length(Karray), 1)]');
xlabel('Nominal coverage rate');
ylabel('Interval size');
legend(cell_legends);

% Interval size figure 2, average sharpness vs. nK
figure()
delta_interval_size = abs(array_average_interval_sizes.*100 - repmat(10:10:90, size(array_average_interval_sizes, 1), 1));
% plot(Karray(:), mean(array_average_interval_sizes, 2));
plot(Karray, mean(delta_interval_size, 2));
xlabel('Number of clusters (K)');
ylabel('Average delta interval size');

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

%% Comprehensive, including reliability and interval size
figure();
% average_coverage_and_size = (delta_interval_size + delta_picp_pinc)./2;
plot(Karray, 0.5.*mean(delta_picp_pinc, 2) - 0.5.*mean(array_average_interval_scores.*100, 2), 'LineWidth', 2, 'Marker', 's', 'Color', 'k');
xlabel('Number of clusters (K)');
ylabel('SEM (%)');
set(findall(gcf,'-property','FontSize'),'FontSize',22);


%% Time analysis
figure();
h = bar(Karray, [time_taken_sc sum(time_taken_red, 1)';]./60, 'stacked');
set(h(1), 'FaceColor',[.5, .5, .5])
set(h(2), 'FaceColor', 'w')
legend('Scenario generation' ,'Scenario reduction')
legend boxoff;
xlabel('Number of clusters (K)')
ylabel('Time (minutes)')
set(findall(gcf,'-property','FontSize'),'FontSize',22);

end

function sem_multiple_figures()
figure();
hold on;
for season = 1: 4

    switch season
        case 1
            load('full_case.S9000_da_45sites_logitn.s1.mat');
        case 2
            load('full_case.S9000_da_45sites_logitn.s2.mat');
        case 3
            load('full_case.S9000_da_45sites_logitn.s3.mat');
        case 4
            load('full_case.S9000_da_45sites_logitn.s4.mat');
    end

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

    delta_picp_pinc = abs(...
        array_average_p_actual - repmat(10:10:90, size(array_average_p_actual, 1), 1)...
        );
    
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
    
    h(season) = plot(Karray, 0.5.*mean(delta_picp_pinc, 2) - 0.5.*mean(array_average_interval_scores.*100, 2));
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
% samples should be a ntime x nsite x nscenarios
% actual should be ntime x nsite
% p_norminal is the nominal coverage rate in percentage, e.g., [10, 20, ..., 90]
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
    
    alpha = (100 - p_norminal)/100;
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
