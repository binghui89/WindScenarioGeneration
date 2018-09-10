%% Data preparation
dirwork = 'test';
dirhome = pwd;
nbins = 50;

sid = [20481, 21055, 22235];
capacity = [16, 14, 16]; % Should write a function to read capacity based on site #
nI = length(sid);
nT1 = 46968 - 24; % Total number of hours for training
nT2 = 24; % Total number of hours for forcasting
nT = nT1 + nT2; % Total number of samples, total hours from 1/1/07 to 5/10/12

%% Read data
xa_MW = nan(nT, nI); % Actual wind power in MW
xf_MW = nan(nT, nI); % Forecasted wind power in MW
xa    = nan(nT, nI); % Actual wind power in p.u.
xf    = nan(nT, nI); % Forecasted wind power in p.u.

cd(dirwork);
for i = 1: nI
    fname = strcat(int2str(sid(i)), '.csv');
%     fname = strcat('Forecast_prob_', int2str(sid(i)), '.csv');
    M = csvread(fname, 1, 0);
    xa_MW(:, i) = M(1:nT, 6);
    xf_MW(:, i) = M(1:nT, 7);
    xa(:, i)    = M(1:nT, 6)./capacity(i);
    xf(:, i)    = M(1:nT, 7)./capacity(i);
end
cd(dirhome);

%% Fit marginal distribution and get marginal CDF

% Fmxa = cell(nI, 1); % Fitted marginal distribution of xa
% Fmxf = cell(nI, 1); % Fitted marginal distribution of xf
uam = nan(nT, nI); % Marginal CDF of actual wind power
ufm = nan(nT, nI); % Marginal CDF of forecasted wind power
for i = 1: nI
    [gmm_a, ~] = fit_GMM(xa(1:nT1, i), nbins); % Only use time 1 to nT1 for training
    [gmm_f, ~] = fit_GMM(xf(1:nT1, i), nbins);
    struct_a(i).Fm = gmm_a; % Fm is fitted marginal distribution function
    struct_f(i).Fm = gmm_f;
    uam(:, i) = cdf(gmm_a, xa(:, i));
    ufm(:, i) = cdf(gmm_f, xf(:, i));
end
% [counts, binedges] = histcounts(xf, nbins);
% bincenter = (binedges(1: end-1) + binedges(2: end))/2;
% binwidth = binedges(2)-binedges(1);
% x_hist = bincenter';
% figure();
% histogram(w_f_pu, 50);
% hold on;
% plot(x_hist, pdf(gmm_model, x_hist).*sum(counts).*binwidth);
% cdf_f_pu = cdf(gmm_model, w_f_pu);

%% Find parameter r for the covariance matrix
ra = nan(nI, 1);
figure();
hold on;
for i = 1: nI
    ra(i) = find_r( uam(:, i), 1 );
end
hold off;

%% Find copula
I = 1: nI;
for i = 1: nI
    % Note in the copula package, the conditioned variable is always the
    % last one
    % C(ufm1, ufm2, ..., ufmI, uam1, uam2, ... uam(i-1), uam(i+1), ...uamI, uami)
    u_copula_input = [ufm(1:nT1, :) uam(1:nT1, I~=i) uam(1:nT1, i)];
    copula_model = func_choose_best_copula_model(u_copula_input);
    struct_a(i).copula = copula_model;
end

%% Generate scenarios
nS = 1000; nS1 = 900; nS2 = nS - nS1; % nS1 is # of burn-in scenarios
tmp = nan(nS, nT2, nI);
for i = 1: nI
    tmp(:, :, i) = cdf_rand_correlated(nT2, ra(i), nS);
end
Pa_conditionedelse = permute(tmp, [2, 3, 1]); % nT2 * nI * nS
uam_new = nan(size(Pa_conditionedelse));

u_sorted = (0:0.001:1)';
nsorted = size(u_sorted, 1);

parfor t = 1: nT2
    tic;
    conditioning_variable_ufm = ufm(nT1+t, :);
    conditioning_variable_ufa = ufm(nT1+t, :); % Use forecasted values as the initial conditioning variables of ufa
    for s = 1: nS
        for i = 1:nI
            % Gibbs sampling
            conditioning_variable = [
                conditioning_variable_ufm conditioning_variable_ufa(I~=i)
                ];
            copula_input = [
                repmat(conditioning_variable, nsorted, 1), u_sorted
                ];
            copula_model = struct_a(i).copula;
            Pa_sorted = copula.cnd(copula_model, copula_input, 2*nI);
            Pa_sorted(u_sorted==0) = 0;
            Pa_sorted(u_sorted==1) = 1;
            Pa_sorted(Pa_sorted<=0) = 0;
            Pa_sorted(Pa_sorted>=1) = 1; % Not sure if this is the correct way, but this is the only way to do it
            i_valid = (Pa_sorted>=0)&(Pa_sorted<=1)&(~isnan(Pa_sorted));
            if sum(isnan(Pa_sorted))
                disp('Pa_sorted has nan');
            end
            tmp = P2x(u_sorted(i_valid), Pa_sorted(i_valid), Pa_conditionedelse(t, i, s));
            uam_new(t, i, s) = tmp;
            conditioning_variable_ufa(i) = tmp;
            if isnan(tmp)
                disp(tmp);
            end
        end
    end
    toc;
end

%% Make a graph
x_ascend = (0:0.001:1)';
for i = 1: nI
    figure();
    hold on;
    for s = nS1+1: nS % The first nS1 scenarios are for burn-in process
        u_tmp = uam_new(:, i, s);
        u_ascend = cdf(struct_a(i).Fm, x_ascend);
        x = interp1(u_ascend, x_ascend, u_tmp, 'linear');
        plot(1: length(x), x, 'Color', [102, 170, 215]./255);
    end
    h1 = plot(1: nT2, xa(nT1+1: nT, i), 'LineWidth', 2, 'Color', 'r');
    h2 = plot(1: nT2, xf(nT1+1: nT, i), 'LineWidth', 2, 'Color', 'm');
    hold off;
    legend([h1, h2], 'Actual', 'DA forecast');
    xlabel('Time');
    ylabel('Power (p.u.)');
end
