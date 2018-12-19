function x_new = gibbs2(data, nT1, nT2, nS1, nS2, selected_id)
% Note selected_id should be a column vector of logical variables having the
% same length as data.sid

% Data preparation
fprintf('Gibbs sampling starts at: %s\n', datestr(now));
% data = read_all();
nbins = 50;
% nT1 = 46968 - 24; % Total number of hours for training
% nT2 = 24; % Total number of hours for forcasting
nT = nT1 + nT2; % Total number of samples, total hours from 1/1/07 to 5/10/12
% nS = 3000; nS1 = 2900; nS2 = nS - nS1; % nS1 is # of burn-in scenarios
nS = nS1 + nS2;

xa = data.xa(:, selected_id);
xf = data.xf(:, selected_id);
sid = data.sid(selected_id);
capacity = data.cap(selected_id);
nI = length(sid);

% Logit transformation
ya = log( xa(1: nT, :)./(1 - xa(1: nT, :)) );
yf = log( xf(1: nT, :)./(1 - xf(1: nT, :)) );

uam = nan(nT, nI); % Marginal CDF of actual wind power
ufm = nan(nT, nI); % Marginal CDF of forecasted wind power
for i = 1: nI
    gm_a = fitdist(ya(~isinf(ya(1: nT1, i)), i),'Normal'); % Only use time 1 to nT1 for training
    gm_f = fitdist(yf(~isinf(yf(1: nT1, i)), i),'Normal');
    struct_a(i).Fm = gm_a; % Fm is fitted marginal distribution function
    struct_f(i).Fm = gm_f;
    uam(:, i) = cdf(gm_a, ya(:, i));
    ufm(:, i) = cdf(gm_f, yf(:, i));
end
uam(uam==1) = 1-eps;
uam(uam==0) = eps;
ufm(ufm==1) = 1-eps;
ufm(ufm==0) = eps;

% Find parameter r for the covariance matrix
ra = nan(nI, 1);
figure();
hold on;
for i = 1: nI
    ra(i) = find_r( uam(:, i), 1 );
end
hold off;

% Find copula
fprintf('Copula models started at: %s\n', datestr(now));
I = 1: nI;
parfor i = 1: nI
    % Note in the copula package, the conditioned variable is always the
    % last one
    % C(ufm1, ufm2, ..., ufmI, uam1, uam2, ... uam(i-1), uam(i+1), ...uamI, uami)
    u_copula_input = [ufm(1:nT1, :) uam(1:nT1, I~=i) uam(1:nT1, i)];
    copula_model = func_choose_best_copula_model(u_copula_input);
    struct_a(i).copula = copula_model;
end
fprintf('Copula models done at: %s\n', datestr(now));

uam_new = scenario_gen(nI, nS, nT1, nT2, ufm, struct_a, ra);
fprintf('Gibbs sampling stops at: %s\n', datestr(now));

% Make a graph
y_new = nan(size(uam_new));
x_new = nan(size(uam_new)); % The actual p.u. output of all scenarios
for i = 1: nI
    figure();
    hold on;
    for s = 1: nS % The first nS1 scenarios are for burn-in process
        y_tmp = norminv(uam_new(:, i, s), struct_a(i).Fm.mean, struct_a(i).Fm.sigma);
%         y_tmp = interp1(u_ascend, x_ascend, rescale(uam_new(:, i, s), u_ascend(1), u_ascend(end)), 'linear');
        x_tmp = 1./(1+exp(-y_tmp));
        y_new(:, i, s) = y_tmp;
        x_new(:, i, s) = x_tmp;
        if s > nS1 % The first nS1 scenarios are for burn-in process
            plot(1: length(x_tmp), x_tmp, 'Color', [102, 170, 215]./255);
        end
    end
    h1 = plot(1: nT2, xa(nT1+1: nT, i), 'LineWidth', 2, 'Color', 'r');
    h2 = plot(1: nT2, xf(nT1+1: nT, i), 'LineWidth', 2, 'Color', 'm');
    hold off;
    legend([h1, h2], 'Actual', 'DA forecast');
    xlabel('Time');
    ylabel('Power (p.u.)');
    title(sid(i));
end

% Total production
x_new_MW = nan(size(x_new));
for i = 1: nI
    x_new_MW(:, i, :) = capacity(i).*x_new(:, i, :);
end
figure();
hold on;
plot(1: nT2, sum(data.xa_MW(nT1+1: nT, selected_id), 2), 'LineWidth', 2, 'Color', 'r');
plot(1: nT2, sum(data.xf_MW(nT1+1: nT, selected_id), 2), 'LineWidth', 2, 'Color', 'm');
for s = nS1+1: nS
    plot(1: size(x_new_MW, 1), squeeze(sum(x_new_MW(:, :, nS1+1: nS), 2)), 'Color', [102, 170, 215]./255);
end

end

function uam_new = scenario_gen(nI, nS, nT1, nT2, ufm, struct_a, ra)
% Generate scenarios
I = 1: nI;
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

end