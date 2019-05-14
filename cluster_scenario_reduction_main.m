function cluster_scenario_reduction_main(write_flag)
% This script is for the 2-cluster case of the complete 45-site run.
if nargin == 0
    write_flag = 0;
end
cluster_scenario_reduction4(write_flag);

end

function cluster_scenario_reduction1()
% Reduce scenarios of multiple clusters directly from their per unit 
% combinations

matfile = {'S3000_da_37sites_logitn.mat', 'S3000_da_8sites_logitn.mat'};
K = 10; % Number of scenarios reduced to
nC = length(matfile);
data = cell(nC, 1);
nI = nan(nC, 1);
for i = 1: nC
    data{i} = load(matfile{i});
    nI(i) = data{i}.nI;
end

nS1 = data{1}.nS1; nS2 = data{1}.nS2;
nT1 = data{1}.nT1; nT2 = data{1}.nT2;

x_c = cell(nC, 1);
for i = 1: nC
    x_c{i} = data{i}.x_new(:, :, nS1+1: nS1+nS2);
end
p = 1/nS2*ones(nC, nS2);

[x_reduced, q] = cluster_scenario_reduction(x_c, p, K);

% Plot
for c = 1: nC
    x_new = data{c}.x_new;
    capacity = data{c}.capacity;
    x_new_MW = nan(size(x_new));
    x_reduced_MW = nan(size(x_reduced{c}));
    xa_MW = nan(size(data{c}.xa));
    xf_MW = nan(size(data{c}.xf));
    
    for i = 1: nI(c)
        x_reduced_MW(:, i, :) = capacity(i).*x_reduced{c}(:, i, :);
        x_new_MW(:, i, :)     = capacity(i).*x_new(:, i, :);
        xa_MW(:, i) = data{c}.xa(:, i).*capacity(i);
        xf_MW(:, i) = data{c}.xf(:, i).*capacity(i);
    end
    
    fig_curve(sum(xa_MW(nT1+1: nT1+nT2, :), 2), sum(xf_MW(nT1+1: nT1+nT2, :), 2), squeeze(sum(x_reduced_MW, 2)))

end

end

function cluster_scenario_reduction2()
% Reduced each cluster to 20 scenarios, then reduce their combinations
matfile = {'S3000_da_37sites_logitn.mat', 'S3000_da_8sites_logitn.mat'};
K = 10; % Number of scenarios reduced to
nC = length(matfile);
data = cell(nC, 1);
nI = nan(nC, 1);
for i = 1: nC
    data{i} = load(matfile{i});
    nI(i) = data{i}.nI;
end

nS1 = data{1}.nS1; nS2 = data{1}.nS2;
nT1 = data{1}.nT1; nT2 = data{1}.nT2;

x_c = cell(nC, 1);
x_reduced_0 = cell(nC, 1);
q_0 = cell(nC, 1);
for i = 1: nC
    p = 1/nS2.*ones(nS2, 1);
    x_c{i} = reshape(data{i}.x_new(:, :, nS1+1: nS1+nS2), nT2*nI(i), nS2);
    [tmp, q_0{i}, ~] = reduction_forward(x_c{i}, p, K);
    x_reduced_0{i} = reshape(tmp, nT2, nI(i), K);
end

[x_reduced, q] = cluster_scenario_reduction(x_reduced_0, [q_0{1}'; q_0{2}'], K);

for c = 1: nC
    x_new = data{c}.x_new;
    capacity = data{c}.capacity;
    x_new_MW = nan(size(x_new));
    x_reduced_MW = nan(size(x_reduced{c}));
    xa_MW = nan(size(data{c}.xa));
    xf_MW = nan(size(data{c}.xf));
    
    for i = 1: nI(c)
        x_reduced_MW(:, i, :) = capacity(i).*x_reduced{c}(:, i, :);
        x_new_MW(:, i, :)     = capacity(i).*x_new(:, i, :);
        xa_MW(:, i) = data{c}.xa(:, i).*capacity(i);
        xf_MW(:, i) = data{c}.xf(:, i).*capacity(i);
    end
    
    fig_curve(sum(xa_MW(nT1+1: nT1+nT2, :), 2), sum(xf_MW(nT1+1: nT1+nT2, :), 2), squeeze(sum(x_reduced_MW, 2)))

end

end

function cluster_scenario_reduction3()
% Use actual wind power for scenario reduction
matfile = {'S3000_da_37sites_logitn.mat', 'S3000_da_8sites_logitn.mat'};
K = 10; % Number of scenarios reduced to
nC = length(matfile);
data = cell(nC, 1);
nI = nan(nC, 1);
for i = 1: nC
    data{i} = load(matfile{i});
    nI(i) = data{i}.nI;
end

nS1 = data{1}.nS1; nS2 = data{1}.nS2;
nT1 = data{1}.nT1; nT2 = data{1}.nT2;

% x_c = cell(nC, 1);
% for i = 1: nC
%     x_c{i} = reshape(data{i}.x_new(:, :, nS1+1: nS1+nS2), nT2*nI(i), nS2);
% end

x_c_MW = cell(nC, 1);
for c = 1: nC
    x_new_MW = nan(nT2, nI(c), nS2);
    for i = 1: nI(c)
        x_new_MW(:, i, :) = data{c}.x_new(:, i, nS1+1: nS1+nS2).*data{c}.capacity(i);
    end
    x_c_MW{c} = x_new_MW(:, :, :);
end
p = 1/nS2*ones(nC, nS2);

[x_reduced, q] = cluster_scenario_reduction(x_c_MW, p, K);

for c = 1: nC
    x_new = data{c}.x_new;
    capacity = data{c}.capacity;
    x_new_MW = nan(size(x_new));
    x_reduced_MW = x_reduced{c};
    xa_MW = nan(size(data{c}.xa));
    xf_MW = nan(size(data{c}.xf));
    for i = 1: nI(c)
        x_new_MW(:, i, :)     = capacity(i).*x_new(:, i, :);
        xa_MW(:, i) = data{c}.xa(:, i).*capacity(i);
        xf_MW(:, i) = data{c}.xf(:, i).*capacity(i);
    end
    fig_curve(sum(xa_MW(nT1+1: nT1+nT2, :), 2), sum(xf_MW(nT1+1: nT1+nT2, :), 2), squeeze(sum(x_reduced_MW, 2)))
end
end

function cluster_scenario_reduction4(write_flag)
% Use aggregated actual power for scenario reduction
matfile = {'S3000_da_37sites_logitn.mat', 'S3000_da_8sites_logitn.mat'};
K = 10; % Number of scenarios reduced to

if nargin == 0
    write_flag = 0;
end

nC = length(matfile);
data = cell(nC, 1);
nI = nan(nC, 1);
for i = 1: nC
    data{i} = load(matfile{i});
    nI(i) = data{i}.nI;
end

nS1 = data{1}.nS1; nS2 = data{1}.nS2;
nT1 = data{1}.nT1; nT2 = data{1}.nT2;

x_c = cell(nC, 1);
for i = 1: nC
    x_c{i} = reshape(data{i}.x_new(:, :, nS1+1: nS1+nS2), nT2*nI(i), nS2);
end

x_c_MW = cell(nC, 1);
for c = 1: nC
    tmp = nan(nT2, nI(c), nS2);
    for i = 1: nI(c)
        tmp(:, i, :) = data{c}.x_new(:, i, nS1+1: nS1+nS2).*data{c}.capacity(i);
    end
    x_c_MW{c} = tmp;
end

x_combine = combvec(...
    squeeze(sum(x_c_MW{1}, 2)), ...
    squeeze(sum(x_c_MW{2}, 2)) ...
    ); % Only supports 2 clusters.
p = 1/size(x_combine, 2)*ones(size(x_combine, 2), 1);
sindex = combvec(1: 100, 1: 100);

%% Heuristic method
[x_reduced, q, b_selected] = reduction_forward(x_combine, p, K);
for i = 1: length(q)
    fprintf('Scenario %g: %f\n', i, q(i))
end

%% Return scenario specific data
x_reduced_MW = cell(nC, 1);
for c = 1: nC
    tmp = data{c}.x_new(:, :, nS1+1: nS1+nS2);
    for i = 1: data{c}.nI
        tmp(:, i, :) = data{c}.x_new(:, i, nS1+1: nS1+nS2).*data{c}.capacity(i);
    end
    x_reduced_MW{c} = tmp(:, :, sindex(c, b_selected));
end

%% Plot figures
for c = 1:nC
    figure();
    hold on;
%     nI = data{c}.nI;
    x_new = data{c}.x_new;
    capacity = data{c}.capacity;
    x_new_MW = nan(size(x_new));
%     x_reduced_MW = nan(size(x_reduced{c}));
    xa_MW = nan(size(data{c}.xa));
    xf_MW = nan(size(data{c}.xf));
    
    % To map the wind site to wind farms in the TX2kB system, we use
    % capacity multiplier here
    TXdata;
    sid_multiplier = nan(nI(c), 1);
    for i = 1: nI(c)
        iselected = (SITE_ID==data{c}.sid(i));
        sid_multiplier(i) = sum(PMAX(iselected)./SITE_CAP(iselected));
    end

    for i = 1: nI(c)
        x_reduced_MW{c}(:, i, :) = x_reduced_MW{c}(:, i, :).*sid_multiplier(i);
        x_new_MW(:, i, :) = x_new(:, i, :).*capacity(i).*sid_multiplier(i);
        xa_MW(:, i) = data{c}.xa(:, i).*capacity(i).*sid_multiplier(i);
        xf_MW(:, i) = data{c}.xf(:, i).*capacity(i).*sid_multiplier(i);
    end
    x_new_MWsum = squeeze(sum(x_new_MW(:, :, nS1+1: nS1+nS2), 2));
    x_reduced_MWsum = squeeze(sum(x_reduced_MW{c}, 2));
    h1 = plot(1: nT2, sum(xa_MW(nT1+1: nT1+nT2, :), 2), 'LineWidth', 2, 'Color', 'r');
    h2 = plot(1: nT2, sum(xf_MW(nT1+1: nT1+nT2, :), 2), 'LineWidth', 2, 'Color', 'm');
%     h = plot(1: size(x_new_MWsum, 1), x_new_MWsum, 'Color', [102, 170, 215]./255);
%     h3 = h(end);
    plot(1: size(x_reduced_MWsum, 1), x_reduced_MWsum, 'Color', 'k');

    up_prcnt = [0.95, 0.75, 0.55];
    dn_prcnt = [0.05, 0.25, 0.45];
    alphas   = [0.1, 0.2, 0.3];
    h_band = nan(size(up_prcnt));
    for j = 1: length(up_prcnt)
        x_new_dn = nan(nT2, nI(c));
        x_new_up = nan(nT2, nI(c));
        for t = 1: 24
            for i = 1: nI(c)
                [f, x] = ecdf(squeeze(x_new_MW(t, i, nS1+1: nS1+nS2)));
                [~, i_dn] = min(abs(f-dn_prcnt(j)));
                x_new_dn(t, i) = x(i_dn);
                [~, i_up] = min(abs(f-up_prcnt(j)));
                x_new_up(t, i) = x(i_up);
            end
        end
        x_new_dn_sum = sum(x_new_dn, 2);
        x_new_up_sum = sum(x_new_up, 2);
        x_fill = [1:24, fliplr(1:24)];
        y_fill = [x_new_dn_sum', fliplr(x_new_up_sum')];
        fh = fill(x_fill, y_fill, 'k', 'FaceAlpha', alphas(j));
        fh.EdgeColor = 'none';
        h_band(j) = fh;
    end
    box on;
    hold off;
    xlim([0, 25]);
    xlabel('Time (h)');
end

% Write results to csv files
csv_sid = [data{1}.sid; data{2}.sid];
csv_xa_MW = nan(nT2, numel(csv_sid));
csv_xf_MW = nan(nT2, numel(csv_sid));
nI_pass = 0;
for c = 1: nC
    for i = 1: data{c}.nI
        csv_xa_MW(:, nI_pass+i) = data{c}.xa(nT1+1: nT1+nT2, i).*data{c}.capacity(i);
        csv_xf_MW(:, nI_pass+i) = data{c}.xf(nT1+1: nT1+nT2, i).*data{c}.capacity(i);
    end
    nI_pass = nI_pass + data{c}.nI;
end
header = {...
    'xa', 'xf', ...
    'N1', 'N2', 'N3', 'N4', 'N5', ...
    'N6', 'N7', 'N8', 'N9', 'N10' ...
    };

if write_flag == 1
    write_scenario(...
        csv_sid, ...
        csv_xa_MW, ...
        csv_xf_MW, ...
        cat(2, x_reduced_MW{1}, x_reduced_MW{2}), ...
        header ...
        );
end
end


function [x_reduced, q] = cluster_scenario_reduction(x_new, p, K)
%% Reduce scenarios of multiple clusters, now only supports 2 clusters.
% x_new: cell array, dim = (nT, nI, nS), scenario data
% p: array, dim = (nC, nS), probability
% K: Number of scenarios that will be reduced to

nC = numel(x_new);
nI = nan(nC, 1);
nT = size(x_new{1}, 1);
nS = size(x_new{1}, 3);
for c = 1: nC
    nI(c) = size(x_new{c}, 2);
end

x_flat_c = cell(nC, 1);
for i = 1: nC
    tmp = x_new{i};
    x_flat_c{i} = reshape(tmp, nT*nI(i), nS);
end

x_combine = combvec(x_flat_c{1}, x_flat_c{2}); % Only supports 2 clusters.
p_combine = prod(combvec(p(1, :), p(2,: )), 1);

%% Heuristic method
[x_reduced_flat, q, b_selected] = reduction_forward(x_combine, p_combine(:), K);

%% Return scenario specific data
x_reduced_together = reshape(x_reduced_flat, nT, sum(nI), K );

x_reduced = cell(nC, 1);
for i = 1: nC
    if i == 1
        x_reduced{i} = x_reduced_together(:, 1: nI(i), :);
    else
        x_reduced{i} = x_reduced_together(:, nI(i-1)+1: sum(nI(1:i)), :);
    end
end

end
