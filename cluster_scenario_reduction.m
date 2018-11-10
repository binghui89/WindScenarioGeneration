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
