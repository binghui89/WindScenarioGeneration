function [x_reduced, q] = forward_reduction(x, p, k)
% Use method from Dupa?ová, Jitka, Nicole Gröwe-Kuska, and Werner Römisch. 
% "Scenario reduction in stochastic programming." 
% Mathematical programming 95.3 (2003): 493-511.
% k is the number of scenarios that will be returned, i.e., we will reduce
% to this number. p is the probability of the original scenario. x includes
% the original scenario data.
[~, nS] = size(x);
c = nan(nS);
for i = 1: nS
    tmp = repmat(x(:, i), 1, nS);
    c(i, :) = vecnorm(x - tmp, 2, 1); % We use l-2 norm as the c function
end

%% Start reduction
b_selected = false(nS, 1); % Boolean variables
iS = 1: nS; % Array of indices of all scenarios
for ik = 1: k
    sigma_pc = inf.*ones(nS, 1); % Eqn. (17) from Dupa?ová et al.
    i_not_selected = iS(~b_selected); % Array of indices of not selected scenarios
    for i = i_not_selected
        b_u = b_selected;
        b_u(i) = true; % Temporarily add scenario i to set U
        sigma_pc(i) = sum(min(c(b_u, ~b_u))*p(~b_u));
    end
    [~, imin] = min(sigma_pc);
    b_selected(imin) = true; % Select scenario imin
end

%% Calculate probability
p_tmp = p;
i_not_selected = iS(~b_selected); % Indices of abandoned ones
for i = i_not_selected
    c_tmp = c(:, i); % The distance of this abondoned one to all
    c_tmp(~b_selected) = inf; % The distance of the abondoned one to selected ones
    [~, imin] = min(c_tmp); % Chose the selected one that is closest to the abondoned ones
    p_tmp(imin) = p_tmp(imin) + p(i); % Give the abondoned one's probability to the chosen selected one
end
q = p_tmp(b_selected);
x_reduced = x(:, b_selected);
end