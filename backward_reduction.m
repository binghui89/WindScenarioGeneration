function [x_reduced, q, b_selected] = backward_reduction(x, p, k)
% Use method from Dupa?ová, Jitka, Nicole Gröwe-Kuska, and Werner Römisch. 
% "Scenario reduction in stochastic programming." 
% Mathematical programming 95.3 (2003): 493-511.
% k is the number of scenarios that will be returned, i.e., we will reduce
% to this number. p is the probability of the original scenario. x includes
% the original scenario data.
tic;
[~, nS] = size(x);
c = nan(nS);
for i = 1: nS
    tmp = repmat(x(:, i), 1, nS);
    c(i, :) = vecnorm(x - tmp, 2, 1); % We use l-2 norm as the c function
end

%% Start reduction
b_not_selected = false(nS, 1); % Boolean variables
iS = 1: nS; % Array of indices of all scenarios
for ik = 1: nS - k
    pc = inf.*ones(nS, 1); % Eqn. (16) from Dupa?ová et al.
    i_selected = iS(~b_not_selected);
    for i = i_selected
        b_l = b_not_selected;
        b_l(i) = true; % Temporarily add scenario i to set L
        pc(i) = p(i)*min(c(~b_l, i));
    end
    [~, imin] = min(pc);
    b_not_selected(imin) = true;
end

%% Calculate probability
p_tmp = p;
i_not_selected = iS(b_not_selected); % Indices of abandoned ones
for i = i_not_selected
    c_tmp = c(:, i); % The distance of this abondoned one to all
    c_tmp(b_not_selected) = inf; % The distance of the abondoned one to selected ones
    [~, imin] = min(c_tmp); % Chose the selected one that is closest to the abondoned ones
    p_tmp(imin) = p_tmp(imin) + p(i); % Give the abondoned one's probability to the chosen selected one
end
q = p_tmp(~b_not_selected);
x_reduced = x(:, ~b_not_selected);
b_selected = ~b_not_selected;
toc;
end