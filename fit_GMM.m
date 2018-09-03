function [GMModel, chi_GMM] = fit_GMM(x, nbins)
%% Use Gaussian Mixture Model, returned as a gmdistribution model.
ncomp = 8; % Maximum number of components
[counts, binedges] = histcounts(x, nbins);
bincenter = (binedges(1: end-1) + binedges(2: end))/2;
binwidth = binedges(2)-binedges(1);
area = length(x)*binwidth;
f_GM_ref = transpose(counts/area);
x_hist = bincenter';
chi_GMM = nan(1, length(ncomp));
for i = 1: ncomp
    GMModel = fitgmdist(x, i);
    f_GM = pdf(GMModel, x_hist);
    chi_GMM(i) = sum((f_GM-f_GM_ref).^2./f_GM);
end

[~, iopt] = min(chi_GMM);
GMModel = fitgmdist(x, iopt);
end
