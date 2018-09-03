function cdf_R = cdf_rand_correlated(dimension, r, nsample)
% Assume expectation MU are a n-dim vector and the covariance matrix SIGMA
% is calculated using r: SIGMA(i, j) = exp(-(i-j)^2/r).
% Given expectation nsample x dimension MU and dimension x dimension 
% covariance matrix SIGMA, first generate n random variables tha follow 
% dimension-dim Gaussian distribution and then return the CDF of them.
MU = zeros(nsample, dimension);
[m, n] = meshgrid(1: dimension);
SIGMA = exp(-(m-n).^2/r);
R = mvnrnd(MU, SIGMA);
cdf_R = normcdf(R);
end