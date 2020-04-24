function sumscore = crps(x, f, x0, deltax)
% x is the 1-d bincenter of random variable, f is 1-d CDF of the random 
% variable, x0 is observation (scaler), deltax is the size of x (scaler)
n = numel(x);
score = nan(n, 1);
for i = 1: n
    x_left = x(i) - deltax/2;
    x_right = x(i) + deltax/2;
    if (x_left>=x0)
        score(i) = deltax*(f(i)-1)^2;
    elseif (x_right<=x0)
        score(i) = deltax*(f(i))^2;
    else
        score(i) = (x0-x_left)*(f(i))^2 + (x_right-x0)*(f(i)-1)^2;
    end
end
sumscore = sum(score);
end