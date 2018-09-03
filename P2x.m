function x0 = P2x(x, y, y0)
% Given a vector x and f(x), return the inverse of f^-1(y), x and fx should
% be the same dimension, y is a scaler
% if abs(y0-0)<=1E-6
%     x0 = 0;
% elseif abs(y0-1)<=1E-6
%     x0 = 1;
% else
    d = y - y0;
    i_plus  = (d>=0);
    i_minus = (d<=0);
    y_plus  = y(i_plus);
    y_minus = y(i_minus);
    x_plus  = x(i_plus);
    x_minus = x(i_minus);
    y2 = y_plus(1);
    y1 = y_minus(end);
    x2 = x_plus(1);
    x1 = x_minus(end);
    deltay = y2 - y1;
    if abs(deltay) <= 1E-6
        x0 = (x1+x2)/2;
    else
        x0 = x1 + (y0-y1)/deltay*(x2-x1);
    end
% end
end