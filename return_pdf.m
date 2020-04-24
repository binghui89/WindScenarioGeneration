function [pdf, bincenter, x_edges] = return_pdf(x, bin_width)
x_min = floor(min(x)./bin_width)*bin_width;
x_max = ceil(max(x)./bin_width)*bin_width;
[x_counts, x_edges] = histcounts(x, x_min:bin_width:x_max);
bincenter = (x_edges(1: end-1) + x_edges(2: end))/2;
pdf = x_counts./sum(x_counts)./bin_width;
end