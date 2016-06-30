function pval = spiegel_test(x)
% compute pvalue under null of x normally distributed;
% x should be a vector;
xm = mean(x);
xs = std(x);
xz = (x - xm) ./ xs;
xz2 = xz.^2;
N = sum(xz2 .* log(xz2));
n = numel(x);
ts = (N - 0.73 * n) / (0.8969 * sqrt(n)); %under the null, ts ~ N(0,1)
pval = 1 - abs(erf(ts / sqrt(2)));    %2-sided test.
