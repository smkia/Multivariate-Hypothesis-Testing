%
% filt = GaussianFilter(fsize, sigma)
%
% create 2D Gaussian weight filter
% Filter is normalised so that the weights sum to one.
% size(filt) = fsize x fsize
% 'sigma' controls the Gaussian width.
% Suggest: fsize should be an odd number (3,5,7,...)
% Suggest: fsize should >= 4*sigma
%
function filt = GaussianFilter(fsize, sigma)

N = fsize; N2 = round((N+1)/2);
filt = zeros(N,N);
for j = 1:N
    for i = 1:N
        u = i-N2;
        v = j-N2;
        filt(j,i) = exp(-(u^2+v^2)/(2*sigma^2));
    end
end
filt = filt/(sum(sum(filt)));


