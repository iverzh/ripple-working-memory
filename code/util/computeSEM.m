
function sem = computeSEM(x, dim)


if ndims(x) > 1 && dim == 1
    x = x(:);
    N = sum(~isnan(x));
else
    N = size(x,1);
end

sem = std(x, 'omitnan') /  sqrt(N);


end











