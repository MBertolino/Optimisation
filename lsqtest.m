rng default % for reproducibility
d = linspace(0,3);
y = exp(-1.3*d) + 0.05*randn(size(d));

fun = @(r)exp(-d*r)-y;

x0 = 4;
x = lsqnonlin(fun,x0)

xp = linspace(0,2);

for ii = 1:100
    yp(ii) = sum(fun(xp(ii)).^2);
end
plot(xp, yp)