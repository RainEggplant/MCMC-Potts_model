function zeta = calc_zeta(u_mean, BETA)
zeta = zeros(1, 5);

for k = 2:5
    u = (u_mean(k) + u_mean(k - 1)) / 2;
    zeta(k) = zeta(k-1) - u * (BETA(k) - BETA(k - 1));
end
end
