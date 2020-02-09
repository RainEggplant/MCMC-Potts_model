function plot_histogram(u, K)
figure;
title('u(x)/K Ö±·½Í¼');
subplot(1, 3, 1);
hold on;
ksdensity(u(:, 1) / K ^ 2);
ksdensity(u(:, 2) / K ^ 2);
ksdensity(u(:, 3) / K ^ 2);
legend('T_1', 'T_2', 'T_3');
xlabel('u(x)/K');
ylabel('Density');
xlim([-2, -0.5]);
ylim([0, 6]);

subplot(1, 3, 2);
ksdensity(u(:, 4) / K ^ 2);
legend('T_4');
xlabel('u(x)/K');
ylabel('Density');
xlim([-2, -0.5]);
ylim([0, 6]);

subplot(1, 3, 3);
ksdensity(u(:, 5) / K ^ 2);
legend('T_5');
xlabel('u(x)/K');
ylabel('Density');
xlim([-2, -0.5]);
ylim([0, 6]);
    %ksdensity(u_rec(l, ceil(N_BURN_IN / SUBSAMPLE_STEP / N_BETA):u_index(l)) / K ^ 2);
end

