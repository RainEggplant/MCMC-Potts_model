clear;
%% 常数、参数定义
K = 20; % 20*20 网格
N_Q = 10; % 10 种取值
BETA_ARRAY = [1.4; 1.4065; 1.413; 1.4195; 1.426]; % 1/T 的值
BURN_IN_FACTOR = 0.3; % burn in 样本数占比

%%
N_MH = 200000;
u_mh = zeros(N_MH, 5);
x_mh = zeros(K, K, 5);
n_skip_mh = floor(N_MH * BURN_IN_FACTOR);

% 并行计算，若要串行请将 parfor 改为 for
parfor k = 1:5
    [u_mh(:, k), x_mh(:, :, k)] = ...
        metropolis(K, N_Q, BETA_ARRAY(k), N_MH);
end

%%
zeta_mh = calc_zeta(mean(u_mh(n_skip_mh:end, :)), BETA_ARRAY);
disp(['MH: Zeta = ', num2str(zeta_mh)]);
plot_histogram(u_mh(n_skip_mh:end,:), K)
plot_sample(x_mh, N_Q);

%%
N_SW = 200000;
u_sw = zeros(N_SW, 5);
x_sw = zeros(K, K, 5);
n_skip_sw = floor(N_MH * BURN_IN_FACTOR);
% parfor k = 1:5
%     [u_sw(:, k), x_sw(:, :, k)] = ...
%         swendsen_wang(K, N_Q, BETA_ARRAY(k), N_SW);
% end

%%
zeta_sw = calc_zeta(mean(u_sw(n_skip_sw:end, :)), BETA_ARRAY);
disp(['SW: Zeta = ', num2str(zeta_sw)]);
plot_histogram(u_sw(n_skip_sw:end,:), K)
plot_sample(x_sw, N_Q);

%%
N_SAMC = 2200000;
N_BURN_IN = 200000;
SUBSAMPLE_STEP = 10;

[zeta_rec, u_index, u_rec] = ...
    samc(K, N_Q, BETA_ARRAY, N_SAMC, N_BURN_IN, SUBSAMPLE_STEP);
disp(['SAMS: Zeta = ', num2str(zeta_rec(:, end)')]);

%%
figure;
hold on;
title('u(x)/K 直方图');
subplot(1, 3, 1);
hold on;
for l = 1:3
    ksdensity(u_rec(l, ceil(N_BURN_IN / SUBSAMPLE_STEP / 5):u_index(l)) / K ^ 2);
end
legend('T_1', 'T_2', 'T_3');
xlabel('u(x)/K');
ylabel('Density');
xlim([-2, -0.5]);
ylim([0, 6]);

subplot(1, 3, 2);
ksdensity(u_rec(4, ceil(N_BURN_IN / SUBSAMPLE_STEP / 5):u_index(4)) / K ^ 2);
legend('T_4');
xlabel('u(x)/K');
ylabel('Density');
xlim([-2, -0.5]);
ylim([0, 6]);

subplot(1, 3, 3);
ksdensity(u_rec(5, ceil(N_BURN_IN / SUBSAMPLE_STEP / 5):u_index(5)) / K ^ 2);
legend('T_4');
xlabel('u(x)/K');
ylabel('Density');
xlim([-2, -0.5]);
ylim([0, 6]);
