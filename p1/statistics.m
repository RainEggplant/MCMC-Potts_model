% 统计不同举荐分布的差别
%% 常数定义
clear;
N = [100, 200, 200000]; % 迭代次数
BURN_IN_FACTOR = 0.3;
mu = [5; 10];
sigma = [1, 1; 1, 4];
n_per_N = 100;
l_N = length(N);

%% 并行计算
rho = zeros(l_N, n_per_N);
rho1 = zeros(l_N, n_per_N);
rho2 = zeros(l_N, n_per_N);
accept_rate1 = zeros(l_N, n_per_N);
accept_rate2 = zeros(l_N, n_per_N);

for k = 1:l_N
    N_k = N(k);
    parfor l = 1:n_per_N
        [rho(k, l), rho1(k, l), accept_rate1(k, l), ...
            rho2(k, l), accept_rate2(k, l)] = ...
            metropolis_hastings(N_k, BURN_IN_FACTOR, mu, sigma, 0);
    end
end

%% 输出
disp(['迭代次数为 ', num2str(N), ' 次时：' ]);
disp(['均匀分布举荐-MH 的平均偏差为：', ...
    num2str(mean(abs(rho1-rho), 2)')]);
disp(['均匀分布举荐-MH 的平均接受率为：', num2str(mean(accept_rate1, 2)')]);
disp(['二维高斯分布举荐-MH 的平均偏差为：', ...
    num2str(mean(abs(rho2-rho), 2)')]);
disp(['二维高斯分布举荐-MH 的平均接受率为：', ...
    num2str(mean(accept_rate2, 2)')]);
