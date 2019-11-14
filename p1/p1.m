clear;
%% 常数、参数定义
N = 100000; % 迭代次数
SKIP_INITIAL_FACTOR = 0.2;
mu = [5; 10];
sigma = [1, 1; 1, 4];
P = @(x) mvnpdf(x, mu, sigma);

n_skip = round(N * SKIP_INITIAL_FACTOR);

%% 计算理论值
s0 = mvnrnd(mu, sigma, N)';
% corr_s0 = corrcoef(s0(1, :), s0(2, :));
% disp(['理论相关系数： ' ,num2str(corr_s0(1, 2))]);
disp('理论相关系数： 0.5');

%% 计算举荐分布为均匀分布
% 这里每一维上只取 n sigma
n = 4;
sigma1 = sqrt(sigma(1, 1));
sigma2 = sqrt(sigma(2, 2));
s1(1, :) = unifrnd(-sigma1, sigma1, 1, N) * n + mu(1);
s1(2, :) = unifrnd(-sigma2, sigma2, 1, N) * n + mu(2);

% 进行迭代
for t = 1:N-1
    if ~isequal(s1(:, t), s1(:, t + 1))
        % 计算接受概率
        alpha = min(P(s1(:, t + 1)) / P(s1(:, t)), 1);
        if rand() >= alpha
            % 如果跳转失败，则下一时刻的状态应与此刻状态一致
           s1(:, t + 1) = s1(:, t); 
        end
    end
end

% 与理论值作图对比
figure;
subplot(1, 2, 1);
histogram2(s0(1, n_skip:end), s0(2, n_skip:end));
title('理论值');
subplot(1, 2, 2);
histogram2(s1(1, n_skip:end), s1(2, n_skip:end));
title('均匀分布举荐-MH');
corr_s1 = corrcoef(s1(1, n_skip:end), s1(2, n_skip:end));
disp(['均匀分布举荐-MH 估计相关系数： ' ,num2str(corr_s1(1, 2))]);

%% 计算举荐分布为二维高斯分布
% 
% Q = @(x, mu) mvnpdf(x, mu, [0.5, 0; 0, 2]);
Q = @(x, mu) mvnpdf(x, mu, sigma);
s2(1, :) = unifrnd(-sigma1, sigma1, 1, N) * n + mu(1);
s2(2, :) = unifrnd(-sigma2, sigma2, 1, N) * n + mu(2);

% 进行迭代
for t = 1:N-1
    if ~isequal(s2(:, t), s2(:, t + 1))
        % 计算接受概率
        Q_forward = Q(s2(:, t + 1), s2(:, t));
        Q_backward = Q(s2(:, t), s2(:, t + 1));
        alpha = min( ...
            P(s2(:, t + 1)) * Q_backward / (P(s2(:, t)) * Q_forward), ...
            1);
        if rand() >= alpha
            % 如果跳转失败，则下一时刻的状态应与此刻状态一致
           s2(:, t + 1) = s2(:, t); 
        end
    end
end

% 与理论值作图对比
figure;
subplot(1, 2, 1);
histogram2(s0(1, n_skip:end), s0(2, n_skip:end));
title('理论值');
subplot(1, 2, 2);
histogram2(s2(1, n_skip:end), s2(2, n_skip:end));
title('二维高斯分布举荐-MH ');
corr_s2 = corrcoef(s2(1, n_skip:end), s2(2, n_skip:end));
disp(['二维高斯分布举荐-MH 估计相关系数： ' ,num2str(corr_s2(1, 2))]);

%% 整体对比
figure;
subplot(1, 2, 1);
hold on;
title('理论分布与均匀分布举荐-MH 采样');
plot(s0(1, n_skip:end), s0(2, n_skip:end), '.b');
plot(s1(1, n_skip:end), s1(2, n_skip:end), '.r');
subplot(1, 2, 2);
hold on;
title('理论分布与二维高斯分布举荐-MH 采样');
plot(s0(1, n_skip:end), s0(2, n_skip:end), '.b');
plot(s2(1, n_skip:end), s2(2, n_skip:end), '.r');
