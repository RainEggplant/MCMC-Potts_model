function [rho1, accept_rate1, rho2, accept_rate2] = ...
    metropolis_hastings(N, BURN_IN_FACTOR, mu, sigma, enable_output)
%% 参数计算、函数定义
P = @(x) mvnpdf(x, mu, sigma);
n_skip = round(N * BURN_IN_FACTOR);
sigma1 = sqrt(sigma(1, 1));
sigma2 = sqrt(sigma(2, 2));
rho = sigma(1, 2) / (sigma1 * sigma2);

%% 生成理论样本
s0 = mvnrnd(mu, sigma, N)';

%% 举荐分布为均匀分布
% 这里每一维上只取 n sigma
n = 4;
s1(1, :) = unifrnd(-sigma1, sigma1, 1, N) * n + mu(1);
s1(2, :) = unifrnd(-sigma2, sigma2, 1, N) * n + mu(2);
n_accept1 = 0;

% 进行迭代
for t = 1:N-1
    if ~isequal(s1(:, t), s1(:, t + 1))
        % 计算接受概率，其中举荐分布项满足对称性可对消
        alpha = min(P(s1(:, t + 1)) / P(s1(:, t)), 1);
        if rand() <= alpha
            n_accept1 = n_accept1 + 1;
        else
            % 如果跳转失败，则下一时刻的状态应与此刻状态一致
            s1(:, t + 1) = s1(:, t); 
        end
    end
end

% 计算相关系数、接受比例
corr_s1 = corrcoef(s1(1, n_skip:end), s1(2, n_skip:end));
rho1 = corr_s1(1, 2);
accept_rate1 = n_accept1 / (N - 1);

%% 举荐分布为二维高斯分布
s2 = zeros(2, N);
s2(:, 1) = mvnrnd(mu, sigma)';
n_accept2 = 0;

% 进行迭代
for t = 1:N-1
    Y = mvnrnd(s2(:, t), sigma)';
    if ~isequal(s2(:, t), Y)
        % 计算接受概率，其中举荐分布项满足对称性可对消
        alpha = min(P(Y) / P(s2(:, t)), 1);
        if rand() <= alpha
            n_accept2 = n_accept2 + 1;
            s2(:, t + 1) = Y;
        else
            % 如果跳转失败，则下一时刻的状态应与此刻状态一致
            s2(:, t + 1) = s2(:, t); 
        end
    end
end

% 计算相关系数、接受比例
corr_s2 = corrcoef(s2(1, n_skip:end), s2(2, n_skip:end));
rho2 = corr_s2(1, 2);
accept_rate2 = n_accept2 / (N - 1);

%% 作图与输出
if enable_output
    disp(['理论相关系数： ', num2str(rho)]);
    disp(['均匀分布举荐-MH 估计相关系数： ' , num2str(rho1)]);
    disp(['均匀分布举荐-MH 接受比例：', num2str(accept_rate1)]);
    disp(['二维高斯分布举荐-MH 估计相关系数： ' , num2str(rho2)]);
    disp(['二维高斯分布举荐-MH 接受比例：', num2str(accept_rate2)]);
    
    % 均匀分布举荐与理论值作图对比
    figure;
    subplot(1, 2, 1);
    histogram2(s0(1, n_skip:end), s0(2, n_skip:end));
    title('理论值');
    subplot(1, 2, 2);
    histogram2(s1(1, n_skip:end), s1(2, n_skip:end));
    title('均匀分布举荐-MH');

    % 二维高斯分布举荐与理论值作图对比
    figure;
    subplot(1, 2, 1);
    histogram2(s0(1, n_skip:end), s0(2, n_skip:end));
    title('理论值');
    subplot(1, 2, 2);
    histogram2(s2(1, n_skip:end), s2(2, n_skip:end));
    title('二维高斯分布举荐-MH ');
    
    % 整体对比
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
end
end