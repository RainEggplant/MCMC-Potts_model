function metropolis(K, N_Q, BETA, N, BURN_IN_FACTOR)
disp('Running Metropolis algorithm:');
tic;
n_skip = round(N * BURN_IN_FACTOR); % 统计时跳过最前面样本的个数
x = zeros(K, K, N); % 样本
x(:, :, 1) = randi([1, N_Q], K, K); % 第一个样本随机取值
n_per_step = K ^ 2; % 生成一个样本所需迭代次数

% 进行迭代
for t = 2:N
    x(:, :, t) = x(:, :, t - 1);
    for n = 1:n_per_step
        row = randi([1, K]);
        col = randi([1, K]);
        cur = x(row, col, t);
        candidate = randi([1, N_Q]);
        if cur ~= candidate
            % 计算能量的改变
            left = x(mod(row - 2, K) + 1, col, t);
            right = x(mod(row, K) + 1, col, t);
            up = x(row, mod(col - 2, K) + 1, t);
            down = x(row, mod(col, K) + 1, t);
            delta_u = 0;
            % 少了数对相同的
            delta_u = delta_u + ...
                (cur == left) + (cur == right) + ...
                (cur == up) + (cur == down);
            % 多了数对相同的
            delta_u = delta_u - ...
                (left == candidate) - (right == candidate) - ...
                (up == candidate) - (down == candidate);
            
            % 计算接受概率
            alpha = min(exp(-BETA * delta_u), 1);
            if rand() <= alpha
                % 成功接受
                x(row, col, t) = candidate;
            end
        end
    end
end

u = calc_u(x(:, :, n_skip:N));
u_mean = mean(u);
toc
disp(['E{u(x)} = ', num2str(u_mean)]);
% ln_z = LN_Z0 - u_mean * BETA;
% disp(['ln Z(T) = ', num2str(ln_z)]);
end