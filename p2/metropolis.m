function [u, x] = metropolis(K, N_Q, BETA, N)
n_per_step = K ^ 2; % 生成一个样本所需迭代次数

x = randi([1, N_Q], K, K); % 第一个样本随机取值
u = zeros(N, 1); % 能量
u(1) = calc_u(x);

% 进行迭代
for t = 2:N
    overall_delta_u = 0;
    for n = 1:n_per_step
        row = randi([1, K]);
        col = randi([1, K]);
        cur = x(row, col);
        candidate = randi([1, N_Q]);
        if cur ~= candidate
            % 计算能量的改变
            left = x(mod(row - 2, K) + 1, col);
            right = x(mod(row, K) + 1, col);
            up = x(row, mod(col - 2, K) + 1);
            down = x(row, mod(col, K) + 1);
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
            % alpha = min(exp(-BETA * delta_u), 1);
            if rand() <= exp(-BETA * delta_u)
                % 成功接受
                x(row, col) = candidate;
                overall_delta_u = overall_delta_u + delta_u;
            end
        end
    end
    u(t) = u(t-1) + overall_delta_u;
end
end
