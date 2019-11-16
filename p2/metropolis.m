clear;
%% 常数、参数定义
N = 100000; % 迭代次数
SKIP_INITIAL_FACTOR = 0.3;
n_skip = round(N * SKIP_INITIAL_FACTOR);

K = 20; % 20*20 网格
q = 10; % 10 种取值
ONE_DIV_T = 1.426; % 1/T 的值
% ONE_DIV_T = [1.4; 1.4065; 1.413; 1.4195; 1.426];
LN_Z0 = 400 * log(10);

%%
x = zeros(K, K, N);
x(:, :, 1) = randi([1, q], K, K);
n_per_step = K ^ 2;

% 进行迭代
for t = 2:N
    x(:, :, t) = x(:, :, t - 1);
    for n = 1:n_per_step
        row = randi([1, K]);
        col = randi([1, K]);
        cur = x(row, col, t);
        candidate = randi([1, q]);
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
            alpha = min(exp(-ONE_DIV_T * delta_u), 1);
            if rand() <= alpha
                % 成功接受
                x(row, col, t) = candidate;
            end
        end
    end
end

u = calc_u(x(:, :, n_skip:N));
u_mean = mean(u);
disp(['E{u(x)} = ', num2str(u_mean)]);
ln_z = LN_Z0 - u_mean * ONE_DIV_T;
disp(['ln Z(T) = ', num2str(ln_z)]);
