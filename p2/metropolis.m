function metropolis(K, N_Q, BETA, N, BURN_IN_FACTOR)
disp('Running Metropolis algorithm:');
tic;
n_skip = round(N * BURN_IN_FACTOR); % ͳ��ʱ������ǰ�������ĸ���
x = zeros(K, K, N); % ����
x(:, :, 1) = randi([1, N_Q], K, K); % ��һ���������ȡֵ
n_per_step = K ^ 2; % ����һ�����������������

% ���е���
for t = 2:N
    x(:, :, t) = x(:, :, t - 1);
    for n = 1:n_per_step
        row = randi([1, K]);
        col = randi([1, K]);
        cur = x(row, col, t);
        candidate = randi([1, N_Q]);
        if cur ~= candidate
            % ���������ĸı�
            left = x(mod(row - 2, K) + 1, col, t);
            right = x(mod(row, K) + 1, col, t);
            up = x(row, mod(col - 2, K) + 1, t);
            down = x(row, mod(col, K) + 1, t);
            delta_u = 0;
            % ����������ͬ��
            delta_u = delta_u + ...
                (cur == left) + (cur == right) + ...
                (cur == up) + (cur == down);
            % ����������ͬ��
            delta_u = delta_u - ...
                (left == candidate) - (right == candidate) - ...
                (up == candidate) - (down == candidate);
            
            % ������ܸ���
            alpha = min(exp(-BETA * delta_u), 1);
            if rand() <= alpha
                % �ɹ�����
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