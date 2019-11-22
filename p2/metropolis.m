function metropolis(K, N_Q, BETA, N, BURN_IN_FACTOR)
disp('Running Metropolis algorithm:');
tic;
n_skip = round(N * BURN_IN_FACTOR); % ͳ��ʱ������ǰ�������ĸ���
n_per_step = K ^ 2; % ����һ�����������������

x = randi([1, N_Q], K, K); % ��һ���������ȡֵ
u = zeros(N, 1); % ����
u(1) = calc_u(x);

% ���е���
for t = 2:N
    overall_delta_u = 0;
    for n = 1:n_per_step
        row = randi([1, K]);
        col = randi([1, K]);
        cur = x(row, col);
        candidate = randi([1, N_Q]);
        if cur ~= candidate
            % ���������ĸı�
            left = x(mod(row - 2, K) + 1, col);
            right = x(mod(row, K) + 1, col);
            up = x(row, mod(col - 2, K) + 1);
            down = x(row, mod(col, K) + 1);
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
            % alpha = min(exp(-BETA * delta_u), 1);
            if rand() <= exp(-BETA * delta_u)
                % �ɹ�����
                x(row, col) = candidate;
                overall_delta_u = overall_delta_u + delta_u;
            end
        end
    end
    u(t) = u(t-1) + overall_delta_u;
end

u_mean = mean(u(n_skip:end));
toc
disp(['E{u(x)} = ', num2str(u_mean)]);
% ln_z = LN_Z0 - u_mean * BETA;
% disp(['ln Z(T) = ', num2str(ln_z)]);
end
