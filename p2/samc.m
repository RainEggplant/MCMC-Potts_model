function samc(K, N_Q, BETA, N, N_BURN_IN, SUBSAMPLE_STEP)
disp('Running Self-adjusted mixture sampling algorithm:');
tic;
N_SUBSAMPLE = ceil(N / SUBSAMPLE_STEP);
N_PER_STEP = K ^ 2; % 生成一个样本所需迭代次数
BETA_FACTOR = 0.8;
N_BETA = 5; % Hard-coded! BETA 的大小，即温度倒数的个数
DIAG_ELEM_BIAS = -N_BURN_IN + N_BURN_IN ^ BETA_FACTOR;
PI_ANY = 1 / N_BETA;

label = 1; % 温度的标号, 初始化为 1

% size of neighboring labels (i.e. N(k))
% also, the proposal distribution 
%   Gamma(label_1, label_2) = 1 / S(label_1)
% WARNING! This is hard-coded!
S = [1, 2, 2, 2, 1]'; 

zeta = zeros(N_BETA, 1); % 对数相对归一化常数，初始化为 0
zeta_rec = zeros(N_BETA, N_SUBSAMPLE);
x = randi([1, N_Q], K, K); % 第一个样本随机取值
u = calc_u(x); % 能量
u_index = zeros(N_BETA, 1);
u_rec = zeros(N_BETA, N_SUBSAMPLE);


for t = 1:N
    % jump
    label_next = random_label(label);
    % 因为 forall k, pi_k = 1/N_BETA, 故该项可抵消
    p_cur = exp(-BETA(label) * u - zeta(label));
    p_next = exp(-BETA(label_next) * u - zeta(label_next));
    alpha = S(label) / S(label_next) * p_next / p_cur;
    if rand() < alpha
        label = label_next;
    end
    
    % move   
    overall_delta_u = 0;
    for n = 1:N_PER_STEP
        row = randi([1, K]);
        col = randi([1, K]);
        cur = x(row, col);
        
        % 从除去当前取值的剩余 N_Q - 1 个取值中随机选取
        candidate = randi([1, N_Q - 1]);
        if (cur <= candidate)
            candidate = candidate + 1;
        end
        
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

        if rand() <= exp(-BETA(label) * delta_u)
            % 成功接受
            x(row, col) = candidate;
            overall_delta_u = overall_delta_u + delta_u;
        end
    end
    
    u = u + overall_delta_u;

    
    % update
    binary_update();
    
    % sub-sample
    if mod(t - 1, SUBSAMPLE_STEP) == 0
        index = u_index(label) + 1;
        u_index(label) = index;
        u_rec(label, index) = u;
        zeta_rec(:, (t - 1) / 10 + 1) = zeta;
    end
        
end

toc;
disp(['zeta = ', num2str(zeta(:)')]);
figure;
hold on;
for l = 1:N_BETA
    ksdensity(u_rec(l, ceil(N_BURN_IN / SUBSAMPLE_STEP / N_BETA):u_index(l)) / K ^ 2);
end
legend('1', '2', '3', '4', '5');
disp('End running.')


    function rnd_label = random_label(label)
        if label == 1
            rnd_label = 2;
        elseif label == N_BETA
            rnd_label = N_BETA - 1;
        else
            if rand() < 0.5
                rnd_label = label - 1;
            else
                rnd_label = label + 1;
            end
        end
    end

    function binary_update()
        if t <= N_BURN_IN
            diag_elem = min(PI_ANY, t ^ (-BETA_FACTOR));
        else
            diag_elem = min(PI_ANY, 1 / (t + DIAG_ELEM_BIAS));
        end
        
        if label == 1
            zeta(1) = zeta(1) + diag_elem / PI_ANY;
            for k = 2:N_BETA
                zeta(k) = zeta(k) - zeta(1);
            end
            zeta(1) = 0;
        else
            zeta(label) = zeta(label) + diag_elem / PI_ANY;
        end
    end
end
