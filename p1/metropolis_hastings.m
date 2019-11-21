function [rho1, accept_rate1, rho2, accept_rate2] = ...
    metropolis_hastings(N, BURN_IN_FACTOR, mu, sigma, enable_output)
%% �������㡢��������
P = @(x) mvnpdf(x, mu, sigma);
n_skip = round(N * BURN_IN_FACTOR);
sigma1 = sqrt(sigma(1, 1));
sigma2 = sqrt(sigma(2, 2));
rho = sigma(1, 2) / (sigma1 * sigma2);

%% ������������
s0 = mvnrnd(mu, sigma, N)';

%% �ټ��ֲ�Ϊ���ȷֲ�
% ����ÿһά��ֻȡ n sigma
n = 4;
s1(1, :) = unifrnd(-sigma1, sigma1, 1, N) * n + mu(1);
s1(2, :) = unifrnd(-sigma2, sigma2, 1, N) * n + mu(2);
n_accept1 = 0;

% ���е���
for t = 1:N-1
    if ~isequal(s1(:, t), s1(:, t + 1))
        % ������ܸ��ʣ����оټ��ֲ�������Գ��Կɶ���
        alpha = min(P(s1(:, t + 1)) / P(s1(:, t)), 1);
        if rand() <= alpha
            n_accept1 = n_accept1 + 1;
        else
            % �����תʧ�ܣ�����һʱ�̵�״̬Ӧ��˿�״̬һ��
            s1(:, t + 1) = s1(:, t); 
        end
    end
end

% �������ϵ�������ܱ���
corr_s1 = corrcoef(s1(1, n_skip:end), s1(2, n_skip:end));
rho1 = corr_s1(1, 2);
accept_rate1 = n_accept1 / (N - 1);

%% �ټ��ֲ�Ϊ��ά��˹�ֲ�
s2 = zeros(2, N);
s2(:, 1) = mvnrnd(mu, sigma)';
n_accept2 = 0;

% ���е���
for t = 1:N-1
    Y = mvnrnd(s2(:, t), sigma)';
    if ~isequal(s2(:, t), Y)
        % ������ܸ��ʣ����оټ��ֲ�������Գ��Կɶ���
        alpha = min(P(Y) / P(s2(:, t)), 1);
        if rand() <= alpha
            n_accept2 = n_accept2 + 1;
            s2(:, t + 1) = Y;
        else
            % �����תʧ�ܣ�����һʱ�̵�״̬Ӧ��˿�״̬һ��
            s2(:, t + 1) = s2(:, t); 
        end
    end
end

% �������ϵ�������ܱ���
corr_s2 = corrcoef(s2(1, n_skip:end), s2(2, n_skip:end));
rho2 = corr_s2(1, 2);
accept_rate2 = n_accept2 / (N - 1);

%% ��ͼ�����
if enable_output
    disp(['�������ϵ���� ', num2str(rho)]);
    disp(['���ȷֲ��ټ�-MH �������ϵ���� ' , num2str(rho1)]);
    disp(['���ȷֲ��ټ�-MH ���ܱ�����', num2str(accept_rate1)]);
    disp(['��ά��˹�ֲ��ټ�-MH �������ϵ���� ' , num2str(rho2)]);
    disp(['��ά��˹�ֲ��ټ�-MH ���ܱ�����', num2str(accept_rate2)]);
    
    % ���ȷֲ��ټ�������ֵ��ͼ�Ա�
    figure;
    subplot(1, 2, 1);
    histogram2(s0(1, n_skip:end), s0(2, n_skip:end));
    title('����ֵ');
    subplot(1, 2, 2);
    histogram2(s1(1, n_skip:end), s1(2, n_skip:end));
    title('���ȷֲ��ټ�-MH');

    % ��ά��˹�ֲ��ټ�������ֵ��ͼ�Ա�
    figure;
    subplot(1, 2, 1);
    histogram2(s0(1, n_skip:end), s0(2, n_skip:end));
    title('����ֵ');
    subplot(1, 2, 2);
    histogram2(s2(1, n_skip:end), s2(2, n_skip:end));
    title('��ά��˹�ֲ��ټ�-MH ');
    
    % ����Ա�
    figure;
    subplot(1, 2, 1);
    hold on;
    title('���۷ֲ�����ȷֲ��ټ�-MH ����');
    plot(s0(1, n_skip:end), s0(2, n_skip:end), '.b');
    plot(s1(1, n_skip:end), s1(2, n_skip:end), '.r');
    subplot(1, 2, 2);
    hold on;
    title('���۷ֲ����ά��˹�ֲ��ټ�-MH ����');
    plot(s0(1, n_skip:end), s0(2, n_skip:end), '.b');
    plot(s2(1, n_skip:end), s2(2, n_skip:end), '.r');
end
end