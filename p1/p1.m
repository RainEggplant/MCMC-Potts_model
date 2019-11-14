clear;
%% ��������������
N = 100000; % ��������
SKIP_INITIAL_FACTOR = 0.2;
mu = [5; 10];
sigma = [1, 1; 1, 4];
P = @(x) mvnpdf(x, mu, sigma);

n_skip = round(N * SKIP_INITIAL_FACTOR);

%% ��������ֵ
s0 = mvnrnd(mu, sigma, N)';
% corr_s0 = corrcoef(s0(1, :), s0(2, :));
% disp(['�������ϵ���� ' ,num2str(corr_s0(1, 2))]);
disp('�������ϵ���� 0.5');

%% ����ټ��ֲ�Ϊ���ȷֲ�
% ����ÿһά��ֻȡ n sigma
n = 4;
sigma1 = sqrt(sigma(1, 1));
sigma2 = sqrt(sigma(2, 2));
s1(1, :) = unifrnd(-sigma1, sigma1, 1, N) * n + mu(1);
s1(2, :) = unifrnd(-sigma2, sigma2, 1, N) * n + mu(2);

% ���е���
for t = 1:N-1
    if ~isequal(s1(:, t), s1(:, t + 1))
        % ������ܸ���
        alpha = min(P(s1(:, t + 1)) / P(s1(:, t)), 1);
        if rand() >= alpha
            % �����תʧ�ܣ�����һʱ�̵�״̬Ӧ��˿�״̬һ��
           s1(:, t + 1) = s1(:, t); 
        end
    end
end

% ������ֵ��ͼ�Ա�
figure;
subplot(1, 2, 1);
histogram2(s0(1, n_skip:end), s0(2, n_skip:end));
title('����ֵ');
subplot(1, 2, 2);
histogram2(s1(1, n_skip:end), s1(2, n_skip:end));
title('���ȷֲ��ټ�-MH');
corr_s1 = corrcoef(s1(1, n_skip:end), s1(2, n_skip:end));
disp(['���ȷֲ��ټ�-MH �������ϵ���� ' ,num2str(corr_s1(1, 2))]);

%% ����ټ��ֲ�Ϊ��ά��˹�ֲ�
% 
% Q = @(x, mu) mvnpdf(x, mu, [0.5, 0; 0, 2]);
Q = @(x, mu) mvnpdf(x, mu, sigma);
s2(1, :) = unifrnd(-sigma1, sigma1, 1, N) * n + mu(1);
s2(2, :) = unifrnd(-sigma2, sigma2, 1, N) * n + mu(2);

% ���е���
for t = 1:N-1
    if ~isequal(s2(:, t), s2(:, t + 1))
        % ������ܸ���
        Q_forward = Q(s2(:, t + 1), s2(:, t));
        Q_backward = Q(s2(:, t), s2(:, t + 1));
        alpha = min( ...
            P(s2(:, t + 1)) * Q_backward / (P(s2(:, t)) * Q_forward), ...
            1);
        if rand() >= alpha
            % �����תʧ�ܣ�����һʱ�̵�״̬Ӧ��˿�״̬һ��
           s2(:, t + 1) = s2(:, t); 
        end
    end
end

% ������ֵ��ͼ�Ա�
figure;
subplot(1, 2, 1);
histogram2(s0(1, n_skip:end), s0(2, n_skip:end));
title('����ֵ');
subplot(1, 2, 2);
histogram2(s2(1, n_skip:end), s2(2, n_skip:end));
title('��ά��˹�ֲ��ټ�-MH ');
corr_s2 = corrcoef(s2(1, n_skip:end), s2(2, n_skip:end));
disp(['��ά��˹�ֲ��ټ�-MH �������ϵ���� ' ,num2str(corr_s2(1, 2))]);

%% ����Ա�
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
