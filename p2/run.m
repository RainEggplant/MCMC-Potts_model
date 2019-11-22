clear;
%% ��������������
N = 2200000; % ��������
N_BURN_IN = 200000;
SUBSAMPLE_STEP = 10;
BURN_IN_FACTOR = 0.3;
K = 20; % 20*20 ����
N_Q = 10; % 10 ��ȡֵ
BETA = 1.4; % 1/T ��ֵ
BETA_ARRAY = [1.4; 1.4065; 1.413; 1.4195; 1.426];
LN_Z0 = 400 * log(10);

% metropolis(K, N_Q, BETA, N, BURN_IN_FACTOR);
% swendsen_wang(K, N_Q, BETA, N, BURN_IN_FACTOR);
samc(K, N_Q, BETA_ARRAY, N, N_BURN_IN, SUBSAMPLE_STEP);
