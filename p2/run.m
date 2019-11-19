clear;
%% 常数、参数定义
N = 10000; % 迭代次数
BURN_IN_FACTOR = 0.3;
K = 20; % 20*20 网格
N_Q = 10; % 10 种取值
BETA = 1.4; % 1/T 的值
% BETA = [1.4; 1.4065; 1.413; 1.4195; 1.426];
LN_Z0 = 400 * log(10);

metropolis(K, N_Q, BETA, N, BURN_IN_FACTOR);
swedensen_wang(K, N_Q, BETA, N, BURN_IN_FACTOR);
