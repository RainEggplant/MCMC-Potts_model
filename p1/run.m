N = 100000; % µü´ú´ÎÊý
BURN_IN_FACTOR = 0.3;
mu = [5; 10];
sigma = [1, 1; 1, 4];

metropolis_hastings(N, BURN_IN_FACTOR, mu, sigma, 1);
