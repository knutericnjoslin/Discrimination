% DISCRIMINATION: GLOBAL DEVIATION FROM SEPARATING
% ================================================
% Profitable deviation by firm setting w_h to w_d below w_l 

%%% Initialize values 

% Function parameters
i=1;
j=1; % REMOVED
% Market parameters
N=3;
B=3;
x=2;
y=2;
% Wages in original equilibrium
w_h = 0.34;
w_l = w_h*(1 - (1 - 1/x)^N ) / ( N / x); 


% White worker behavior
theta_l = @(theta_h, x, y) (1 - (x-1)*theta_h)/y ; % White application to low firm written as function of theta_h

theta_h_optimal = fsolve( @(theta_h) w_l*binom_sum_constructor(N, i, theta_l(theta_h, x, y))...
    - w_h*binom_sum_constructor(N, i, theta_h), 1 );

% Deviant wage
w_d = w_h*(1 - (1 - theta_h_optimal)^N)/(N*theta_h_optimal);


% Black worker behavior
gamma_l = @(gamma_d, y) (1 - gamma_d)/y;

gamma_d_optimal = fsolve( @(gamma_d) w_l*((1 - (1 - (x-1)*theta_h_optimal)/y)^N) *binom_sum_constructor(B, i, gamma_l(gamma_d, y))...
    - w_d*binom_sum_constructor(B, i, gamma_d), 0.5 );


% Profit original equilibrium

pi_o = (1 - w_h) * ( 1 - ( 1 - 1/x)^N );


% Profit deviant

pi_d = (1 - w_d) * ( 1 - ( 1 - gamma_d_optimal)^B );


