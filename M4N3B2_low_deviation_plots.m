% ### Examining shape of dtheta_d__dw_d for various values of w_d

i=1;

% Market parameters
N=3;
B=2;
x=2;
y=2;
% Wages in original equilibrium
w_h = 0.3651;
w_l = w_h*(1 - (1 - 1/x)^N ) / ( N / x); 

w_d = (w_l+0.01):0.001:(w_h-0.01);


% Consider a deviation by one of the firms setting w_l upward to some 
% wage w_d such that w_l < w_d < w_h. 


% ### White application behavior - Never apply to firms offering w_l

% theta_h as a function of theta_d
theta_h = @(theta_d, x) (1 - theta_d)/x ;

% ## Solving for theta_d given market constellation using equation 3.3

theta_d_optimal = zeros(1,133);
for j = 1:133;
    theta_d_optimal(j) = fsolve( @(theta_d) w_h*binom_sum_constructor(N, i, theta_h(theta_d, x)) - ...
        w_d(j) * binom_sum_constructor(N, i, theta_d) , 0.25 );
end;

% # Plot of the application probabilities of whites

plot(w_d, theta_d_optimal)
hold on






% ### Black application behavior - 
% gamma_l = @(gamma_d, y) (1 - gamma_d)/y;
% 
% gamma_d_optimal = zeros(1,133);
% for j = 1:133;
%    gamma_d_optimal(j) = lsqnonlin( @(gamma_d) w_d(j)*( 1 - theta_d_optimal(j) )^N * binom_sum_constructor(B, i, gamma_d)...
%     - (w_l/(y-1))*binom_sum_constructor(B, i, gamma_l(gamma_d,y)), 0.5,0,1) ;
% end;

% Black app probs go to zero beyond w(60)=0.282
gamma_d_optimal = zeros(1,60);
for j = 1:133;
   gamma_d_optimal(j) = lsqnonlin( @(gamma_d) w_d(j)*( 1 - theta_d_optimal(j) )^N * binom_sum_constructor(B, i, gamma_d)...
    - (w_l/(y-1))*binom_sum_constructor(B, i, gamma_l(gamma_d,y)), 0.5,0,1) ;
end;

% # Plot of the application probabilities of blacks
plot(w_d, gamma_d_optimal)




% Profits
% -------

pi_d = (1 - w_d).*(1 - (1 - theta_d_optimal).^N .* (1 - gamma_d_optimal).^B );

% Plot of profits
plot(w_d, pi_d)
hold on



% ### derivative of profits

% dtheta_d__dw_d
dtheta_d__dw_d = binom_sum_constructor( N, i, theta_d_optimal ) ./ ( w_d.*deriv_binom_sum_constructor( N, i, theta_d_optimal ) + (w_h/x)*deriv_binom_sum_constructor( N, i, theta_h(theta_d_optimal, x)) );   

% dgamma_d__dw_d
dgamma_d__dw_d = ( (1 - theta_d_optimal).^N .* binom_sum_constructor(B, i, gamma_d_optimal) - w_d.*N.*(1 - theta_d_optimal).^(N-1).*binom_sum_constructor(B, i, gamma_d_optimal).*dtheta_d__dw_d  )...
    ./ ( w_d.*( 1 - theta_d_optimal ).^N .* deriv_binom_sum_constructor(B, i, gamma_d_optimal) + (w_l/(y-1)).*deriv_binom_sum_constructor(B, i, gamma_l(gamma_d_optimal, y) ) ) ;
% Notice that beyond w(60) = 0.282 no blacks apply, so changes in this
% region do not affect the application probabilities
dgamma_d__dw_d(1, 59:133) = 0;


d_pi__d_wd = - (1 - (1 - theta_d_optimal).^N .* (1 - gamma_d_optimal).^B ) + ...
    (1 - w_d).*N.*(1 - theta_d_optimal).^(N-1).*(1 - gamma_d_optimal).^B .* dtheta_d__dw_d + ...
    (1 - w_d).*(1 - theta_d_optimal).^N .*B.* (1 - gamma_d_optimal).^(B-1) .* dgamma_d__dw_d ;

plot(w_d, d_pi__d_wd)

