% DISCRIMINATION: LOCAL DEVIATION FROM SEPARATING BY FIRM SETTING w_l
% ===================================================================
% Profitable deviation by firm setting w_h to w_d below w_l 

%%% Initialize values 

% Function parameters
i=1;
% j=1;
% Market parameters
N=3;
B=2;
x=2;
y=2;
% Wages in original equilibrium
w_h = 0.3651;
w_l = w_h*(1 - (1 - 1/x)^N ) / ( N / x); 


% White worker behavior
theta_l = @(theta_h, x, y) (1 - (x-1)*theta_h)/y ; % White application to low firm written as function of theta_h

theta_h_optimal = fsolve( @(theta_h) w_l*binom_sum_constructor(N, i, theta_l(theta_h, x, y))...
    - w_h*binom_sum_constructor(N, i, theta_h), 0.5 );

% Deviant wage
w_d = w_h*(1 - (1 - theta_h_optimal)^N)/(N*theta_h_optimal);


% Black worker behavior
gamma_l = @(gamma_d, y) (1 - gamma_d)/y;

gamma_d_optimal = fsolve( @(gamma_d) w_l*((1 - (1 - (x-1)*theta_h_optimal)/y)^N) *binom_sum_constructor(B, i, gamma_l(gamma_d, y))...
    - w_d*binom_sum_constructor(B, i, gamma_d), 0.5 );








% LOCAL OPTIMALITY LOW DEVIATION 
% ------------------------------

% *Original equilibrium, with x firm setting w_h and y firms setting w_l
% *N white workers and B black workers

w_h = 0.3651 ;
w_l = w_h*(1 - ( 1 - 1/x)^N )/(N/x) ;
w_d = w_l*1.001

% Consider a deviation by one of the firms setting w_l upward to some 
% wage w_d such that w_l < w_d < w_h. 


% ### White application behavior - Never apply to firms offering w_l

% theta_h as a function of theta_d
theta_h = @(theta_d, x) (1 - theta_d)/x ;


% ## Solving for theta_d given market constellation using equation 3.3

theta_d_optimal = fsolve( @(theta_d) w_h*binom_sum_constructor(N, i, theta_h(theta_d, x))...
    - w_d*binom_sum_constructor(N, i, theta_d), 0.5 );

% ## dtheta_d__dw_d
dtheta_d__dw_d = binom_sum_constructor( N, i, theta_d_optimal ) / ( w_d*deriv_binom_sum_constructor( N, i, theta_d_optimal ) + (w_h/x)*deriv_binom_sum_constructor( N, i, theta_h(theta_d_optimal, x)) );   













% ### Black application behavior - "Never" apply to firms offering w_h 


% gamma_l as a function of gamma_d

gamma_l = @(gamma_d, y) ( 1 - gamma_d )/(y - 1) ;


% ## Solving for gamma_d given market constellation using equation 3.6

gamma_d_optimal = fsolve( @(gamma_d) w_d*( 1 - theta_d_optimal )^N * binom_sum_constructor(B, i, gamma_d)...
    - (w_l)*binom_sum_constructor(B, i, gamma_l(gamma_d,y)), 0.5) ;

% ## dgamma_d__dw_d

dgamma_d__dw_d = ( (1 - theta_d_optimal)^N * binom_sum_constructor(B, i, gamma_d_optimal) - w_d*N*(1 - theta_d_optimal)^(N-1)*binom_sum_constructor(B, i, gamma_d_optimal)*dtheta_d__dw_d  )/...
    ( w_d*( 1 - theta_d_optimal )^N * deriv_binom_sum_constructor(B, i, gamma_d_optimal) + (w_l/(y-1))*deriv_binom_sum_constructor(B, i, gamma_l(gamma_d_optimal, y) ) ) ;





% ### Profits

% Original equilibrium

pi_l = (1 - w_l)*(1 - (1 - 1/y)^B );

% After deviation

pi_d = (1 - w_d)*(1 - (1 - theta_d_optimal)^N * (1 - gamma_d_optimal)^B );


% ### Evaluate the derivative of the profits just above w_l: This
%     derivative is negative

d_pi_d_wd = - (1 - (1 - theta_d_optimal)^N * (1 - gamma_d_optimal)^B ) + ...
    (1 - w_d)*N*(1 - theta_d_optimal)^(N-1)*(1 - gamma_d_optimal)^B * dtheta_d__dw_d + ...
    (1 - w_d)*(1 - theta_d_optimal)^N *B* (1 - gamma_d_optimal)^(B-1) * dgamma_d__dw_d ;









% ### Evaluate the derivative of the profits just at w_l: This derivative
%     is positive 2.0215

w_d = w_l
theta_d_optimal = 0 % Since no whites apply when w_d = w_l
gamma_d_optimal = 1/y % Since this is separating, blacks appyl to all y firms setting w_l with equal probability

dtheta_d__dw_d = 0

dgamma_d__dw_d = ( (1 - theta_d_optimal)^N * binom_sum_constructor(B, i, gamma_d_optimal) - w_d*N*(1 - theta_d_optimal)^(N-1)*binom_sum_constructor(B, i, gamma_d_optimal)*dtheta_d__dw_d  )/...
    ( w_d*( 1 - theta_d_optimal )^N * deriv_binom_sum_constructor(B, i, gamma_d_optimal) + (w_l/(y-1))*deriv_binom_sum_constructor(B, i, gamma_l(gamma_d_optimal, y) ) ) ;


% Profits 

d_pi_d_wd_initial = - (1 - (1 - theta_d_optimal)^N * (1 - gamma_d_optimal)^B ) + ...
    (1 - w_d)*N*(1 - theta_d_optimal)^(N-1)*(1 - gamma_d_optimal)^B * dtheta_d__dw_d + ...
    (1 - w_d)*(1 - theta_d_optimal)^N *B* (1 - gamma_d_optimal)^(B-1) * dgamma_d__dw_d ;























% LOCAL OPTIMALITY HIGH DEVIATION 
% ------------------------------

% *Original equilibrium, with x firm setting w_h and y firms setting w_l
% *N white workers and B black workers

% Function parameters
i=1;

% Market parameters
N=3;
B=2;
x=2;
y=2;

w_h = 0.3651 ;
w_l = w_h*(1 - ( 1 - 1/x)^N )/(N/x) ;
w_d = w_h*0.99


% Consider a deviation by one of the firms setting w_h downward to some 
% wage w_d such that w_l < w_d < w_h. Hence x-1 firms setting w_h, y firms 
% setting w_l, and one firm setting w_d.




% ### White application behavior - Whites apply to all firms

% ## Stage 1. Solve for theta_h as a function of theta_d

% Eliminate theta_l by writing as a function of theta_d and theta_h
theta_l = @(theta_d, theta_h, x, y ) ( 1 - theta_d - (x-1)*theta_h )/y ;


% % ## EXAMPLE: Solving for theta_d and theta_h given market constellation using equation 3.3
% 
% [theta_optimal] = lsqnonlin( @(theta) [w_d*binom_sum_constructor(N, i, theta(1))...
%     - w_l*binom_sum_constructor(N, i, theta_l(theta(1), theta(2), x, y)); ...
%     w_d*binom_sum_constructor(N, i, theta(1))...
%     - w_h*binom_sum_constructor(N, i, theta(2))]...
%     , [0.5 0.5] , [0 0], [1 1] );


% ## Solving for theta_d and theta_h as w_d varies from w_l to w_h

w_d = (w_l+0.001):0.001:(w_h-0.001) ; % 1x151

theta_optimal = zeros( 151, 3); % Container with dimensions 151x50

for l=1:151;
    theta_optimal(l, 1:2) = lsqnonlin( @(theta) [w_d(l)*binom_sum_constructor(N, i, theta(1))...
    - w_l*binom_sum_constructor(N, i, theta_l(theta(1), theta(2), x, y)); ...
    w_d(l)*binom_sum_constructor(N, i, theta(1))...
    - w_h*binom_sum_constructor(N, i, theta(2))]...
    , [0.5 0.5] , [0 0], [1 1] );
end;

theta_optimal( :, 3) = theta_l(theta_optimal(:,1), theta_optimal(:,2), x, y);

plot(w_d, theta_optimal(:,1))
hold on
plot(w_d, theta_optimal(:,2))
plot(w_d, theta_optimal(:,3))

% ### Approximation of theta_d as a function of w_d
% (given w_l and w_h)

% ## Estimate the polynomial approximation for theta_d as a function of w_D
p = polyfit(transpose(w_d), theta_optimal(:,1), 3) % Polynomial approximation
p_fit = polyval(p, w_d) % Estimate values of theta_d from a polynomial approximation 

% compare the data and approximation
plot(w_d, p_fit)
hold on
plot(w_d, theta_optimal(:,1))

diff_approx = theta_optimal(:,1)-transpose(p_fit); % Difference in errors between simulated values and the approximation
plot(w_d, diff_approx)



% % Checking the approximation - Appears to be quite good
% test = @(w_d) p(1)*w_d.^3 + p(2)*w_d.^2 + p(3)*w_d + p(4)
% test_eval = test(w_d)
% 
% plot(w_d, test_eval)
% hold on
% plot(w_d, theta_optimal(:,1))

% % ## Approximation of theta_h as a function of w_d
% p_theta_h = polyfit(transpose(w_d), theta_optimal(:,2), 3); % Polynomial approximation
% p_fit_theta_h = polyval(p_theta_h, w_d); % Estimate values of theta_d from a polynomial approximation 
% 
% plot(w_d, p_fit_theta_h)
% hold on
% plot(w_d, theta_optimal(:,2))
% 
% 
% % Difference in errors between simulated values and the approximation
% diff_approx_theta_h = theta_optimal(:,2)-transpose(p_fit_theta_h);
% plot(w_d, diff_approx_theta_h)

% Derivative dtheta_d__dwd approximation based on polyfit
dtheta_d__dw_d = @(w_d) (3*p(1)).*w_d.^2 + (2*p(2)).*w_d + p(3);

dthetad_dwd = dtheta_d__dw_d(w_d);

plot(w_d, dthetad_dwd)


% ### Black application behavior 

% Eliminate gamma_l by writing as a function of gamma_d and gamma_h
gamma_l = @(gamma_d, gamma_h, x, y ) ( 1 - gamma_d - (x-1)*gamma_h )/y ;

% % ## Solving for theta_d and theta_h as w_d varies from w_l to w_h
% 
% gamma_optimal = zeros( 151, 3); % Container with dimensions 151x50
% 
% for l=1:151;
%     gamma_optimal(l, 1:2) = lsqnonlin( @(gamma) [w_d(l)*(1 - theta_optimal(l, 1))^N*binom_sum_constructor(B, i, gamma(1))...
%     - w_l*(1 - theta_optimal(l, 3))^N*binom_sum_constructor(B, i, gamma_l(gamma(1), gamma(2), x, y)); ...
%     w_d(l)*(1 - theta_optimal(l, 1))^N*binom_sum_constructor(B, i, gamma(1))...
%     - w_h*(1 - theta_optimal(l, 2))^N*binom_sum_constructor(B, i, gamma(2))]...
%     , [0.5 0.5] , [0 0], [1 1] );
% end;
% 
% gamma_optimal( :, 3) = gamma_l(gamma_optimal(:,1), gamma_optimal(:,2), x, y)




% ALTERNATIVELY SOLVING SYSTEM WITHOUT SUBSTITUTION

gamma_optimal = zeros(151, 3); % Container with dimensions 151x3

w_d_hat = ((1 - theta_optimal(:,3))./(1 - theta_optimal(:,1))).^N * w_l *(binom_sum_constructor(B,i,(1/y))/2);
index = w_d_hat < transpose(w_d);

% Can not include indifference relation when gamma_h=0!
for j=1:151;
    if index(j) == 1
        gamma_optimal(j, [1,3]) = lsqnonlin( @(gamma) [ w_d(j)*(1 - theta_optimal(j, 1))^N*binom_sum_constructor(B, i, gamma(1))...
                                                            - w_l*(1 - theta_optimal(j, 3))^N*binom_sum_constructor(B, i, gamma(2)); ...
                                                        gamma(1) + y*gamma(2) - 1]...
                                    , [0.25 0.25] , [0 0], [1 1] );
        gamma_optimal(j, 2) = 0; 
    else     
        gamma_optimal(j, :) = [0, 0, 1/y];
    end;
end;



% gamma_optimal = zeros( 151, 3); % Container with dimensions 151x3
% 
% for j=1:151;
%     gamma_optimal(j, :) = lsqnonlin( @(gamma) [w_d(j)*(1 - theta_optimal(j, 1))^N*binom_sum_constructor(B, i, gamma(1))...
%     - w_l*(1 - theta_optimal(j, 3))^N*binom_sum_constructor(B, i, gamma(3)); ...
%     w_d(j)*(1 - theta_optimal(j, 1))^N*binom_sum_constructor(B, i, gamma(1))...
%     - w_h*(1 - theta_optimal(j, 2))^N*binom_sum_constructor(B, i, gamma(2)); ...
%     gamma(1) + (x-1)*gamma(2) + y*gamma(3) - 1]...
%     , [0.25 0.25 0.25] , [0 0 0], [1 1 1] );
% end;

plot(w_d, gamma_optimal(:,1))
hold on 
plot(w_d, gamma_optimal(:,2))
plot(w_d, gamma_optimal(:,3))

% ## Estimate the polynomial approximation for gamma_d as a function of w_d
% - Only use the approximation where gamma_d is positive
p_gamma_d = polyfit(transpose(w_d(1:30)), gamma_optimal(1:30,1), 1) % Polynomial approximation
% p_gamma_d = polyfit(transpose(w_d(1:30)), gamma_optimal(1:30,1), 5) % Higher order polynomial approximation
p_fit_gamma_d = polyval(p_gamma_d, w_d(1:30)) % Estimate values of theta_d from a polynomial approximation 

% compare the data and approximation
plot(w_d(1:30), p_fit_gamma_d)
hold on
plot(w_d(1:30), gamma_optimal(1:30,1))
hold off

diff_approx_gamma_d = gamma_optimal(1:30,1)-transpose(p_fit_gamma_d) % Difference in errors between simulated values and the approximation
plot(w_d(1:30), diff_approx_gamma_d)

% % In the vicinity of the hat{w}
% 
% p_gamma_d = polyfit(transpose(w_d(25:31)), gamma_optimal(25:31,1), 4) % Polynomial approximation
% p_fit_gamma_d = polyval(p_gamma_d, w_d(25:31)) % Estimate values of theta_d from a polynomial approximation 
% 
% % compare the data and approximation
% plot(w_d(25:31), p_fit_gamma_d)
% hold on
% plot(w_d(25:31), gamma_optimal(25:31,1))
% hold off
% 
% diff_approx_gamma_d = gamma_optimal(25:31,1)-transpose(p_fit_gamma_d) % Difference in errors between simulated values and the approximation
% plot(w_d(25:31), diff_approx_gamma_d)


% % Computing dgamma__dwd
% 
% d_gammad__dw_d = @(w_d) (p_gamma_d(1)*5)*w_d.^4 + (p_gamma_d(2)*4)*w_d.^3 + ...
%     (p_gamma_d(3)*3)*w_d.^2 + (p_gamma_d(4)*2).*w_d + p_gamma_d(5) % dgamma_wd over the entire interval

dgammad_dwd = zeros(1,151);
% dgammad_dwd(1, 1:30) = d_gamma__dw_d(w_d(1:30))
dgammad_dwd(1, 1:30) = p_gamma_d(1);

plot(w_d, dgammad_dwd)
hold on
plot(w_d, dthetad_dwd)

% ### Profits

% Original equilibrium

pi_h = (1 - w_h)*(1 - (1 - 1/x)^N );

% After deviation
pi_d = (1 - w_d).*(1 - (1 - transpose(theta_optimal(:,1))).^N .* (1 - transpose(gamma_optimal(:, 1))).^B) ;

plot(w_d, pi_d) % Looks like pi_d declines smoothly from the separating equilibrium
hold on
plot(w_d, pi_h)


% ### Evaluate the derivative of the profits below w_h: This
%     derivative is positive (I HOPE)

d_pi__d_wd = - (1 - (1 - transpose(theta_optimal(:,1))).^N .* (1 - transpose(gamma_optimal(:,1))).^B )...
    + N*(1 - w_d).*(1 - transpose(theta_optimal(:,1))).^(N-1) .* (1 - transpose(gamma_optimal(:, 1))).^B .* dthetad_dwd...
    + B*(1 - w_d).*(1 - transpose(theta_optimal(:,1))).^N .* (1 - transpose(gamma_optimal(:,1))).^(B-1) .* dgammad_dwd ;

plot(w_d, d_pi__d_wd)
plot(transpose(w_d), theta_optimal(:,1))
hold on
plot(transpose(w_d), gamma_optimal(:,1))

% ### Evaluate the derivative of the profits just at w_h: This derivative
%     should be NEGATIVE (it is, -0.3044)

w_d = w_h
theta_d_sep = 1/x % Since no whites apply when w_d = w_l
gamma_d_sep = 0 % Since no blacks apply when w_d = w_h

dtheta_d__dw_d = binom_sum_constructor(N,i,theta_d_sep)/( w_d*deriv_binom_sum_constructor(N, i, theta_d_sep) + (w_h/(x-1))*deriv_binom_sum_constructor(N, i, (1 - theta_d_sep)/(x-1) ) )

dgamma_d__dw_d = 0;

% Profits 

d_pi_d_wd_initial_high = - (1 - (1 - theta_d_optimal)^N * (1 - gamma_d_sep)^B ) + ...
    (1 - w_d)*N*(1 - theta_d_sep)^(N-1)*(1 - gamma_d_sep)^B * dtheta_d__dw_d + ...
    (1 - w_d)*(1 - theta_d_sep)^N *B* (1 - gamma_d_sep)^(B-1) * dgamma_d__dw_d ;


































