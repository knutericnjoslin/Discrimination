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

w_d_above = transpose((w_l+0.001):0.001:(w_h-0.01));
len_w_d_above = length(w_d_above);
w_d_below = (w_l-0.05):0.001:(w_l);
len_w_d_below = length(w_d_below);

% Consider a deviation by one of the firms setting w_l upward to some 
% wage w_d such that w_l < w_d < w_h. 


% ### White application behavior - Never apply to firms offering w_l

% theta_h as a function of theta_d
theta_h = @(theta_d, x) (1 - theta_d)/x ;

% ## Solving for theta_d given market constellation using equation 3.3

theta_d_optimal = zeros(1,len_w_d_above);
for j = 1:len_w_d_above;
    theta_d_optimal(j) = lsqnonlin( @(theta_d) w_h*binom_sum_constructor(N, i, theta_h(theta_d, x)) - ...
        w_d_above (j) * binom_sum_constructor(N, i, theta_d) , 0.25, 0, 1 );
end;

% # Plot of the application probabilities of whites

plot(w_d_above, theta_d_optimal)
hold on


% appoximating the derivative and then comparing it with the calculated
% value

p = polyfit(w_d_above, theta_d_optimal, 4)
p_fit = polyval(p, w_d_above)
dp_dw = @(w_d) 4*p(1)*w_d.^3 + 3*p(2)*w_d.^2 + 2*p(3)*w_d + p(4)
test = dp_dw(w_d_above);

plot(w_d_above, p_fit)
hold on
plot(w_d_above, theta_d_optimal)

% dtheta_d__dw_d computed exactly
dtheta_d__dw_d = binom_sum_constructor( N, i, theta_d_optimal ) ./ ( w_d_above.*deriv_binom_sum_constructor( N, i, theta_d_optimal ) + (w_h/x)*deriv_binom_sum_constructor( N, i, theta_h(theta_d_optimal, x)) );   

plot(w_d_above, test)
hold on
plot(w_d_above, dtheta_d__dw_d) % These derivatives match very closely. 

compare = dtheta_d__dw_d - test;
plot(w_d_above, compare)















% ### Black application behavior - 
gamma_l = @(gamma_d, y) (1 - gamma_d)/(y-1);
% 
% gamma_d_optimal = zeros(1,133);
% for j = 1:133;
%    gamma_d_optimal(j) = lsqnonlin( @(gamma_d) w_d(j)*( 1 - theta_d_optimal(j) )^N * binom_sum_constructor(B, i, gamma_d)...
%     - (w_l/(y-1))*binom_sum_constructor(B, i, gamma_l(gamma_d,y)), 0.5,0,1) ;
% end;

% gamma_d_optimal = zeros(1, len_w_d_below + len_w_d_above);

for j = 1:len_w_d_below;
   gamma_d_optimal_below(j) = lsqnonlin( @(gamma_d) w_d_below(j) * binom_sum_constructor(B, i, gamma_d)...
    - (w_l)*binom_sum_constructor(B, i, gamma_l(gamma_d,y)), 0.5,0,1) ;
end;


for k = 1:len_w_d_above;
   gamma_d_optimal_above(k) = lsqnonlin( @(gamma_d) w_d_above(k)*( 1 - theta_d_optimal(k) )^N * binom_sum_constructor(B, i, gamma_d)...
    - (w_l)*binom_sum_constructor(B, i, gamma_l(gamma_d,y)), 0.5,0,1) ;
end;

plot([ w_d_below, w_d_above ], [gamma_d_optimal_below, gamma_d_optimal_above])

% plot( w_d_below, gamma_d_optimal)
% plot( w_d_above, gamma_d_optimal )


% Derivative 
dgamma_d__dw_d = ( (1 - theta_d_optimal).^N .* binom_sum_constructor(B, i, gamma_d_optimal_above) - w_d_above.*N.*(1 - theta_d_optimal).^(N-1).*binom_sum_constructor(B, i, gamma_d_optimal_above).*dtheta_d__dw_d  )...
    ./ ( w_d_above.*( 1 - theta_d_optimal ).^N .* deriv_binom_sum_constructor(B, i, gamma_d_optimal_above) + (w_l/(y-1)).*deriv_binom_sum_constructor(B, i, gamma_l(gamma_d_optimal_above, y) ) ) ;

% approximated derivative 
p_gamma = polyfit(w_d_above, gamma_d_optimal_above, 4);
p_gamma_fit = polyval(p_gamma, w_d_above);
dgamma_dwd = @(w_d) 4*p_gamma(1)*w_d.^3 + 3*p_gamma(2)*w_d.^2 + 2*p_gamma(3)*w_d + p_gamma(4);
test_gamma = dgamma_dwd(w_d_above);

plot(w_d_above, test_gamma)
hold on
plot(w_d_above, dgamma_d__dw_d) % These derivatives are interesting, and curved, and match, but deviate a bit near the boundaries. 





% Profits
% -------

% Original equilibrium

pi_l = (1 - w_l)*(1 - (1 - 1/y)^B );

% After deviation - POSSIBLE TO SHOW THAT THE PROFITS MUST FALL BELOW THE
% EQUILIBRIUM BY LOGIC?

pi_d_below = (1 - w_d_below).*(1 - (1 - gamma_d_optimal_below).^B );
pi_d_above = (1 - w_d_above).*(1 - (1 - theta_d_optimal).^N .* (1 - gamma_d_optimal_above).^B );

% Plot of profits
plot([ w_d_below, w_d_above ], [pi_d_below, pi_d_above])
hold on
plot([ w_d_below, w_d_above ], pi_l)


% ### derivative of profits

% dtheta_d__dw_d
dtheta_d__dw_d = binom_sum_constructor( N, i, theta_d_optimal ) ./ ( w_d_above.*deriv_binom_sum_constructor( N, i, theta_d_optimal ) + (w_h/x)*deriv_binom_sum_constructor( N, i, theta_h(theta_d_optimal, x)) );   

% dgamma_d__dw_d
dgamma_d__dw_d = ( (1 - theta_d_optimal).^N .* binom_sum_constructor(B, i, gamma_d_optimal_above) - w_d_above.*N.*(1 - theta_d_optimal).^(N-1).*binom_sum_constructor(B, i, gamma_d_optimal_above).*dtheta_d__dw_d  )...
    ./ ( w_d_above.*( 1 - theta_d_optimal ).^N .* deriv_binom_sum_constructor(B, i, gamma_d_optimal_above) + (w_l/(y-1)).*deriv_binom_sum_constructor(B, i, gamma_l(gamma_d_optimal_above, y) ) ) ;
% Notice that beyond w(60) = 0.282 no blacks apply, so changes in this
% region do not affect the application probabilities


d_pi__d_wd = - (1 - (1 - theta_d_optimal).^N .* (1 - gamma_d_optimal).^B ) + ...
    (1 - w_d_above).*N.*(1 - theta_d_optimal).^(N-1).*(1 - gamma_d_optimal).^B .* dtheta_d__dw_d + ...
    (1 - w_d_above).*(1 - theta_d_optimal).^N .*B.* (1 - gamma_d_optimal).^(B-1) .* dgamma_d__dw_d ;

plot(w_d_above, d_pi__d_wd)





























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
w_d_below = (w_l+0.001):0.001:(w_h-0.001);
len_w_d_below = length(w_d_below);
w_d_above = (w_h):0.001:(w_h+0.05);
len_w_d_above = length(w_d_above);

% Consider a deviation by one of the firms setting w_h downward to some 
% wage w_d such that w_l < w_d < w_h. Hence x-1 firms setting w_h, y firms 
% setting w_l, and one firm setting w_d.




% ### White application behavior - Whites apply to all firms

% ## Stage 1. Solve for theta_h as a function of theta_d

% Eliminate theta_l by writing as a function of theta_d and theta_h
theta_l = @(theta_d, theta_h, x, y ) ( 1 - theta_d - (x-1).*theta_h )/y ;

theta_optimal = zeros(len_w_d_below + len_w_d_above, 3);

for l=1:len_w_d_below;
    theta_optimal(l,1:2) = lsqnonlin( @(theta) [  w_d_below(l)*binom_sum_constructor(N, i, theta(1))...
                                                    - w_l*binom_sum_constructor(N, i, theta_l(theta(1), theta(2), x, y)); ...
                                                w_d_below(l)*binom_sum_constructor(N, i, theta(1))...
                                                    - w_h*binom_sum_constructor(N, i, theta(2))]...
                                                , [0.5 0.5] , [0 0], [1 1] );
    theta_optimal(l, 3) = theta_l( theta_optimal(l,1), theta_optimal(l,2), x, y);
end;

for l=1:len_w_d_above;
    theta_optimal(l+len_w_d_below,1:2) = lsqnonlin( @(theta) [  w_d_above(l)*binom_sum_constructor(N, i, theta(1))...
                                                    - w_h*binom_sum_constructor(N, i, theta(2));...
                                                theta(1) + (x-1)*theta(2) - 1 ]...
                                                , [0.5 0.5] , [0 0], [1 1] );
    theta_optimal(l+len_w_d_below, 3) = theta_l( theta_optimal(l+len_w_d_below, 1), theta_optimal(l+len_w_d_below, 2), x, y);
end;

plot([w_d_below, w_d_above], theta_optimal(:,1))


% Fitting an approximation to the derivative

% appoximating the derivative and then comparing it with the calculated
% value

p_theta_below = polyfit(transpose(w_d_below), theta_optimal(1:len_w_d_below, 1), 4);
p_theta_below_fit = polyval(p_theta_below, w_d_below);
dtheta_below_dwd = @(w_d) 4*p_theta_below(1)*w_d.^3 + 3*p_theta_below(2)*w_d.^2 + 2*p_theta_below(3)*w_d + p_theta_below(4);
test = dtheta_below_dwd(w_d_below);

plot(w_d_below, p_theta_below_fit)
hold on
plot(w_d_below, theta_optimal(1:len_w_d_below,1))



% dtheta_d__dw_d computed exactly
dtheta_d__dw_d = binom_sum_constructor( N, i, theta_optimal(1:len_w_d_below,1) ) ./... 
                    (... 
                        transpose(w_d_below).*deriv_binom_sum_constructor( N, i, theta_optimal(1:len_w_d_below,1) ) -...
                                (... 
                                w_h*deriv_binom_sum_constructor( N, i, theta_optimal(1:len_w_d_below,2)).*...
                                    ( (w_l/y)*deriv_binom_sum_constructor( N, i, theta_l(theta_optimal(1:len_w_d_below,1), theta_optimal(1:len_w_d_below,2), x, y ) ) )...  
                                        )./ ...
                                            (... 
                                                - w_h*deriv_binom_sum_constructor( N, i, theta_optimal(1:len_w_d_below, 2) )...
                                                - (w_l*((x-1)/y)*deriv_binom_sum_constructor(N, i, theta_l(theta_optimal(1:len_w_d_below,1), theta_optimal(1:len_w_d_below,2), x, y ) )...
                                            )...
                                    )...
                            );   


dtheta_h__dw_d = - dtheta_d__dw_d .*(...
                                        (w_l/y)*deriv_binom_sum_constructor( N, i, theta_l(theta_optimal( 1:len_w_d_below,1), theta_optimal(1:len_w_d_below,2), x, y ) )...
                                            ./(...
                                                w_h*deriv_binom_sum_constructor( N, i, theta_optimal(1:len_w_d_below, 2) )...
                                                + (w_l*((x-1)/y)*deriv_binom_sum_constructor(N, i, theta_l(theta_optimal(1:len_w_d_below,1), theta_optimal(1:len_w_d_below,2), x, y ) ) )... 
                                            )...
                                    );

dtheta_l__dw_d = (1/y).*( - dtheta_d__dw_d - (x-1)*dtheta_h__dw_d) ;                                



plot(w_d_below, test)
hold on
plot(transpose(w_d_below), dtheta_d__dw_d) 

compare = transpose(dtheta_d__dw_d) - test;
plot(transpose(w_d_below), compare)


plot(w_d_below, dtheta_d__dw_d)
hold on
plot(w_d_below, dtheta_h__dw_d)
plot(w_d_below, dtheta_l__dw_d )



%test the dtheta_h_dw_d 

% p_thetah_below = polyfit(transpose(w_d_below), theta_optimal(1:len_w_d_below, 2), 4);
% p_thetah_below_fit = polyval(p_thetah_below, w_d_below);
% dthetah_below_dwd = @(w_d) 4*p_thetah_below(1)*w_d.^3 + 3*p_thetah_below(2)*w_d.^2 + 2*p_thetah_below(3)*w_d + p_thetah_below(4);
% test_h = dthetah_below_dwd(w_d_below);
% 
% plot(w_d_below, p_thetah_below_fit)
% hold on
% plot(w_d_below, theta_optimal(1:len_w_d_below,2))
% 
% Fit of derivative appears to be close:
% plot(w_d_below, test_h)
% hold on
% plot(transpose(w_d_below), dtheta_h__dw_d)  

% ALSO THE FIT OF THE DERIVATIVE OF dtheta_l__dw_d appears close 
% p_thetal_below = polyfit(transpose(w_d_below), theta_optimal(1:len_w_d_below, 3), 4);
% p_thetal_below_fit = polyval(p_thetal_below, w_d_below);
% dthetal_below_dwd = @(w_d) 4*p_thetal_below(1)*w_d.^3 + 3*p_thetal_below(2)*w_d.^2 + 2*p_thetal_below(3)*w_d + p_thetal_below(4);
% test_l = dthetal_below_dwd(w_d_below);
% 
% plot(w_d_below, p_thetal_below_fit)
% hold on
% plot(w_d_below, theta_optimal(1:len_w_d_below,3))
% 
% % Fit of derivative appears to be close:
% plot(w_d_below, test_l)
% hold on
% plot(transpose(w_d_below), dtheta_l__dw_d)  


% ### Black application behavior 

% Eliminate gamma_l by writing as a function of gamma_d and gamma_h
% gamma_l = @(gamma_d, gamma_h, x, y ) ( 1 - gamma_d - (x-1)*gamma_h )/y ;

gamma_optimal = zeros(151, 3); % Container with dimensions 151x3

w_d_hat = ((1 - theta_optimal(1:len_w_d_below,3))./(1 - theta_optimal(1:len_w_d_below,1))).^N * w_l *(binom_sum_constructor(B,i,(1/y))/B);
index = w_d_hat < transpose(w_d_below);

plot(w_d_below, w_d_below)
hold on
plot(w_d_below, w_d_hat)

% Can not include indifference relation when gamma_h=0!
for j=1:len_w_d_below;
    if index(j) == 1
        gamma_optimal(j, [1,3]) = lsqnonlin( @(gamma) [ w_d_below(j)*(1 - theta_optimal(j, 1))^N*binom_sum_constructor(B, i, gamma(1))...
                                                            - w_l*(1 - theta_optimal(j, 3))^N*binom_sum_constructor(B, i, gamma(2)); ...
                                                        gamma(1) + y*gamma(2) - 1]...
                                    , [0.25 0.25] , [0 0], [1 1] );
        gamma_optimal(j, 2) = 0; 
    else     
        gamma_optimal(j, :) = [0, 0, 1/y];
    end;
end;

for j=1:len_w_d_below;
    if index(j) == 1
        gamma_optimal(j, 1) = lsqnonlin( @(gamma) w_d_below(j)*(1 - theta_optimal(j, 1))^N*binom_sum_constructor(B, i, gamma)...
                                                            - w_l*(1 - theta_optimal(j, 3))^N*binom_sum_constructor(B, i, (1-gamma)/y)...
                                    , [0.25] , [0], [1] );
        gamma_optimal(j, 2) = 0 ;
        gamma_optimal(j, 3) = (1 - gamma_optimal(j,1))/y ; 
    else     
        gamma_optimal(j, :) = [0, 0, 1/y];
    end;
end;

plot(w_d_below, gamma_optimal(:,1))
hold on 
plot(w_d_below, gamma_optimal(:,2))
plot(w_d_below, gamma_optimal(:,3))




% dgamma_d__dw_d computed exactly: 

dgamma_d__dw_d = ( ... 
                    (1 - theta_optimal(1:len_w_d_below, 1) ).^N .* binom_sum_constructor( B, i, gamma_optimal(1:len_w_d_below, 1))...
                        - transpose(w_d_below).* N.*(1 - theta_optimal(1:len_w_d_below, 1) ).^(N-1) .* binom_sum_constructor( B, i, gamma_optimal(1:len_w_d_below, 1)).*dtheta_d__dw_d...
                            + w_l * N*(1 - theta_optimal(1:len_w_d_below, 3) ).^(N-1) .* binom_sum_constructor( B, i, (1 - gamma_optimal(1:len_w_d_below, 1))./y ).*dtheta_l__dw_d...
                  )...
                 ./ ...
                 (...  
                     transpose(w_d_below).*(1 - theta_optimal(1:len_w_d_below, 1) ).^N .* deriv_binom_sum_constructor( B, i, gamma_optimal(1:len_w_d_below, 1))...
                        + (w_l/y) * (1 - theta_optimal(1:len_w_d_below, 3) ).^N .* deriv_binom_sum_constructor( B, i, (1 - gamma_optimal(1:len_w_d_below, 1))./y )...
                 );

for j=1:len_w_d_below;
    if index(j) == 1
        dgamma_d__dw_d(j) = ( ... 
                    (1 - theta_optimal(j, 1) ).^N .* binom_sum_constructor( B, i, gamma_optimal(j, 1))...
                        - w_d_below(j).* N.*(1 - theta_optimal(j, 1) ).^(N-1) .* binom_sum_constructor( B, i, gamma_optimal(j, 1)).*dtheta_d__dw_d(j)...
                            + w_l * N*(1 - theta_optimal(j, 3) ).^(N-1) .* binom_sum_constructor( B, i, (1 - gamma_optimal(j, 1))./y ).*dtheta_l__dw_d(j)...
                  )...
                 ./ ...
                 (...  
                     w_d_below(j).*(1 - theta_optimal(j, 1) ).^N .* deriv_binom_sum_constructor( B, i, gamma_optimal(j, 1))...
                        + (w_l/y) * (1 - theta_optimal(j, 3) ).^N .* deriv_binom_sum_constructor( B, i, (1 - gamma_optimal(j, 1))./y )...
                 );

    else     
        dgamma_d__dw_d(j) = 0;
    end;
    
end;



% Approximating the derivative             
% p_gamma_below = polyfit(transpose(w_d_below(1:30)), gamma_optimal(1:30, 1), 4);
p_gamma_below = polyfit(transpose(w_d_below(1:30)), gamma_optimal(1:30, 1), 2);
p_gamma_below_fit = polyval(p_gamma_below, w_d_below(1:30));
% dgamma_below_dwd = @(w_d) 4*p_gamma_below(1)*w_d.^3 + 3*p_gamma_below(2)*w_d.^2 + 2*p_gamma_below(3)*w_d + p_gamma_below(4);
dgamma_below_dwd = @(w_d) 2*p_gamma_below(1)*w_d + p_gamma_below(2)
test2 = dgamma_below_dwd(w_d_below(1:30));

lin_slope = (gamma_optimal(2:len_w_d_below, 1)-gamma_optimal(1:len_w_d_below-1, 1))./.001;

% fitted model -- Not so good
plot(transpose(w_d_below(1:30)), p_gamma_below_fit(1:30))
hold on
plot(transpose(w_d_below(1:30)), gamma_optimal(1:30, 1))


plot(w_d_below, dgamma_d__dw_d)
hold on
% plot(w_d_below, test2)
plot(w_d_below , [lin_slope;0])

























% ### Profits

% Original equilibrium

pi_h = (1 - w_h)*(1 - (1 - 1/x)^N );

% After deviation
pi_d_below = (1 - transpose(w_d_below)).*(1 - (1 - theta_optimal(1:len_w_d_below, 1)).^N .* (1 - gamma_optimal(1:len_w_d_below, 1)).^B) ;
pi_d_above = (1 - transpose(w_d_above)).*(1 - (1 - theta_optimal(1+len_w_d_below: len_w_d_below+len_w_d_above, 1)).^N ) ;

plot([ transpose(w_d_below); transpose(w_d_above) ], [pi_d_below; pi_d_above])
hold on
plot([ w_d_below, w_d_above ], pi_h)
















% ### Evaluate the derivative of the profits below w_h: This
%     derivative is positive (I HOPE)

d_pi__d_wd = - (1 - (1 - theta_optimal(1:len_w_d_below,1)).^N .* (1 - gamma_optimal(:,1)).^B )...
    + N*(1 - transpose(w_d_below)).*(1 - theta_optimal(1:len_w_d_below,1)).^(N-1) .* (1 - gamma_optimal(:, 1)).^B .*  dtheta_d__dw_d...
    + B*(1 - transpose(w_d_below)).*(1 - theta_optimal(1:len_w_d_below,1)).^N .* (1 - gamma_optimal(:,1)).^(B-1) .*  transpose(dgamma_d__dw_d) ;

plot(w_d_below, d_pi__d_wd)
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














































