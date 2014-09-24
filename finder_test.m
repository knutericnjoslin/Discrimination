% LOCAL OPTIMALITY HIGH DEVIATION 
% ------------------------------

% *Original equilibrium, with x firm setting w_h and y firms setting w_l
% *N white workers and B black workers

% Function parameters
i=1;

% Market parameters
N=4;
B=2;
x=3;
y=2;

w_h_value=0.10:0.01:0.45;
w_h_val_length = length(w_h_value);

for k=1:w_h_val_length;
    
w_h = w_h_value(k);     
w_l = w_h*(1 - ( 1 - 1/x)^N )/(N/x) ;
w_d_below = (w_h-0.01):0.001:(w_h+0.01);
len_w_d_below = length(w_d_below);

bf = 0;    
apply_all = 1;
index = [];

% ### White application behavior - Whites apply to all firms

% Eliminate theta_l by writing as a function of theta_d and theta_h
theta_l = @(theta_d, theta_h, x, y ) ( 1 - theta_d - (x-1).*theta_h )/y ;

theta_optimal = zeros(len_w_d_below, 3);

for l=1:len_w_d_below;
     if apply_all == 1;
        theta_optimal(l,1:2) = lsqnonlin( @(theta) [  w_d_below(l)*binom_sum_constructor(N, i, theta(1))...
                                                        - w_l*binom_sum_constructor(N, i, theta_l(theta(1), theta(2), x, y)); ...
                                                    w_d_below(l)*binom_sum_constructor(N, i, theta(1))...
                                                        - w_h*binom_sum_constructor(N, i, theta(2))]...
                                                    , [0.5 0.5], [0 0], [1 1] );
        theta_optimal(l, 3) = theta_l( theta_optimal(l,1), theta_optimal(l,2), x, y);
        if theta_optimal(l,1) + theta_optimal(l,2) > .9995, apply_all = 0; end; 
    else
        theta_optimal(l, 1:2) = lsqnonlin( @(theta) [  w_d_below(l)*binom_sum_constructor(N, i, theta(1))...
                                                    - w_h*binom_sum_constructor(N, i, theta(2));...
                                                theta(1) + (x-1)*theta(2) - 1 ]...
                                                , [0.5 0.5], [0 0], [1 1] );
        theta_optimal(l, 3) = 0;
    end;
end;


% dtheta_d__dw_d computed exactly

dtheta_d__dw_d = zeros(len_w_d_below, 1);

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

dtheta_h__dw_d = zeros(len_w_d_below);
                        
dtheta_h__dw_d = - dtheta_d__dw_d .*(...
                                        (w_l/y)*deriv_binom_sum_constructor( N, i, theta_l(theta_optimal( 1:len_w_d_below,1), theta_optimal(1:len_w_d_below,2), x, y ) )...
                                            ./(...
                                                w_h*deriv_binom_sum_constructor( N, i, theta_optimal(1:len_w_d_below, 2) )...
                                                + (w_l*((x-1)/y)*deriv_binom_sum_constructor(N, i, theta_l(theta_optimal(1:len_w_d_below,1), theta_optimal(1:len_w_d_below,2), x, y ) ) )... 
                                            )...
                                    );

dtheta_l__dw_d = (1/y).*( - dtheta_d__dw_d - (x-1)*dtheta_h__dw_d) ;                                


                        


w_d_hat = ((1 - theta_optimal(1:len_w_d_below,3))./(1 - theta_optimal(1:len_w_d_below,1))).^N * w_l *(binom_sum_constructor(B,i,(1/y))/B);
index = w_d_hat < transpose(w_d_below);



% Computing gamma and dgamma_dwd

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


% ### Profits

% Original equilibrium

d_pi__d_wd = - (1 - (1 - theta_optimal(:,1)).^N .* (1 - gamma_optimal(:,1)).^B )...
    + N*(1 - transpose(w_d_below)).*(1 - theta_optimal(:,1)).^(N-1) .* (1 - gamma_optimal(:, 1)).^B .*  dtheta_d__dw_d...
    + B*(1 - transpose(w_d_below)).*(1 - theta_optimal(:,1)).^N .* (1 - gamma_optimal(:,1)).^(B-1) .*  transpose(dgamma_d__dw_d) ;

    for t=2:len_w_d_below;
        hi_cand = w_d_below(t);
        low_cand = hi_cand*(1 - ( 1 - 1/x)^N )/(N/x);
        candidate = [hi_cand, low_cand];
        if d_pi__d_wd(t) < 0 && d_pi__d_wd(t-1)>0; %  && (w_d_below(t)==w_h(k) | w_d_below(t+1)==w_h(k) | w_d_below(t-1)==w_h(k))
            bf = 1;
            break;
        end;
    end;

    if bf==1;
        disp(candidate); 
        pi_h = (1 - w_h)*(1 - (1 - 1/x)^N );
        pi_l = (1 - w_l)*(1 - (1 - 1/y)^B );
        disp([pi_h, pi_l]);
        break;
    end;

end;











































