% FINDER TEST
% ===========

% ### FUNCTION HANDLES

% White worker function handles

theta_l = @(theta_d, theta_h, x, y ) ( 1 - theta_d - (x-1).*theta_h )/y ;

% # Indifference relations

% Writing the indifference equations as function handles

white_indiff_wd_wl = @(w_d, w_l, x, y, N, i, theta_d, theta_h) w_d*binom_sum_constructor(N, i, theta_d)...
    - w_l*binom_sum_constructor(N, i, theta_l(theta_d, theta_h, x, y));

white_indiff_wd_wh = @(w_d, w_h, x, y, N, i, theta_d, theta_h) w_d*binom_sum_constructor(N, i, theta_d)...
    - w_h*binom_sum_constructor(N, i, theta_h);

% # Derivatives

% dtheta_d__dw_d for wages below and above w_h
fxn_dtheta_d__dw_d = @(w_d, w_h, w_l, x, y, N, i, theta_d, theta_h) binom_sum_constructor( N, i, theta_d ) ./... 
                    (... 
                        w_d.*deriv_binom_sum_constructor( N, i, theta_d ) - ...
                                (... 
                                w_h*deriv_binom_sum_constructor( N, i, theta_h).*...
                                    ( (w_l/y)*deriv_binom_sum_constructor( N, i, theta_l(theta_d, theta_h, x, y ) ) )...  
                                        )./ ...
                                            (... 
                                                - w_h*deriv_binom_sum_constructor( N, i, theta_h )...
                                                - (w_l*((x-1)/y)*deriv_binom_sum_constructor(N, i, theta_l(theta_d, theta_h, x, y ) )...
                                            )...
                                    )...
                            );   

fxn_dtheta_d__dw_d_above = @( w_d, w_h, x, N, i, theta_d) binom_sum_constructor(N,i,theta_d)/...
                            (...
                                w_d*deriv_binom_sum_constructor(N,i,theta_d)...
                                + (w_h/(x-1))*deriv_binom_sum_constructor(N,i,(1-theta_d)/(x-1))...
                            );
                        
% dtheta_h__dw_d for wages below and above w_h                   
fxn_dtheta_h__dw_d = @(w_d, w_h, w_l, x, y, N, i, theta_d, theta_h, dtheta_d__dw_d) - dtheta_d__dw_d .*(...
                                        (w_l/y)*deriv_binom_sum_constructor( N, i, theta_l(theta_d, theta_h, x, y ) )...
                                            ./(...
                                                w_h*deriv_binom_sum_constructor( N, i, theta_h )...
                                                + (w_l*((x-1)/y)*deriv_binom_sum_constructor(N, i, theta_l(theta_d, theta_h, x, y ) ) )... 
                                            )...
                                    );

fxn_dtheta_h__dw_d_above = @(x, dtheta_d__dw_d) -(1/(x-1))*dtheta_d__dw_d; 

% dtheta_l__dw_d, only relevant for wages below w_h
fxn_dtheta_l__dw_d = @( x, y, dtheta_d__dw_d, dtheta_h__dw_d) (1/y).*( - dtheta_d__dw_d - (x-1)*dtheta_h__dw_d) ;  

% ### Black function handles

% Black indifference condition
black_indiff_wd_wl = @(w_d, w_l, x, y, N, i, theta_d, theta_l, gamma_d) w_d * (1 - theta_d )^N * binom_sum_constructor(B, i, gamma_d)...
    - (w_l)* (1 - theta_l)^N *binom_sum_constructor(B, i, (1-gamma_d)/y) ;

% Function handle for dgamma_d__dw_d
fxn_dgamma_d__dw_d = @(w_d, w_l, x, y, N, B, i, theta_d, dtheta_d__dw_d, gamma_d)( (1 - theta_d).^N .* binom_sum_constructor(B, i, gamma_d)... 
                    - w_d.*N.*(1 - theta_d).^(N-1).*binom_sum_constructor(B, i, gamma_d).*dtheta_d__dw_d  )...
                    ./... 
                    ( w_d.*( 1 - theta_d ).^N .* deriv_binom_sum_constructor(B, i, gamma_d)... 
                    + (w_l/y).*deriv_binom_sum_constructor(B, i, (1 - gamma_d)/y ) ) ;

% MAIN SCRIPT
% -----------

i=1;

% Market parameters
N=3;
B=2;
x=3;
y=2;

w_h_value=0.20:0.01:0.55;
w_h_val_length = length(w_h_value);

for k=1:w_h_val_length;
    
w_h = w_h_value(k);     
w_l = w_h*(1 - ( 1 - 1/x)^N )/(N/x) ;
w_d = transpose((w_h-0.15):0.001:(w_h+0.15));
len_w_d = length(w_d);

bf = 0;    
apply_all = 1;
index = [];


% APPLICATION PROBABILITIES

% Computing the optimal application probabilities for whites

theta_optimal = zeros(len_w_d, 3);
dtheta_dwd = zeros(len_w_d, 3);


for l=1:len_w_d;
    
    if ( theta_optimal(l,3) < 0.0001);
        theta_optimal(l, 1:2) = lsqnonlin( @(theta) [  white_indiff_wd_wh( w_d(l), w_h, x, y, N, i, theta(1), theta(2));...
                                                theta(1) + (x-1)*theta(2) - 1 ]...
                                                , [0.5 0.5], [0 0], [1 1] );
        theta_optimal(l, 3) = 0;
    
        dtheta_dwd( l, 1) = fxn_dtheta_d__dw_d_above( w_d(l), w_h, x, N, i, theta_optimal(l,1) );
        dtheta_dwd( l, 2) = fxn_dtheta_h__dw_d_above( x, dtheta_dwd( l, 1) );
        dtheta_dwd( l, 3) = 0;
    else 
        theta_optimal(l,1:2) = lsqnonlin( @(theta) [ white_indiff_wd_wl( w_d(l), w_l, x, y, N, i, theta(1), theta(2) ); ...
                                                     white_indiff_wd_wh( w_d(l), w_h, x, y, N, i, theta(1), theta(2) ) ]...
                                                    , [0.5 0.5], [0 0], [1 1] );
        theta_optimal(l, 3) = theta_l( theta_optimal(l,1), theta_optimal(l,2), x, y);
    
        dtheta_dwd( l, 1) = fxn_dtheta_d__dw_d(w_d(l), w_h, w_l, x, y, N, i, theta_optimal(l,1), theta_optimal(l,2));
        dtheta_dwd( l, 2) = fxn_dtheta_h__dw_d(w_d(l), w_h, w_l, x, y, N, i, theta_optimal(l,1), theta_optimal(l,2), dtheta_dwd(l,1));
        dtheta_dwd( l, 3) = fxn_dtheta_l__dw_d( x, y, dtheta_dwd(l,1), dtheta_dwd(l,2));        
    end;
    
end;

for l=1:len_w_d;
     
    if apply_all == 1;
        theta_optimal(l,1:2) = lsqnonlin( @(theta) [ white_indiff_wd_wl( w_d(l), w_l, x, y, N, i, theta(1), theta(2) ); ...
                                                     white_indiff_wd_wh( w_d(l), w_h, x, y, N, i, theta(1), theta(2) ) ]...
                                                    , [0.5 0.5], [0 0], [1 1] );
        theta_optimal(l, 3) = theta_l( theta_optimal(l,1), theta_optimal(l,2), x, y);
    
        dtheta_dwd( l, 1) = fxn_dtheta_d__dw_d(w_d(l), w_h, w_l, x, y, N, i, theta_optimal(l,1), theta_optimal(l,2));
        dtheta_dwd( l, 2) = fxn_dtheta_h__dw_d(w_d(l), w_h, w_l, x, y, N, i, theta_optimal(l,1), theta_optimal(l,2), dtheta_dwd(l,1));
        dtheta_dwd( l, 3) = fxn_dtheta_l__dw_d( x, y, dtheta_dwd(l,1), dtheta_dwd(l,2));
        
        if (theta_optimal(l,3) < 0.0005), apply_all = 0; end; 
    
    else
        theta_optimal(l, 1:2) = lsqnonlin( @(theta) [  white_indiff_wd_wh( w_d(l), w_h, x, y, N, i, theta(1), theta(2));...
                                                theta(1) + (x-1)*theta(2) - 1 ]...
                                                , [0.5 0.5], [0 0], [1 1] );
        theta_optimal(l, 3) = 0;
    
        dtheta_dwd( l, 1) = fxn_dtheta_d__dw_d_above( w_d(l), w_h, x, N, i, theta_optimal(l,1) );
        dtheta_dwd( l, 2) = fxn_dtheta_h__dw_d_above( x, dtheta_dwd( l, 1) );
        dtheta_dwd( l, 3) = 0;
        
    end;
    
end;


% Application index
w_d_hat = ((1 - theta_optimal(:,3))./(1 - theta_optimal(:,1))).^N * w_l *(binom_sum_constructor(B,i,(1/y))/B);
index = w_d_hat < w_d;


% ### BLACK BEHAVIOR

% initialize vectors
gamma_optimal = zeros(len_w_d, 3);
dgamma_dwd = zeros(len_w_d, 1);

for j=1:len_w_d;
    if index(j) == 1
        gamma_optimal( j, 1) = lsqnonlin( @(gamma) black_indiff_wd_wl( w_d(j), w_l, x, y, N, i, theta_optimal(j,1), theta_optimal(j,3), gamma),... 
            0.25 , 0, 1 );
        gamma_optimal( j, 2) = 0 ;
        gamma_optimal( j, 3) = (1 - gamma_optimal(j,1))/y ; 
        dgamma_dwd( j, 1) = fxn_dgamma_d__dw_d( w_d(j), w_l, x, y, N, B, i, theta_optimal(j,1), dtheta_dwd(j,1), gamma_optimal(j,1));
    else     
        gamma_optimal( j, :) = [0, 0, 1/y];
        dgamma_dwd( j, 1) = 0;
    end;
end;


% PROFITS

% # Profits
pi_h = (1 - w_h)*(1 - (1 - 1/x)^N ); % original equilibrium
pi_d = (1 - w_d).*(1 - (1 - theta_optimal(:, 1)).^N .* (1 - gamma_optimal(:, 1)).^B) ;

% # Derivative of profits
d_pi__d_wd = - (1 - (1 - theta_optimal(:,1)).^N .* (1 - gamma_optimal(:,1)).^B ) + ...
    (1 - w_d).*N.*(1 - theta_optimal(:,1)).^(N-1).*(1 - gamma_optimal(:,1)).^B .* dtheta_dwd(:,1) + ...
    (1 - w_d).*(1 - theta_optimal(:,1)).^N .*B.* (1 - gamma_optimal(:,1)).^(B-1) .* dgamma_dwd ;

    
    for t=(149):1:(152);
        hi_cand = w_d(t-1);
        low_cand = hi_cand*(1 - ( 1 - 1/x)^N )/(N/x);
        candidate = [hi_cand, low_cand];
        if d_pi__d_wd(t) < 0 && d_pi__d_wd(t-2)>0 ; %  && (w_d_below(t)==w_h(k) | w_d_below(t+1)==w_h(k) | w_d_below(t-1)==w_h(k))
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





plot(w_d, theta_optimal(:,1))
hold on
plot(w_d, theta_optimal(:,2))
plot(w_d, theta_optimal(:,3))

plot(w_d, dtheta_dwd(:,1))
hold on
plot(w_d, dtheta_dwd(:,2))
plot(w_d, dtheta_dwd(:,3))

plot( w_d, gamma_optimal(:,1))
hold on
plot( w_d, gamma_optimal(:,3))

plot(w_d, pi_d)
hold on
plot(w_d, pi_h)

plot(w_d, d_pi__d_wd)


















