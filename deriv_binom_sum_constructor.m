function [ sum ] = deriv_binom_sum_constructor( T, i, app_prob)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Function to construct binomial sums 

% Inputs
% ------
% ### Market parameters:
%   * The number of white firms, x, and the number of black firms, y.
%   * The wages offered in the market, w_l, w_h, w_d
%   * 
% ### Other paramters:
%   * Starting index for the sum, i.
%   * Number to be subtracted from exponent in binomial expansion


sum = 0;
for k=i:T
    sum = sum + nchoosek(T,k)*(k-1)*(- app_prob ).^(k-2);
end

end
