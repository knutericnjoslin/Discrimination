function [ theta_l ] = app_prob( x , y, theta_h )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

theta_l = (1 - ( x - 1) * theta_h )/y;

end

