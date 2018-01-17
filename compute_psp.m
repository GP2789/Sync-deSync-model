function [ EPSP ] = compute_psp( tau_syn )
%COMPUTE_PSP Summary of this function goes here
%   Detailed explanation goes here
T = 1:1:10000;
EPSP = (exp(1)*T./tau_syn).*exp(-T./tau_syn).*heaviside(T); 
EPSP = transpose(EPSP(EPSP > 0.01)); 

end

