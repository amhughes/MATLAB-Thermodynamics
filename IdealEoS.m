function [ val ] = IdealEoS(P,V,n,T,U)
%IDEALEOS This function computes the missing value of the Ideal Gas Law
%   To use replace the missing value with 0
%   The value of R is set by the U parameter
%   Defaults to 0.08206 L*atm/mol*K
%   For molar volume, replace n with 1
%   U =
%       0 - 0.08206 L*atm/mol*K
%       1 - 0.08314 L*bar/mol*K
%       2 - 10.73 ft^3*psi/lb-mol*R
%       3 - 0.7302 ft^3*atm/lb-mol*R
if nargin == 4
    U = 0;
end

switch U
    case 0
        R = 0.08206;
    case 1
        R = 0.08314;
    case 2
        R = 10.73;
    case 3
        R = 0.7302;
    otherwise
        R = 0.08206;
end

if P == 0
    val = n*R*T/V;
elseif V == 0
    val = n*R*T/P;
elseif n == 0
    val = P*V/R/T;
else
    val = P*V/n/R;
end
    
end

