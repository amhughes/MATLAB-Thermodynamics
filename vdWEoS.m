function [ V,it ] = vdWEoS( T,Tc,P,Pc,R,w,type )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
es = 0.00001;
itstop = 1000;

Tr = T/Tc;
Pr = P/Pc;

sigma = 0;
epsilon = 0;
omega = 1/8;
psi = 27/64;
alpha = @(Tr) 1;
beta = omega*Pr/Tr;
q = psi*alpha(Tr)/omega/Tr;

if type == 'v'
    Z = 1;
    n = 0;
    while n < itstop
        Zold = Z;
        Z = 1 + beta - q*beta*(Z-beta)/((Z+epsilon*beta)*(Z+sigma*beta));
        ea = abs((Z-Zold)/Z);
        n = n + 1;
        if ea <= es
            break
        end
    end
else
    Z = 1;
    n = 0;
    while n < itstop
        Zold = Z;
        Z = beta + (Z+epsilon*beta)*(Z+sigma*beta)*(1+beta-Z)/(q*beta);
        ea = abs((Z-Zold)/Z);
        n = n + 1;
        if ea <= es
            break
        end
    end
end
if n>=itstop
    disp('Convergence Warning')
end
V = Z*R*T/P;
it = n;
end

