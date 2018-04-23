function [ Vs ] = integrate_V(Vs,Rs,Fs,dR,N,Vn,Vnm1,Xs,dXs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Vs(1,N)=Vn;
Vs(1,N-1)=Vnm1;

for n=N-1:-1:2
    %Vs(n-1)=(-Vs(n+1)*(1+dR/Rs(n))+2*Vs(n)*(1+(dR/Rs(n))^2)-Fs(n)*dR^2)/(1-dR/Rs(n));
     Vs(n-1)=(-Vs(n+1)*(1+dR/Rs(n)+0.5*dXs(n)*dR/Xs(n))+2*Vs(n)*(1+(dR/Rs(n))^2)-(Fs(n)*dR^2)/Xs(n))/(1-dR/Rs(n)-0.5*dXs(n)*dR/Xs(n));
end

end

