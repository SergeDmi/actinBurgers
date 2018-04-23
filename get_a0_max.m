function [ amin ] = get_a0_max( alph0,DELTA_T,Tmax )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[~,k]=max(alph0);
%k=k-1;
if 10>k
    k=numel(alph0);
end
if k>ceil(1+(Tmax/DELTA_T))
    k=ceil(1+(Tmax/DELTA_T));
end
tt=DELTA_T*(0:(k-1))/100;
aa=alph0(1:k)*100;
%a1=min(aa);
a1=aa(10);
amin=fminsearch(@(a) error_a0(tt,aa,a),a1)/100;

end

function [er]=error_a0(T,A,a0)
er=A-a0*1./(1-a0*T);
%er=sum(abs(er));
er=sum(er.^2);
end
