function [ ool,pts,res,unres ] = out_of_lineness2D(PTS,theta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

O=mean(PTS,2);
s=size(PTS);
pts=PTS-O*ones(1,s(2));
R=[cos(theta) -sin(theta); sin(theta) cos(theta)];
pts=R*pts;
res=pts(2,:);
unres=pts(1,:);
ool=sum(abs(res))/s(2);
end

