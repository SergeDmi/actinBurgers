function [ slope,offset,score,t,vars] = ortho_robust_coeff( PTS )
% Gives the slope fitting the cloud of points PTS
% Pts is as :
%   x1 x2 x3 ... xn
%   y1 y2 y3 ... yn
% Uses robust orthogonal fitting
% S. Dmitrieff 2015 

t0=pi/2;
t=fminsearch(@(t) out_of_lineness2D(PTS,t),t0);
[score,~,res,unres]=out_of_lineness2D(PTS,t);
slope=-tan(t);
x=mean(PTS(1,:));
y=mean(PTS(2,:));
offset=y-x*slope;
vars=var(res)/(size(PTS,2)*var(unres)); 
end

