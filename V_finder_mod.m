function [ Vs,Vn,dvg ] = V_finder_mod(Rs,As,Es,Ls,dR,N,Vn,a,m,c,AN,A0,Pa,Pm,Pv,Pl,THRESH,NVMAX,EDV)
%Finds the Vn such that the velocity Vs has the correct behaviour
%   More precisely, V is integrated from the stress balance
%   Additionally, V(0)=0 and the normal stress is 0 at R
Vs=zeros(1,N);
dVs=zeros(1,N);
b=1000;
Xs=(1+b*As.^Pv)/(1+b*A0^Pv);
dXs(1,2:N-1)=(Xs(1,3:N)-Xs(1,1:N-2))/(2*dR);

% Special function depending on the stress
Fs=stress_grad(As,Es,Ls,N,dR,a,m,Pa,Pm,Pl,AN);
% Normal stress free
%Vnm1=Vn-grad_V(As(1,N),Es(1,N),a,m,Pa,Pm,c)*dR;
Vnm1=Vn-grad_V(AN,Es(1,N),a,m,Pa,Pm,Pl,c,Xs(N),Vn,Ls(N))*dR;

% Getting the speed
Vs=integrate_V(Vs,Rs,Fs,dR,N,Vn,Vnm1,Xs,dXs);

% Converging towards good Vn
error=errV(Vs,dR);
dvg=0;

ntry=0;


while abs(error)>THRESH
    % Counding
    ntry=ntry+1;
    % Variational
    dVn=Vn+EDV;
    dVnm1=dVn-grad_V(As(1,N),Es(1,N),a,m,Pa,Pm,Pl,c,Xs(N),dVn,Ls(N))*dR;
    dVs(1,:)=integrate_V(Vs,Rs,Fs,dR,N,dVn,dVnm1,Xs,dXs);
    dVc1=errV(dVs,dR);
    %error
    dVdE=(dVc1-error)/EDV;
    %Following gradient
	Vn=Vn-(0.999+0.0001*rand)*error/dVdE;
    %Vnm1=Vn-grad_V(As(1,N),Es(1,N),a,m,Vn,Rs(1,N),A0)*dR;
    Vnm1=Vn-grad_V(As(1,N),Es(1,N),a,m,Pa,Pm,Pl,c,Xs(N),Vn,Ls(N))*dR;
    Vs=integrate_V(Vs,Rs,Fs,dR,N,Vn,Vnm1,Xs,dXs);
    error=errV(Vs,dR);
    if ntry>NVMAX
        dvg=1;
        error=0;
        disp('Number of tries exceeded')
    end
end


if isnan(sum(Vs))
    dvg=1
end



end



function [dV]=grad_V(An,En,a,m,Pa,Pm,Pl,c,Xn,Vn,Ln)
    
    %% working one
    dV=(0-c*m*An^Pm-a*(An^Pa)*(1-En))/Xn;
    
end

function [err]=errV(V,dR)
    err=1.5*V(1)-V(2)/2.0;
    %err=V(1);
end

