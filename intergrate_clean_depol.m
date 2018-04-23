function [ sA,sE,sV,sL,sT,dvg ] = intergrate_clean_depol(data,cstes,pars)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%disp(pars)
% Rs,As,Es,dR,Vn,a,m,G,Ge,dt,Tobs,A0,E0
%% Physical parameters
%% Physical parameters
a=pars(1);  % Elasticity
m0=pars(2); % Contraction
c=pars(3); % Adhesion (epsilon)
Pa=pars(4);%  Power law for elasticity
Pm=pars(5);%  Power law for contraction
Pv=pars(6);%  Power law for viscosity
fe=pars(7);%  elastictiy relaxation rate
Pl=pars(8);%   1 if contraction depends on filament end density (model D)
%              0 otherwise (model M)


%% From the exp data
Ks=data.Ks;
G=Ks(1);
m=m0*G;
Ge=G;
ARs=data.ARs;
A0=ARs(1);
AN=ARs(1);
AA=data.datas;
Ai=AA(1,:);
li=length(Ai);
%expo=0.276;
expo=data.expo;

%% Other parameters
Nobs=data.Nobs;
Tobs=data.Tobs;
Tend=Tobs(end);

%% Snuff
Rs=cstes.Rs;
dR=cstes.dR;
Vn=cstes.Vn;
dt=cstes.dt;
Ns=length(Rs);
%As=spline((1:li)/li,Ai,Rs);
%As=ones(1,Ns)*mean(Ai);
%% ZDLIZLDJLZDJ
% qzddzqdq
As=ones(1,Ns)*AN;
%dqzdk,z:qd,
%% jldlqjidzljqzldj
Ls=ones(1,Ns);
dLdR=ones(1,Ns-1);
VDS=ones(1,Ns);
Gs=ones(1,Ns);
Ges=ones(1,Ns);
Es=ones(1,Ns);
E0=cstes.E0;
L0=1.0;


R2=(Rs-dR/2.0).^2;
R2(1)=dR*dR;

MAS=As.*R2;
MES=Es.*R2;
MAVS=As.*R2;
MEVS=Es.*R2;

%disp(Tend)
%disp(dt)
%% Annoying discrete parameters
THRESH=0.01*dR;
NVMAX=1000;
EDV=0.000000001;

dts2dR=0.5*dt/dR;

%% Making sure time step is adequate
if abs(m*dt/dR)>1
    disp('WARNING : most likely unstable !!!!!!!!')
end


%Nobs=size(Tobs,2);
%Tend=Tobs(Nobs);
sA=zeros(Nobs,Ns);
sL=zeros(Nobs,Ns);
sE=zeros(Nobs,Ns);
sV=zeros(Nobs,Ns);
sR=zeros(1,Nobs);
sT=ones(1,Nobs);
%FsA=zeros(1,Ns+1);
%FsE=zeros(1,Ns+1);
[ Vs,Vn,~] = V_finder_mod(Rs,As,Es,Ls,dR,Ns,Vn,a,m,c,AN,A0,Pa,Pm,Pv,Pl,THRESH,NVMAX,EDV);

dvg=0;
t=0;
nobs=1;
nextT=Tobs(nobs);
t=nextT;
Tpos=1;


gg=0.1;

Vmax=max(1.0,abs(m*A0+a*A0*(1-E0)));
dtmax=abs(0.4*dR/Vmax);

vdep=Ks(nobs);
%vdep=Ks(nobs);
vrel=vdep*fe;

%% Pre-run to equilibrate the simulation
while t<60*data.K0 && dvg==0
    t=t+dt;
  if dvg==0
        %% Applying conservation laws
        % Decomposing the equation between semi-implicit and Lax-Friedrich
        VDS(1,:)=vdep*(1+exp(expo))./(1+exp(expo./Ls(1,:)));
        Gs(1,:)=vdep./Ls(1,:); 
        Ges(1,:)=vrel./Ls(1,:); 
        if max(Gs)*dt>0.5
            disp('Depolymerization rate too high, decrease time step');
        end
        
       
        a_2=As(1,2);
        e_2=Es(1,2);
        l_2=Ls(1,2);
         
       
      
       
        MAVS(:)=As.*R2.*Vs;
        MEVS(:)=Es.*R2.*Vs;
       
        
        Es(1,2:Ns-1)=((Es(1,1:Ns-2)+Es(1,3:Ns)+gg*Es(1,2:Ns-1))/(2.0+gg)-dts2dR*(MEVS(1,3:Ns)-MEVS(1,1:Ns-2))./R2(1,2:Ns-1)+Ges(1,2:Ns-1)*dt)./(1+Ges(1,2:Ns-1)*dt);
        Es(1,1)=((2.0*e_2+gg*Es(1,1))/(2.0+gg)-dts2dR*(MEVS(1,2)-MEVS(1,1))./R2(1)+Ges(1,1)*dt)/(1+Ges(1,1)*dt);
        Es(1,Ns)=E0;
        
  end
  
end
dvg=0;
t=0;
nobs=1;
nextT=Tobs(nobs);
t=nextT;
Tpos=1;

%% Main run
while t<Tend && dvg==0
    %Saving
    if t>=nextT
        Tobs(nobs)=t;
        sA(nobs,:)=As(1,:);
        sE(nobs,:)=Es(1,:);
        sV(nobs,:)=Vs(1,:);
        sL(nobs,:)=Ls(1,:);
        [~,n]=max(As);
        sR(nobs)=n*dR;
        nobs=nobs+1;
        nextT=Tobs(nobs);
        sT(nobs)=Tpos;
        vdep=Ks(nobs);
        m=m0*Ks(nobs);
        vrel=vdep*fe;
        AN=ARs(nobs);
    end
    t=t+dt;
    % Computing the velocities
    [ Vs,Vn,dvg ] = V_finder_mod(Rs,As,Es,Ls,dR,Ns,Vn,a,m,c,AN,A0,Pa,Pm,Pv,Pl,THRESH,NVMAX,EDV);
    %Vs=-(1:Ns)/Ns;
    if dvg==0
        %% Applying conservation laws
        % Decomposing the equation between semi-implicit and Lax-Friedrich
        
        VDS(1,:)=vdep*(1+exp(expo))./(1+exp(expo./Ls(1,:)));
        Gs(1,:)=vdep./Ls(1,:); 
        Ges(1,:)=vrel./Ls(1,:); 
        if max(Gs)*dt>0.5
            disp('Depolymerization rate too high, decrease time step');
        end
        
        a_2=As(1,2);
        e_2=Es(1,2);
        l_2=Ls(1,2);
           
          
       
        MAVS(:)=As.*R2.*Vs;
        MEVS(:)=Es.*R2.*Vs;
        
        As(1,2:Ns-1)=((As(1,1:Ns-2)+As(1,3:Ns)+gg*As(1,2:Ns-1))/(2.0+gg)-dts2dR*(MAVS(1,3:Ns)-MAVS(1,1:Ns-2))./R2(1,2:Ns-1))./(1+Gs(1,2:Ns-1)*dt);
        As(1,1)=((2.0*a_2+gg*As(1,1))/(2.0+gg)-dts2dR*(MAVS(1,2)-MAVS(1,1))/R2(1))./(1+Gs(1,1)*dt);
        As(1,Ns)=AN;
        
        Es(1,2:Ns-1)=((Es(1,1:Ns-2)+Es(1,3:Ns)+gg*Es(1,2:Ns-1))/(2.0+gg)-dts2dR*(MEVS(1,3:Ns)-MEVS(1,1:Ns-2))./R2(1,2:Ns-1)+Ges(1,2:Ns-1)*dt)./(1+Ges(1,2:Ns-1)*dt);
        Es(1,1)=((2.0*e_2+gg*Es(1,1))/(2.0+gg)-dts2dR*(MEVS(1,2)-MEVS(1,1))./R2(1)+Ges(1,1)*dt)/(1+Ges(1,1)*dt);
        Es(1,Ns)=E0;
        
        %% Upstream method for the size of filaments, because conservation equation on length is *messy*
        Ls(1,Ns)=L0;
        Ls(1,2:Ns-1)=(Ls(1,1:Ns-2)+Ls(1,3:Ns)+gg*Ls(1,3:Ns))/(2.0+gg)-0.5*dt*Vs(1,2:Ns-1).*(Ls(1,3:Ns)-Ls(1,1:Ns-2))/dR-(VDS(2:Ns-1)*dt);
        Ls(1,1)=(2.0*l_2+gg*Ls(1,1))/(2.0+gg)-VDS(1,1)*dt;
        
        kp=ceil(Tpos/dR);
        if kp>0 && kp<=Ns
            Tpos=Tpos+dt*Vs(kp);
        end
    end

    
end


end

