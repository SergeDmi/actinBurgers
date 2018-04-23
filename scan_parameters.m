%% Integ Script
Ns=300;
dR=1.0/Ns;
Rs=dR*(1:Ns);
model='M';

%dataname='mean_data_march2017.mat';
%dataname='may_LATA_control.mat';
%dataname='may_FILA_control.mat';
%dataname='june_SMIFH2_control.mat';
dataname='march_UTRO_control.mat';


%% Loading the data
% we will save the system state every 5 seconds
data=get_data(5,dataname);
cstes.Nmin=5;
cstes.Nmax=100;
cstes.Rs=Rs;
cstes.dR=dR;

%Phys pars
E0=1.0;
cstes.E0=E0;

%data.Tobs=data.Tobs(1:(cstes.Nmax+1));
data.Tobs=data.Tobs(1:(cstes.Nmax));
data.Tend=data.Tobs(end);
data.Nobs=length(data.Tobs);

par0=[0.4,2,0.995];
a=par0(1); % Elasticity
m=par0(2); % Contraction
c=par0(3); % Adhesion (epsilon)
% par0(4): Power law for elasticity
% par0(5): Power law for contraction
% par0(6): Power law for viscosity
% par0(7): elastictiy relaxation rate
% par0(8): 1 if contraction depends on filament end density (model D)
%          0 otherwise (model M)

% Preparing the constants
cstes.Vn=max(1.0,abs(m*1+a*1*(1-E0)));
cstes.dt=abs(0.4*dR/cstes.Vn);

% Getting ready to run many simulations
NRUNS=400000;
ERRS=zeros(5,NRUNS);
PARS=zeros(8,NRUNS);


n0=1
for n=n0:NRUNS
    n
    if strcmp(model,'D')
        PARS(:,n)=[5*rand 0.3+rand*5 0.2+0.8*rand ceil(rand*2) 1 ceil(rand*2) 1+20*rand 1];
    else
        PARS(:,n)=[0.0+6*rand 0.7+rand*2.5 0.2+rand*0.8 ceil(rand*2) ceil(rand*2) ceil(rand*2) 1+20*rand 0];
    end
    e=errors_params(data,cstes,PARS(:,n));
    ERRS(:,n)=e;
end
%[pars,score]=fminsearch(@(pars) error_params(data,cstes,pars),par0);
%[ sA,sE,sV,sR,dvg ] = intergrate_clean(data,cstes,pars);