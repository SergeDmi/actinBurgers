function [ alphs,mold ] = get_alphas( data,cstes,sA,sV,sT )
Nmin=cstes.Nmin;
Nmax=cstes.Nmax;
sS=size(sA);
rr=cstes.Rs-cstes.dR/2;
Vr=4*pi*cstes.dR*rr.^2;
Vr(1)=(4.0/3.0)*pi*cstes.dR.^3;
%figure
%hold all
Nmax=min(Nmax,size(sV,1));
for t=1:Nmax 
    p=floor(sS(2)*sT(t));
    p=min(p,size(sV,2));
    PTS=[1:p;sV(t,1:p)];
    [ slope,offset] = ortho_robust_coeff( PTS );
    alphs(t)=-slope*data.K0*sS(2);
    sp(t)=p;
    mold(t)=sum(Vr(1:p).*sA(t,1:p));
    
end
mold=mold/max(mold);

end