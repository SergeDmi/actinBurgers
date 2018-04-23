function [dS]=stress_grad(A,E,L,N,dR,a,m,Pa,Pm,Pl,A0)
    

    S=(m*A./(L.^Pl)).^Pm+a*(A.^Pa).*(1-E);
    
    dS=zeros(1,N);
    dS(1,2:N-1)=(S(3:N)-S(1:N-2))/(2*dR);
end