function [dPi]=dPi_sub(P,dP,Pi,n)

I=eye(n);
temp_A=I-P;
temp_A(:,n)=ones(n,1);
temp_c=Pi*dP;
temp_c(n)=0;
dPi=temp_c*inv(temp_A);