function [Pi]=Pi_sub(P,n)

I=eye(n);
temp_A=I-P;
temp_A(:,n)=ones(n,1);
temp_b=zeros(1,n);
temp_b(n)=1;
Pi=temp_b*inv(temp_A);


