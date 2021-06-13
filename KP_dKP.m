function [KP,dKP]=KP_dKP(mu,weight,n)

KP=zeros(n,n);
dKP=zeros(n,n,n,n);
I=eye(n);
P=zeros(n,n);
for i=1:1:n
   for j=(i+1):1:n
       P(i,j)=mu(j,i)*weight(i,j);
       P(j,i)=mu(i,j)*weight(j,i);
   end
end
for i=1:1:n
   P(i,i)=1-sum(P(i,:));
end
[Pi,class,class_num,DP]=Pi_DP(P,n);
dP=zeros(n,n,n,n);
for i=1:1:n
    for j=(i+1):1:n
        dP(i,i,i,j)=weight(i,j);
        dP(i,j,i,j)=-weight(i,j);
        dP(j,i,i,j)=weight(j,i);
        dP(j,j,i,j)=-weight(j,i);
    end
end
dPi=zeros(n,n,n,n);
for i=1:1:n
    for j=(i+1):1:n
        for k=1:1:class_num
            Pi_k=Pi_sub(P(class{k},class{k}),length(class{k}));
            if ismember(i,class{k})&&ismember(j,class{k})
               dPi_ij=dPi_sub(P(class{k},class{k}),dP(class{k},class{k},i,j),Pi_k,length(class{k}));
               dPi(class{k},class{k},i,j)=ones(length(class{k}),1)*dPi_ij;
            end
        end
    end
end
dDP=zeros(n,n,n,n);
for i=1:1:n
    for j=(i+1):1:n
        dDP(:,:,i,j)=inv(I-P+Pi)*(dP(:,:,i,j)-dPi(:,:,i,j))*inv(I-P+Pi)-dPi(:,:,i,j);
    end
end
   
for i=1:1:n
    KP(i,i)=Inf;
    for j=(i+1):1:n
        ei=I(:,i);
        ej=I(:,j);
        if weight(i,j)>0
            KP(i,j)=trace((ei*ej'-ei*ei'*P)*DP*DP);
            for k=1:1:n
                for l=(k+1):1:n
                    dKP(i,j,k,l)=trace(-ei*ei'*dP(:,:,k,l)*DP*DP+(ei*ej'-ei*ei'*P)*(dDP(:,:,k,l)*DP+DP*dDP(:,:,k,l)));
                end
            end
        else
            KP(i,j)=Inf;
        end
        if weight(j,i)>0
            KP(j,i)=trace((ej*ei'-ej*ej'*P)*DP*DP);
            for k=1:1:n
                for l=(k+1):1:n
                    dKP(j,i,k,l)=trace(-ej*ej'*dP(:,:,k,l)*DP*DP+(ej*ei'-ej*ej'*P)*(dDP(:,:,k,l)*DP+DP*dDP(:,:,k,l)));
                end
            end
        else
            KP(j,i)=Inf;
        end
    end
end