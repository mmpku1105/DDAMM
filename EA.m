function [CS,class,Best_set,Best_true]=EA(CS,n,p,weight,T,alpha)

T0=0;
N=zeros(n,n);
mu=zeros(n,n);
y=zeros(n,n);
dpi=zeros(n,n,n);
for i=1:1:n
    for j=i+1:1:n
        N(i,j)=2*T0/n/(n-1);
        for k=1:1:N(i,j)
            temp=unifrnd(0,1,1);
            if temp<p(i,j)
                y(i,j)=y(i,j)+1;
            end
        end
        mu(i,j)=(y(i,j)+1)/(N(i,j)+2);
        mu(j,i)=1-mu(i,j);
    end
end

set=zeros(1,n);
for i=1:1:n
    set(i)=i;
end

t=0;
while t<T
    for i=1:1:n
        for j=(i+1):1:n
            if weight(i,j)>0
                t=t+1;
                N(i,j)=N(i,j)+1;
                temp=unifrnd(0,1,1);
                if temp<p(i,j)
                   y(i,j)=y(i,j)+1;
                end
                mu(i,j)=(y(i,j)+1)/(N(i,j)+2);
                mu(j,i)=1-mu(i,j);
            end
        end
    end
    G=zeros(n,n);
    for i=1:1:n
       for j=(i+1):1:n
           G(i,j)=mu(j,i)*weight(i,j);
           G(j,i)=mu(i,j)*weight(j,i);
       end
    end
    for i=1:1:n
       G(i,i)=1-sum(G(i,:));
    end
    [~,class,class_num,~]=Pi_DP(G,n);
    Best_set=zeros(1,class_num);
    Pi_n=zeros(1,n);
    for k=1:1:class_num
       Pi_n(class{k})=Pi_sub(G(class{k},class{k}),length(class{k}));
    end    
    for i=1:1:class_num
        [~,index]=max(Pi_n(class{i}));
        Best_set(i)=class{i}(index);
    end
    
    [KP,dKP]=KP_dKP(mu,weight,n);
    minVal=min(min(KP));
    [row,col]=find(KP==minVal,1);
    M=sum(sum(weight~=0));
    mark=0;
    for i=1:1:n
        for j=(i+1):1:n
            if weight(i,j)>0
                if (row~=i)||(col~=j)
                    tau=0;
                    for k=1:1:n
                        for l=(k+1):1:n
                            tau=tau+(dKP(row,col,k,l)-dKP(i,j,k,l))^2*(y(k,l)+1)*(N(k,l)+1-y(k,l))/(N(k,l)+2)/(N(k,l)+2)/(N(k,l)+3);
                        end
                    end
                    if normcdf((KP(row,col)-KP(i,j))/sqrt(tau))>alpha/(M-1)
                        mark=1;
                        break;
                    end
                end
            end
            if weight(j,i)>0
                if (row~=j)||(col~=i)
                    tau=0;
                    for k=1:1:n
                        for l=(k+1):1:n
                            tau=tau+(dKP(row,col,k,l)-dKP(j,i,k,l))^2*(y(k,l)+1)*(N(k,l)+1-y(k,l))/(N(k,l)+2)/(N(k,l)+2)/(N(k,l)+3);
                        end
                    end
                    if normcdf((KP(row,col)-KP(j,i))/sqrt(tau))>alpha/(M-1)
                        mark=1;
                        break;
                    end                    
                end
            end
        end
        if mark==1
            break;
        end
    end
    if (mark==0)||(class_num<=T/500)
       weight(row,col)=0;
       weight(col,row)=0;
    end
end

G=zeros(n,n);
for i=1:1:n
   for j=(i+1):1:n
       G(i,j)=p(j,i)*weight(i,j);
       G(j,i)=p(i,j)*weight(j,i);
   end
end
for i=1:1:n
   G(i,i)=1-sum(G(i,:));
end
[~,class,class_num,~]=Pi_DP(G,n);
Best_true=zeros(1,class_num);
Pi_true=zeros(1,n);
for k=1:1:class_num
   Pi_true(class{k})=Pi_sub(G(class{k},class{k}),length(class{k}));
end  
for i=1:1:class_num
    [~,index]=max(Pi_true(class{i}));
    Best_true(i)=class{i}(index);
end
if length(Best_set)==length(Best_true)
    if all(Best_set==Best_true)
       CS=CS+1;
    end
end
