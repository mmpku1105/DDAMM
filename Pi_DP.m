function [Pi,class,class_num,DP]=Pi_DP(P,n)

Pi=zeros(n,n);
temp=P^0;
for i=1:1:1000
    Pi=Pi+temp;
    temp=temp*P;
end
Pi=Pi/1000;

set=1:n;
class=cell(1,n);
class_num=0;
while ~isempty(set)
    class_num=class_num+1;
    maxVal=max(max(Pi(set,set)));
    [row,~]=find(Pi==maxVal,1);
    if maxVal>0.005
        class{class_num}=set(find(Pi(row,set)>0.0025));
    else
        class{class_num}=set;
    end
    set=setdiff(set,class{class_num});
end
I=eye(n);
DP=inv(I-P+Pi)-Pi;


