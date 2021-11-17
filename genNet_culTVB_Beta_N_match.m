%=====Simulate the relationship between TV,IV,SV and U1,U2==========
Ns=[50,100,200,500,1000]; %Number of nodes  [50,100,200,500,1000,5000]
c=0.05; %cost
ps=[0.3,0.5,0.7];  %network density

TVs=zeros(); IVs=zeros(); % where data is stored
N=Ns(3);
p=ps(2);
TVB=zeros(); IVB=zeros(); %where data is stored
I=eye(N);
one=ones(N,1);
for si=1:500        %simulate 500 groups
  
    [G,Beta]=get_G(N,p); %generate a network
    alpha=rand(N,1);
    
    key=Weighted_Bonacich_centrality((G+G')/2,(Beta)^-1,((Beta)^-1)*one);  %proposition 5
    alpha=sort_alpha(key,alpha);
    TVB(si)=TVU1Base(G,alpha,Beta,c);  %work out the result
    
    
    key=Bonacich_centrality(G',(Beta)^-1);  %proposition 4
    alpha=sort_alpha(key,alpha);
    IVB(si)=TVU1Base(G,alpha,Beta,c);  %work out the result
    
    for ti=1:1000             %match 1000 times
        key=randperm(N); %randomly arrange
        alpha=alpha(key);   %change the order


         TVs(si,ti)=TVU1Base(G,alpha,Beta,c);  %work out the result
         IVs(si,ti)=IVU2Base(G,alpha,Beta,c);  %work out the result

        fprintf('execute:(si:%d,ti=%d)\n',si,ti);
    end
end


%--------plot-------------------------------
% figure();
% boxplotUU(TVU1);
% 
% figure();
% boxplotUU(IVU2);
disp("Result:")
key=0;
for si=1:length(TVB)
    if TVB(si)>=max(TVs(si,:))
        key=key+1;
    else
%         print(key)
    end
end
fprintf('(TV:best match>=random match:%4.2f %%)\n',(key/length(TVB))*100);

key=0;
for si=1:length(IVB)
    if IVB(si)>=max(IVs(si,:))
        key=key+1;
    else
%         print(key)
    end
end
fprintf('(IV:best match>=random match:%4.2f %%)\n',(key/length(IVB))*100);



function B=Weighted_Bonacich_centrality(G,Y,w) %proposition 5
I=eye(length(G));   
B= ((I-Y*G)^-1)*w;
end


function B=Bonacich_centrality(G,Y) %proposition 4
I=eye(length(G)); 
one=ones(length(G),1);
B= ((I-Y*G)^-1)*one;
end


function [TV]=TVU1Base(G,alpha,Beta,c)   %step calculation
one=ones(length(alpha),1);
     X3=(Beta^-1)*((alpha-c*one)/2);
%     Buy=find(X3>=0);
     NoBuy=find(X3<0);  %find people who don't buy,delete them.
     
         alpha1=alpha;
         Beta1=Beta;

         alpha1(NoBuy)=[];
         Beta1(NoBuy,:)=[];Beta1(:,NoBuy)=[];  
         one=ones(length(alpha1),1);
     %calculate B3
     B3=(1/4)*((alpha1-c*one)')*(Beta1^-1)*(alpha1-c*one);
     
     NoBuy1=NoBuy;
     for i =length(NoBuy):-1:1
         temp=NoBuy;    %find people who don't buy,try to keep them one by one,delete others.
         G1=G;
         alpha1=alpha;
         Beta1=Beta;
         
         temp(i)=[];
         G1(temp,:)=[];G1(:,temp)=[];
         alpha1(temp)=[];
         Beta1(temp,:)=[];Beta1(:,temp)=[];
         one=ones(length(alpha1),1);
         
         X1=((Beta1-((G1+G1')/2))^-1)*((alpha1-c*one)/2);
         if all(X1>=0)
             NoBuy1(i)=[]; %Those who make X1>=0 can be kept. 
         end
     end
     %calculate B1
         G1=G;
         alpha1=alpha;
         Beta1=Beta;
         
         G1(NoBuy1,:)=[];G1(:,NoBuy1)=[];
         alpha1(NoBuy1)=[];
         Beta1(NoBuy1,:)=[];Beta1(:,NoBuy1)=[];  %delete costomer in NoBuy
         one=ones(length(alpha1),1);
         
         %verify whether X1>=0
         X1=((Beta1-((G1+G1')/2))^-1)*((alpha1-c*one)/2);
         if ~all(X1>=0)
             disp("Erro in TV");
             pause;
         end
     
     A=(Beta1-G1)^-1;
try
   B1=(1/2)*((alpha1-c*one)')*(((A^-1)+((A')^-1))^-1)*(alpha1-c*one);
catch 
   disp("Erro!")
end
     
     TV=B1-B3;
     
%      temp=(alpha1-c*one);
%      U1=(temp'/2)*(((A^-1)+((A')^-1))^-1)*Beta1*(((A^-1)+((A')^-1))^-1)*temp;
     
end



function [IV]=IVU2Base(G,alpha,Beta,c)   %%step calculation

    % calculate B4
    [P,PU,X]=getP4(Beta,alpha,c);
    alpha1=alpha;
    Beta1=Beta;
    while min(PU)<P
        NoBuy=find(P>PU);

         alpha1(NoBuy)=[];
         Beta1(NoBuy,:)=[];Beta1(:,NoBuy)=[];
         one=ones(length(alpha1),1);
         [P,PU,X]=getP4(Beta1,alpha1,c);
    end
    
    %verify whether X>=0 
    if ~all(X>=0)
        disp("B4 ERRO!")
        pause;
    end
    
    B4=(1/4)*((one'*(Beta1^-1)*alpha1-c*one'*(Beta1^-1)*one)^2)/(one'*(Beta1^-1)*one);
 
    % calculate B2
    G1=G;
    [P,PU,X]=getP2(G,Beta,alpha,c);
    alpha1=alpha;
    Beta1=Beta;
    while min(PU)<P  %if min(PU)<P,means there are negative values in X,we still need to delete customer
        NoBuy=find(min(PU)==PU);

         alpha1(NoBuy)=[];
         Beta1(NoBuy,:)=[];Beta1(:,NoBuy)=[];
         G1(NoBuy,:)=[];G1(:,NoBuy)=[];  %delete the customer whose PU=min(PU) each time in alpha,beta,G
         [P,PU,X]=getP2(G1,Beta1,alpha1,c);  %update P,PU,X
    end
    
    %verify whether X>=0
    if ~all(X>=0)
        disp("B2 ERRO!")
        pause;
    end
    A=(Beta1-G1)^-1;
    one=ones(length(alpha1),1);
try
   B2=(1/4)*((one'*A*alpha1-c*one'*A*one)^2)/(one'*A*one);
catch 
   disp("Erro!")
end
    

    IV=B2-B4;
     
% Delta=(1/2)*(((one')*A*alpha1)/((one')*A*one));
% U2=(1/2)*(alpha1-c*0.5*one-Delta*one)'*(A')*Beta1*A*(alpha1-c*0.5*one-Delta*one);

end

function [P,PU,X]=getP4(Beta,alpha,c) %get P,PU,X in case 4
    one=ones(length(alpha),1);
    P=(1/2)*((one')*(Beta^-1)*(alpha+c*one))/((one')*(Beta^-1)*one);
    PU=((Beta^-1)*(alpha))./((Beta^-1)*one);
    
    X=(Beta^-1)*alpha-P*(Beta^-1)*one;
end 

function [P,PU,X]=getP2(G,Beta,alpha,c) %get P,PU,X in case 2
    one=ones(length(alpha),1);
    A=(Beta-G)^-1;
    P=(1/2)*((one')*A*(alpha+c*one))/((one')*A*one);
    PU=(A*(alpha))./(A*one);
    
    X=A*alpha-P*A*one;
end 

function boxplotUU(Data) %plot
[Ns,ps,si]=size(Data);
datas=zeros();
for i=1:Ns
    for j=1:ps
       datas(1:si,(i-1)*ps+j)=Data(i,j,1:si);
    end
end
boxplot(datas);
end


function [G,Beta]=get_G(N,p)
%------generate a network-----------------------------------------
I=eye(N);
G=zeros(N);
for i=1:N-1
    for j=i+1:N
        if rand(1)<p  % limit the number of edges
            G(i,j)=1;
            G(j,i)=1;
        end
    end
end
% Beta=sum(G,2)+(N*0.1);  %sum of the rows of G plus a little
Beta=N;
Beta=Beta.*I;
end


function temp=sort_alpha(key,alpha) %sort alpha
temp=zeros(length(alpha),1);
    [~,index1] = sort(key);
    [~,index2] =sort(alpha);
    for i =1:length(alpha)
        temp(index1(i))=alpha(index2(i));
    end

end
