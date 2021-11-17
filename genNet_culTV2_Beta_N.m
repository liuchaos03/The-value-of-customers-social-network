%=====Simulate the relationship between TV,IV,SV and U1,U2==========
Ns=[50,100,200,500,1000]; %Number of nodes [50,100,200,500,1000]
c=0.05; %cost
ps=[0.3,0.5,0.7];  %network density

TVU1=zeros(); IVU2=zeros(); % where data is stored
for i = 1:length(Ns)  % traverse every number of nodes
    N=Ns(i);
    for j=1:length(ps)  % traverse every network density
        p=ps(j);
        
        for si=1:500          %simulate 500 groups
            TVs=zeros(); IVs=zeros(); U1s=zeros(); U2s=zeros();
            [G,alpha,Beta]=test(N,p,c);   %Find a workable G and Alpha
            for ti=1:100            %random match 100 times
                key=randperm(N); %randomly arrange
                alpha=alpha(key);   %change the order
                Beta=sum(Beta,2);      
                Beta=Beta(key);
                Beta=Beta.*eye(N);
             
                key=Assumption(G,alpha,Beta,c);  %judge whether the conditions are met
                if all(key==1)
                    
                    [TV,IV,SV,U1,U2,X1]=Profit_utility(G,alpha,Beta,c);  %work out the result
                    TVs(ti)=TV;  %save the result
                    IVs(ti)=IV; 
                    U1s(ti)=U1; 
                    U2s(ti)=U2;
                else
                    disp("Erro!")
                    pause; %terminate the program
                end
                fprintf('execute:(N:%d,p=%d,si=%d,ti=%d)\n',i,j,si,ti);
            end
            TVU1(i,j,si)=corr(TVs',U1s','type','Kendall');  %calculate the similarity
            IVU2(i,j,si)=corr(IVs',U2s','type','Kendall');
        end
    end
end

%--------plot-------------------------------
figure();
boxplotUU(TVU1);

figure();
boxplotUU(IVU2);




function boxplotUU(Data)  %plot
[Ns,ps,si]=size(Data);
datas=zeros();
for i=1:Ns
    for j=1:ps
       datas(1:si,(i-1)*ps+j)=Data(i,j,1:si);
    end
end
boxplot(datas);
end



function [G,alpha,Beta]=test(N,p,c)  
%--------
[G,Beta]=get_G(N,p);  %generate a network
ti=0;
key=0;
while ~all(key==1)
    alpha=customer(N,c);      %generate Alpha
    key=Assumption(G,alpha,Beta,c);  %judge whether the conditions are met
    ti=ti+1;
    if ti>N
        disp(["find:",string(ti)])
        pause;
    end
end
% 
end

function [TV,IV,U1,U2]=tau_si(N,p,c)
%--------
[G,Beta]=get_G(N,p);  %generate a network
ti=0;
key=0;
while ~all(key==1)
    alpha=customer(N,c);      %generate Alpha
    key=Assumption(G,alpha,Beta,c);  %judge whether the conditions are met
    ti=ti+1;
end
disp(["execute:",string(ti)])

[TV,IV,SV,U1,U2,X1]=Profit_utility(G,alpha,Beta,c);  %work out the result
if all(X1<0)
    disp("Erro!")
end

end

function alpha=customer(N,c)
 alpha=rand(N,1)*0.5+c*10;
end

function [TV,IV,SV,U1,U2,X1]=Profit_utility(G,alpha,Beta,c)
N=length(G);
one=ones(N,1);
I=eye(N);
A=(Beta-G)^-1;
%------calculate profit and utility--------------------------------------------
temp=(alpha-c*one);
a=Beta^-1;

TV=(temp'/2)*(((Beta-((G+(G'))/2))^-1)-a)*(temp/2);
IV=(temp'/2)*((((A')*one*(one')*A)/((one')*A*one))-(((a')*one*(one')*a)/((one')*a*one)))*(temp/2);
SV=(temp'/2)*(FM(A)-FM(a))*(temp/2);

U1=(temp'/2)*(((A^-1)+((A')^-1))^-1)*Beta*(((A^-1)+((A')^-1))^-1)*temp;
Delta=(1/2)*(((one')*A*alpha)/((one')*A*one));
U2=(1/2)*(alpha-c*0.5*one-Delta*one)'*(A')*Beta*A*(alpha-c*0.5*one-Delta*one);

% P1=((A+A')^-1)*(A*alpha+(A')*one*c);
X1=(((A^-1)+((A^-1)'))^-1)*temp;
% P2=((one')*A*alpha)/(2*(one')*A*one)+(0.5*c);
% X2=A*(alpha-P2*one);
% P3=0.5*(alpha+c*one);
% X3=0.5*(Beta^-1)*(alpha-c*one);
% P4=((one')*a*alpha)/(2*(one')*a*one)+(0.5*c);
% X4=a*(alpha-P4*one);

end


function key=Assumption(G,alpha,Beta,c)
N=length(G);
one=ones(N,1);
%------judge whether the conditions are met-------------------------------------
key=[0,0,0,0];   %key is used to represent the satisfaction of the four inequalities
%Assumption 1.
if sum(sum(Beta,2)>sum(G,2))==N
    key(1)=1;
else
    disp("Erro in As1")
end
    
%Assumption 2.
if min(alpha)>=c
    key(2)=1;
else
    disp("Erro in As2")
end

%Assumption 3.
A=(Beta-G)^-1;
C=(A*alpha)./(A*one);
condition2=(1/2)*((one')*A*(alpha+c*one))/((one')*A*one);
if condition2<=min(C)
    key(3)=1;
else
     disp("Erro in As3-1")
end
condition3=(1/2)*((one')*(Beta^-1)*(alpha+c*one))/((one')*(Beta^-1)*one);
if condition3<=min(((Beta^-1)*alpha)./((Beta^-1)*one))
    key(4)=1;
else
   disp("Erro in As3-2")
end

end


function [G,Beta]=get_G(N,p)
%------generate a network------------------------------------------
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
Beta=N+N*rand(N,1);
Beta=Beta.*I;
end

function f=FM(M)
one=ones(length(M),1);
  f=((((M^-1)+((M')^-1))/2)^-1)-((M')*one*(one')*M)/((one')*M*one);
end

