clear all
close all
clc

q=10;
p=20;
k=2;
epsilon = 10^(-8);
sigma=10^(-2);
media=0;
tempi=zeros(100,1);
massimi=zeros(20,1);
minimi=zeros(20,1);
medie=zeros(20,1);
%while tau>0
   %tau = tau-0.001
    
    x_tilde=zeros(p,100);
    for j=1:100
        S=randi(20,2,1);
        while S(1)==S(2)
            S(2)=randi(20);
        end
        idx=randi(4);
        valori(1,1)= 1 + (2-1)*rand();
        valori(2,1)= 1 + (2-1)*rand();
        valori(1,2)= -2 + (-1+2)*rand();
        valori(2,2)= -2 + (-1+2)*rand();
        valori(1,3)= 1 + (2-1)*rand();
        valori(2,3)= -2 + (-1+2)*rand();
        valori(1,4)= -2 + (-1+2)*rand();
        valori(2,4)= 1 + (2-1)*rand();
        x_tilde(S(1),j)=valori(1,idx);
        x_tilde(S(2),j)=valori(2,idx);
    end
    delta=10^(-12);
        %for q=10:1:25
            ok=0;
            eta=sigma*randn(q,1);
            C=randn(q,p);
            tau = (norm(C,2))^(-2) - epsilon;        %%togliere 0.001
            lambda=10^(-2)/tau -0.5;               %%togliere 0.5
            LAMBDA=lambda*(ones(p,1));
            
            %ok_vec=zeros(100,1);
        while lambda<8
            ok=0;
            lambda=lambda+0.5
            LAMBDA=lambda*(ones(p,1));
           for j=1:100
            y=C*x_tilde(:,j) + eta;
            x_t1 = zeros(p,1);
            x_0 = zeros(p,1);
            
            %% algoritmo
            
            for t = 1:intmax
                if t==1
                  x_t=x_0;
                elseif norm((x_t1-x_t), 2)>delta
                    x_t=x_t1;
                else 
                    break
                end
                v=x_t+tau*C'*(y-C*x_t); 
                for i=1:p
                    if v(i)>tau*LAMBDA(i)
                        x_t1(i)= v(i)-tau*LAMBDA(i);
                    elseif v(i)<-tau*LAMBDA(i)
                        x_t1(i)= v(i)+tau*LAMBDA(i);
                    else
                        x_t1(i)=0;
                    end
                end
            end
            ok=ok+calcolo_supporto(x_tilde(:,j), x_t, 0.3);
            %ok_vec(j)=calcolo_supporto(x_tilde(:,j), x_t, 0.3);
           end
           ok
           %ok_vec
        end
        %%check di x
        tempi(j)=t-1;
        
    
    tempi;
    %massimi(k)=max(tempi(:,k));
    %minimi(k)=min(tempi(:,k));
    %medie(k)=mean(tempi(:,k));
    ok;
    %media=media+ok;
    %massimi
    %minimi
    %medie
%media=media/20