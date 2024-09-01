clc
close all
clear all

load ('stochasticmatrices.mat');

n=10;
q=20;
h=2;

ok_vec=zeros(1,4);
norma=zeros(1,4);
iterazioni=zeros(1,4);

for indice=1:20
    C=randn(q,n);
    x=randn(n,1);
    
    a=zeros(q,1);
    S=randi(q,2,1);             %scelta randomica degli indici
    while S(1)==S(2)            %ripete finchè i due indici non sono diversi
        S(2)=randi(20);
    end
    idx=randi(4);                       %scelta randomica nella matrice
    valori(1,1)= 1 + (2-1)*rand();      %creazione di una matrice con le possibili combinazioni di valori sulle due componenti non nulle di a
    valori(2,1)= 1 + (2-1)*rand();
    valori(1,2)= -2 + (-1+2)*rand();
    valori(2,2)= -2 + (-1+2)*rand();
    valori(1,3)= 1 + (2-1)*rand();
    valori(2,3)= -2 + (-1+2)*rand();
    valori(1,4)= -2 + (-1+2)*rand();
    valori(2,4)= 1 + (2-1)*rand();
    a(S(1))=valori(1,idx);
    a(S(2))=valori(2,idx);
    
    sigma=10^-2;
    eta=sigma*randn(q,1);
    
    tau=0.03;
    lambda=(2*10^-5)/tau;
    LAMBDA=lambda*[zeros(n,1);ones(q,1)];
    
    y=C*x+eta+a;
    z=[x;a];
    z_0=zeros(n+q,1);
    G=[C eye(q)];
    
    T=10^5;
    delta=10^-7;
    
    z_k1=zeros(n+q,q);
    z_k=zeros(n+q,q);

    for qu=1:4
        switch qu
            case 1
                Q=Q1;
            case 2
                Q=Q2;
            case 3
                Q=Q3;
            case 4
                Q=Q4;
            otherwise 
                break;
        end
    
        for k=1:T
            for i=1:q
                s=0;
                if k==1
                    z_k(:,i)=z_0;
                else
                    z_k(:,i)=z_k1(:,i);
                end   
                for j=1:q
                    s=s+Q(i,j)*z_k1(:,j);
                end
                v=s+tau*G(i,:)'*(y(i)-G(i,:)*z_k(:,i));
                for w=1:(n+q)
                    if v(w)>tau*LAMBDA(w)
                        z_k1(w,i)= v(w)-tau*LAMBDA(w);
                    elseif v(w)<-tau*LAMBDA(w)
                        z_k1(w,i)= v(w)+tau*LAMBDA(w);
                    else
                        z_k1(w,i)=0;
                    end
                end
            end
            s=0;
            for w=1:q
                s=s+norm(z_k1(:,w)-z_k(:,w), 2);
            end
            if s<delta
                break;
            end
        end
        
        iterazioni(1,qu)=iterazioni(1,qu)+k;
        for j=1:q
            for i=1:q
                if abs(z_k1(n+i,j))<0.2
                    z_k1(n+i,j)=0;
                end
            end
         end
        %controllo=zeros(n,q);
        x_medio=zeros(n,1);
        ok=0;
        for i=1:q
            %controllo(:,i)=z_k1(1:n,i);
            ok=ok+calcolo_supporto(a, z_k1(n+1:end, i), 0.6);
            x_medio=x_medio+z_k1(1:n,i);
        end
        x_medio=x_medio/q;
        norma(1,qu)=norma(1,qu)+norm(x-x_medio, 2)^2;
        %controllo
        ok_vec(1,qu)=ok_vec(1,qu)+ok;
    end
end

iterazioni=iterazioni/20
norma=norma/20
ok_vec=ok_vec/20

eig(Q1)
eig(Q2)
eig(Q3)
eig(Q4)