clear all
close all
clc

n=10;
q=20;
h=2;
epsilon = 10^(-8);
sigma=10^(-2);
delta=10^(-12);
EA=0;
RAD=0;
x_a_plot=zeros(q*100,1);
x_a_t_plot=zeros(q*100,1);
i_plot=0;

for j=1:100 % The experiment is repeated 100 times
    
    %CASE 1 Unaware attack
    C=randn(q,n);
    x=randn(n,1);
    I = eye(q);
    G = [ C I ];
    tau = (norm(G,2))^(-2) - epsilon;
    lambda=2*10^(-3)/tau;

    LAMBDA=lambda*[zeros(n,1);ones(q,1)];
    %eta=0; %case 1.1 with no noise
    eta=sigma*randn(q,1); %case 1.2 with noise

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

    x_0 = zeros(n,1);
    a_0 = zeros(q,1);
    x_a_0 = [x_0 ; a_0];
   

    x_a = [x ; a]; %x_a e x_t defined in this way
    y = G * x_a + eta;

    x_a_t = zeros(n+q, 1);
    x_a_t1 = zeros(n+q, 1); 
    %% algoritmo

    for t = 1:intmax
        if t==1
          x_a_t=x_a_0;
        elseif norm((x_a_t1-x_a_t), 2)>delta
            x_a_t=x_a_t1;
        else 
            break
        end
        v=x_a_t+tau*G'*(y-G*x_a_t); 
        for i=1:(n+q)
            if v(i)>tau*LAMBDA(i)
                x_a_t1(i)= v(i)-tau*LAMBDA(i);
            elseif v(i)<-tau*LAMBDA(i)
                x_a_t1(i)= v(i)+tau*LAMBDA(i);
            else
                x_a_t1(i)=0;
            end
        end
    end

    %%check di a
    controllo(:,1)=x_a;
    controllo(:,2)=x_a_t;

    x_a_plot(i_plot+1:q+i_plot)=x_a_plot(i_plot+1:q+i_plot)+x_a(n+1:end);
    x_a_t_plot(i_plot+1:q+i_plot)=x_a_t_plot(i_plot+1:q+i_plot)+x_a_t(n+1:end);
    i_plot=i_plot+q;

    RAD=RAD+calcolo_supporto(x_a(n+1:end), x_a_t(n+1:end), 0.3); 
    %Times when the support of a is correctly estimated. To do this
    %automatically, we use a function, called calcolo_supporto, in which it
    %returns 1 if the real a is correctly estimated (with an uncertainty of
    %0.3), otherwise it returns 0.

    norma = (norm(x_a(1:n) - x_a_t(1:n), 2))^2;
    EA=EA+norma;
    %The variable "EA" aim to estimate the accuracy of x, calculating
    %the l2 norm of each attempt, and in the end the global result is divided by 
    %100 to have the mean result.
        
end

RAD %Rate of attack detection

EA=EA/100 %Estimation accuracy

axis([0 2000 -3 3]);
hold on
plot(x_a_plot,'-bs','LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])
plot(x_a_t_plot,'-rs','LineWidth',1,'MarkerEdgeColor','r','MarkerFaceColor',[0.5,0.5,0.5])
title('CASE 1.2','FontSize',50);
legend({'Real attack','Estimated attack'},'FontSize',20);
xlabel(['RAD: ',num2str(RAD),'  EA: ',num2str(EA)],'FontSize',40);
