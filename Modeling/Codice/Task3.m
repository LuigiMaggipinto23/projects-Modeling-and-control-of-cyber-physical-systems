clear all
close all
clc

D = -[46 58 76 69 67 80 61;
      54 49 72 63 56 65 59;
      52 50 56 58 58 62 42;
      61 62 49 61 60 65 44;
      73 78 65 69 55 57 61;
      72 65 69 47 53 44 62];

y= -[62 58 41 46 64 63]';

q=6;
p=7;


%D=normalize(D);


tau=norm(D,2)^(-2)-10^(-8);
lambda = 10^(-3)/tau;
%lambda=1/tau;
LAMBDA=lambda*(ones(p,1));
x_t1 = zeros(p,1);
x_0 = zeros(p,1);
delta=10^(-12);

%% algoritmo

 for t = 1:intmax
    if t==1
      x_t=x_0;
    elseif norm((x_t1-x_t), 2)>delta
        x_t=x_t1;
    else 
        break
    end
    v=x_t+tau*D'*(y-D*x_t); 
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
 
 
 
 x_t
 %for i=1:p
  %   norm(D(:,i)-y,2)
 %end
 