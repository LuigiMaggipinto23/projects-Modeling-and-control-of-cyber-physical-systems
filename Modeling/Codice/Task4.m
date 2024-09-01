clear all
close all
clc

load('task4data.mat');

load('task4_sensors_positions.mat');

plot(sensors_pos(1,:), sensors_pos(2,:), '-go', LineStyle='none', LineWidth=2, MarkerFaceColor='g');
grid on;
hold on;

j=4;
p=100;
q=25;
h=2;
sigma=0.2;
heta=sigma*randn(q, 1);

%Unaware attacks
a=zeros(q,1); %choose randomly the 2 attacks (set to 30) in the 25 sensors
for c=1:h
    r=randi(q);
    while a(r)==30
        r=randi(q);
    end
    a(r)=30;
end

x=zeros(p,1); %choose randomly the 4 targets (set to 1) in the 100 cells
for c=1:j
    r=randi(100);
    while x(r)==1
        r=randi(100);
    end
    x(r)=1;
end

y=D*x+eta;

sensore_x=zeros(h,2);
sensore_y=zeros(h,2);
jdx=1;
for i=1:q
    if a(i)==30
        sensore_x(jdx,1)=sensors_pos(1,i);
        sensore_y(jdx,1)=sensors_pos(2,i);
        jdx=jdx+1;
    end
end
plot(sensore_x(:,1), sensore_y(:,1), '-rx', LineStyle='none', LineWidth=1, MarkerFaceColor='r');
G=[D eye(q)];
G=normalize(G);

z=[x; a];
z_0=zeros(p+q,1);
T=50;
tau= norm(G,2)^-2 - 10^-8;
LAMBDA=[10*ones(p,1); 20*ones(q,1)];
y=y+a;
z_k_hat=zeros(p+q,1);
z_k_mezzo_hat=zeros(p+q,1);

target_x=zeros(j,2);
target_y=zeros(j,2);

plot(sensore_x(:,2), sensore_y(:,2), '-ko', LineStyle='none', LineWidth=0.8, Marker='o', MarkerSize=10, MarkerEdgeColor='k');
plot(target_x(:,1), target_y(:,1),'-msquare', LineStyle='none', LineWidth=1, MarkerFaceColor='m');
plot(target_x(:,2), target_y(:,2), '-ksquare', LineStyle='none', LineWidth=0.8, Marker='square', MarkerSize=10, MarkerEdgeColor='k');
plot(target_x(:,1), target_y(:,1),'-wsquare', LineStyle='none', LineWidth=1, MarkerFaceColor='w');
plot(target_x(:,2), target_y(:,2), '-wsquare', LineStyle='none', LineWidth=0.8, Marker='square', MarkerSize=10, MarkerEdgeColor='w');
plot(sensore_x(:,2), sensore_y(:,2), '-wo', LineStyle='none', LineWidth=0.8, Marker='o', MarkerSize=10, MarkerEdgeColor='w');

for k=1:T
    if k==1
        z_k_hat=z_0;
    else 
        z_k_hat=[x_k1_hat; a_k1_hat];
    end   
    v=z_k_hat+tau*G'*(y-G*z_k_hat); 
    for i=1:(p+q)
        if v(i)>tau*LAMBDA(i)
            z_k_mezzo_hat(i)= v(i)-tau*LAMBDA(i);
        elseif v(i)<-tau*LAMBDA(i)
            z_k_mezzo_hat(i)= v(i)+tau*LAMBDA(i);
        else
            z_k_mezzo_hat(i)=0;
        end
    end
    x_k1_hat=A*z_k_mezzo_hat(1:p);
    a_k1_hat=z_k_mezzo_hat(p+1:p+q);
    x=A*x;
    y=D*x+eta+a;

    if k>1
        plot(target_x(:,1), target_y(:,1),'-wsquare', LineStyle='none', LineWidth=1, MarkerFaceColor='w');
        plot(target_x(:,2), target_y(:,2), '-wsquare', LineStyle='none', LineWidth=0.8, Marker='square', MarkerSize=10, MarkerEdgeColor='w');
        plot(sensore_x(:,2), sensore_y(:,2), '-wo', LineStyle='none', LineWidth=0.8, Marker='o', MarkerSize=10, MarkerEdgeColor='w');
    end
    jdx=1;
    idx=1;
    copia=x_k1_hat;
    %We implement 2 different function called trova_massimi and 
    % trova_massimi_abs, to find ,respectively, the max j estimated targets
    % and the max h estimated attacks, for each iteration. We need this to 
    % track the real targets and  the real attacks, and then depict it on 
    % the graph.                
    max_x=trova_massimi(copia,j); 
    %put to 1 all of the j max values of x_k_hat in terms of j estimated 
    % targets
    for i=1:p
        if x(i)==1
            target_x(jdx,1)=50+mod(i-1,10)*100;
            target_y(jdx,1)=50+floor((i-1)/10)*100;
            jdx=jdx+1;
        end
        if max_x(i)==1
            target_x(idx,2)=50+mod(i-1,10)*100;
            target_y(idx,2)=50+floor((i-1)/10)*100;
            idx=idx+1;
        end
    end
    plot(sensors_pos(1,:), sensors_pos(2,:), '-go', LineStyle='none', LineWidth=2, MarkerFaceColor='g');
    plot(target_x(:,1), target_y(:,1),'-msquare', LineStyle='none', LineWidth=1, MarkerFaceColor='m');
    plot(target_x(:,2), target_y(:,2), '-ksquare', LineStyle='none', LineWidth=0.8, Marker='square', MarkerSize=10, MarkerEdgeColor='k');
    plot(sensore_x(:,1), sensore_y(:,1), '-rx', LineStyle='none', LineWidth=1.5, Marker='x', MarkerSize=8, MarkerFaceColor='r');
    jdx=1;
    copia=a_k1_hat;
    max_a=trova_massimi_abs(copia,h);
    %put to 1 all of the h max values of a_k_hat in terms of h estimated
    %attack
    for i=1:q
        if max_a(i)==1
            sensore_x(jdx,2)=sensors_pos(1,i);
            sensore_y(jdx,2)=sensors_pos(2,i);
            jdx=jdx+1;
        end
    end
    plot(sensore_x(:,2), sensore_y(:,2), '-ko', LineStyle='none', LineWidth=0.8, Marker='o', MarkerSize=10, MarkerEdgeColor='k');
    title('Sparse Observer','FontSize',20);
    legend({'Sensor','Sensor under attack','Estimated attacks','Target','Estimated target'},'FontSize',15,'Location','eastoutside');
end
%Final refinement
for j=1:q
    if abs(a_k1_hat(j))<2
        a_k1_hat(j)=0;
    end
end



