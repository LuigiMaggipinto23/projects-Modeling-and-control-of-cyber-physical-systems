clear all 
close all

Ai = [ 0    1               %x=2x1
     880.87 0];
Bi = [ 0 ; -9.9453];        %u=1x1
Ci = [708.27 0];    
D0 = 0;

Ad=[0 0 0 0 0 0
    1 0 0 0 0 0
    0 1 0 0 0 0
    0 0 1 0 0 0
    0 0 0 1 0 0
    0 0 0 0 1 0];

D=diag([0 1 1 1 1 1]);

% Ad=[0 0 0 0 0 0
%     0 0 0 0 0 0
%     0 1 0 0 0 0
%     5 1 0 0 0 0
%     5 0 0 0 0 0
%     0 0 4 2 8 0];
% 
% D=diag([0 0 1 6 5 14]);

% Ad=[0 0 0 0 0 0
%     0 0 0 0 0 0
%     0 1 0 0 0 0
%     1 1 0 0 0 0
%     1 0 0 0 0 0
%     0 0 1 1 1 0];
% 
% D=diag([0 0 1 2 1 3]);

% Ad=[0 0 0 0 0 1
%     1 0 0 0 0 0
%     1 1 0 0 0 0
%     0 1 1 0 0 0
%     0 0 1 1 0 0
%     0 0 0 1 1 0];
% 
% D=diag([1 1 2 2 2 2]);

L=D-Ad;

G=diag([1 0 0 0 0 0]);

%G=diag([4 4 0 0 0 0]);

%G=diag([1 1 0 0 0 0]);

%G=diag([1 0 0 0 0 0]);

lambda_g=eig(L+G);

c=10/2*min(real(lambda_g));      

q=1;
Q=q*eye(2);
r=q/10;
%r=q*100;
R=r;

x0=[3,1];
Di=zeros(1,2);

K0=place(Ai,Bi,[-2j 2j]);
%K0=acker(Ai,Bi,[0 0]);
A0=Ai-Bi*K0;
Ai=A0;

P=are(Ai,Bi*R^-1*Bi', Q);
K=R^-1*Bi'*P;

P=are(Ai',Ci'*R^-1*Ci, Q);
F=P*Ci'*R^-1;

L0=place(Ai',Ci',[0 -3])';

sim('global_estim.slx');
sim("global_noise.slx");