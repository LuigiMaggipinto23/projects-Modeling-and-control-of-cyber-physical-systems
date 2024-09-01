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

P=are(Ai,Bi*R^-1*Bi', Q);
K=R^-1*Bi'*P;

x0=[3,1];
Di=zeros(1,2);

K0=place(Ai,Bi,[-2j 2j]);
%K0=acker(Ai,Bi,[0 0]);
A0=Ai-Bi*K0;
Ai=A0;

F1=place(Ai',-c*Ci',[-1 -2])';

s=tf('s');
H=inv(-(s*eye(2)-(A0+c*F1*Ci)))*c*F1;
bode(H(1));

sim('local_estim.slx');
sim("local_noise.slx");