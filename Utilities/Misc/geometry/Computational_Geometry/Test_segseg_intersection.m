close all
clear all
clc

global IAXIS JAXIS NOD1 NOD2 EDG1 EDG2 GEOMEPS

IAXIS = 1;
JAXIS = 2;
NOD1  = 1;
NOD2  = 2;
EDG1  = 1;
EDG2  = 2;
GEOMEPS=1.e-12;

P1 = [0.5 0.5]; D1 = [.5 .5]; T1 = D1/norm(D1);
P2 = [.85 .07]; D2 = [-1. 3.]; T2 = D2/norm(D2);

figure
hold on
plot([P1(1) P1(1)+D1(1)],[P1(2) P1(2)+D1(2)],'-xk','Linewidth',4)
plot([P2(1) P2(1)+D2(1)],[P2(2) P2(2)+D2(2)],'-ob','Linewidth',2)
grid on
axis equal 
box on

[SVARV,SLENV,INT_FLG]=GET_SEGSEG_INTERSECTION(P1,D1,P2,D2)

if(INT_FLG == 1)
plot([P1(1)+SVARV(1,1)*T1(1)],[P1(2)+SVARV(1,1)*T1(2)],'-xk','MarkerSize',12)
plot([P2(1)+SVARV(1,2)*T2(1)],[P2(2)+SVARV(1,2)*T2(2)],'-ob','MarkerSize',12)
end

if(INT_FLG == 2)
plot([P1(1)+SVARV(1,1)*T1(1)],[P1(2)+SVARV(1,1)*T1(2)],'-xk','MarkerSize',12)
plot([P2(1)+SVARV(1,2)*T2(1)],[P2(2)+SVARV(1,2)*T2(2)],'-ob','MarkerSize',12)

plot([P1(1)+SVARV(2,1)*T1(1)],[P1(2)+SVARV(2,1)*T1(2)],'-xk','MarkerSize',12)
plot([P2(1)+SVARV(2,2)*T2(1)],[P2(2)+SVARV(2,2)*T2(2)],'-ob','MarkerSize',12)
end
    

