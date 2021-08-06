% simulate the axis travelling wave and tranverse standing wave in larvae
% refer: <Mechanics of exploration in drosophila melanogaster>

clc;
clear all;
set(0,'defaultfigurecolor','w');
% first step: simulate for the aixs travelling wave for locomotion;
% parameter space
i = 1;  % i ranges between 1 to 5 for transverse wave;
N = 12;
mass = 1;
l = 1;
ka = (2*pi)^2;
w_a = sqrt(ka/mass);
w_a_i = 2*w_a*abs(sin((pi*i)/(N-1)));

% the segmental boundary displacement
% generate the forward locomotion
t = 0:0.01:5;
F= zeros(12,501);
for m =1:12
    for n = 1:501
        F(m,n) = (13-m)+cos(w_a_i*(n/100)+2*pi*1*m/(N-1));
    end
end
figure(1);
ylabel('Segmental actuation length');
plot(t,F);

% generate the backward locomotion
t = 0:0.01:5;
P = zeros(12,501);
P(:,1) = [11;10;9;8;7;6;5;4;3;2;1;0];
for m =1:12
    for n = 2:501
        P(m,n) = P(m,n-1)+(F(m,n)-(12-m))*0.01;
    end
end
figure(2);
ylabel('Segmental location');
xlabel('time/s');
plot(t,P);
 
% generate the locomotion animation
% forwared
% P1 = F(:,1); 
P1 = [11;10;9;8;7;6;5;4;3;2;1;0];
Y1 = zeros(12,1);
figure(3);
h = plot(P1,Y1,'d-');
% line(P,Y);
axis([0 16 -0.5 0.5])
axis off
set(h,'EraseMode','xor','MarkerSize',30);
while 1
    for i = 1:500
        P1 = P(:,i+1);
        set(h,'XData',P1);
        drawnow
    end
end
