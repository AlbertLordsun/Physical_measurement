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
% generate the backward locomotion
t = 0:0.01:5;
B= zeros(12,501);
for m =1:12
    for n = 1:501
        B(m,n) = (13-m)+cos(w_a_i*(n/100)-2*pi*1*m/(N-1));
    end
end
figure(1);
plot(t,B);

% generate the locomotion animation
% P2 = B(:,1); % backward
P2 = [11;10;9;8;7;6;5;4;3;2;1;0];

Y2 = zeros(12,1);
figure(2);
h = plot(P2,Y2,'.-');
axis([-20 20 -1 1])
axis off
set(h,'EraseMode','xor','MarkerSize',30);
while 1
    for i = 1:500
    drawnow
    for j=1:12
        P2(j) = P2(j)+ (B(j,i+1)-(14-j))*0.01;
    end
    set(h,'XData',P2);
    end
end