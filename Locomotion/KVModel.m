% kinematic model of purely Kelvin-voigt model for peristalsis motion
% case: five-segment larval structure
% Mechanical model:    
%	Kelvin-voigt model (spring and dashpot in parallel)
%   Mass is represented for the cuticle circle

clc;
clear all;
% configure for the model
% refer: mechanics of exploration in drosophila.M
% refer website: https://www.biorxiv.org/content/biorxiv/suppl/2018/06/24/354795.DC1/354795-8.pdf
m = 1; % standard mass
l = 1;
k = (2*pi)^2; 

% solution for the kinematic model
% p0 = [1;2;3;4;5;6;1;0;0;0;0;1];   % position and velocity
p0 = [1;2;2;1;3;1;4;1;5;1;6;0];
[t,y] = ode45(@eq,[0,10],p0);
figure(1)
% z = y(:,1:6);
z = [y(:,1),y(:,3),y(:,5),y(:,7),y(:,9),y(:,11)];
plot(t,z);

% plot for the kinematic energy
z2 = [y(:,2) y(:,4) y(:,6) y(:,8) y(:,10) y(:,12)];% cuticle velocity
K = 1/2*m*z2.^2;
figure(2)
plot(t,K);
ylabel('kinamatic energy');
xlabel('time/s');legend('1','2','3','4','5','6');

% plot for the elastic energy
z3 = [y(:,2)-y(:,1) y(:,3)-y(:,2) y(:,4)-y(:,3) y(:,5)-y(:,4) y(:,6)-y(:,5)]; % segment strain
E = 1/2*k*z3.^2;
figure(3)
plot(t,E);
ylabel('Elastic energy');
xlabel('time/s');legend('1','2','3','4','5');
 
% Newton's second law for the cuticle masses
function dp = eq(~,p)
% rule: transform the second-order equation into two one-order equation
% differential equation transformation: p--vector of positions and velocities
k = (2*pi)^2; 
n = 2*pi;
% dp = [     % velocity and acceleration 
%      p(7);p(8);p(9);p(10);p(11);p(12);
%     (k*p(2)+k*p(6)-2*k*p(1))+(n*p(8)+n*p(12)-2*n*p(7));
%     (k*p(1)+k*p(3)-2*k*p(2))+(n*p(7)+n*p(9)-2*n*p(8));
%     (k*p(2)+k*p(4)-2*k*p(3))+(n*p(8)+n*p(10)-2*n*p(9));
%     (k*p(3)+k*p(5)-2*k*p(4))+(n*p(9)+n*p(11)-2*n*p(10));
%     (k*p(4)+k*p(6)-2*k*p(5))+(n*p(10)+n*p(12)-2*n*p(11));
%     (k*p(5)+k*p(1)-2*k*p(6))+(n*p(11)+n*p(7)-2*n*p(12));
%     ];
dp = [     % velocity and acceleration 
     p(2);    (k*p(3)+k*p(11)-2*k*p(1))+(n*p(4)+n*p(12)-2*n*p(2));
     p(4);    (k*p(1)+k*p(5)-2*k*p(3))+(n*p(2)+n*p(6)-2*n*p(4));
     p(6);    (k*p(3)+k*p(7)-2*k*p(5))+(n*p(4)+n*p(8)-2*n*p(6));
     p(8);    (k*p(5)+k*p(9)-2*k*p(7))+(n*p(6)+n*p(10)-2*n*p(8));
     p(10);   (k*p(7)+k*p(11)-2*k*p(9))+(n*p(8)+n*p(12)-2*n*p(10));
     p(2);   (k*p(9)+k*p(1)-2*k*p(11))+(n*p(10)+n*p(2)-2*n*p(12));
    ];

end


