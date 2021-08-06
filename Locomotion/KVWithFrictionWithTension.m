% kinematic model of Kelvin-voigt model with friction with tension for peristalsis motion
% case: five-segment larval structure
% Mechanical model: scaled model of Newton's second law
%	Kelvin-voigt model (spring and dashpot in parallel)
%   Mass is represented for the cuticle circle
%   friction -- dependent on the motion states and direction
%   tension -- considering the power flow in each segments

clc;
clear all;
% configure for the model
% refer: mechanics of exploration in drosophila.M
% refer website: https://www.biorxiv.org/content/biorxiv/suppl/2018/06/24/354795.DC1/354795-8.pdf
m = 1;
l = 1;
k = (2*pi)^2; 
n = 2*pi;
mu = 1; % friction coefficient;

% solution for the kinematic model
p0 = [1;1;2;0;3;0;4;0;5;0;6;1];
[t,y] = ode45(@eq,[0,10],p0);
z1 = [y(:,1) y(:,3) y(:,5) y(:,7) y(:,9) y(:,11)]; % cuticle position

% plot for the cuticle position change
figure(1)
plot(t,z1);
ylabel('Cuticle position');
xlabel('time/s');
legend('1','2','3','4','5','6');

% plot for the cuticle kinematic energy
z2 = [y(:,2) y(:,4) y(:,6) y(:,8) y(:,10) y(:,12)]; % cuticle velocity
K = 1/2*m*z2.^2;
figure(2)
plot(t,K);
ylabel('Kinamatic energy');
xlabel('time/s');
legend('1','2','3','4','5','6');

% plot for the elastic energy
z3 = [y(:,2)-y(:,1) y(:,3)-y(:,2) y(:,4)-y(:,3) y(:,5)-y(:,4) y(:,6)-y(:,5)]; % segment strain
E = 1/2*k*z3.^2;
figure(3)
plot(t,E);
ylabel('Elastic energy');
xlabel('time/s');
legend('1','2','3','4','5');

% plot for the tension change
D1 = [1 0 0 0 0 -1;-1 1 0 0 0 0; 0 -1 1 0 0 0; 0 0 -1 1 0 0; 0 0 0 -1 1 0; 0 0 0 0 -1 1];
beta = 0.1;
u = beta* D1'*z2';
z4 = u';
figure(4)
plot(t,z4);
ylabel('Tension');
xlabel('time/s');
legend('1','2','3','4','5','6');

% function of Newton's second law for the cuticle masses
function dp = eq(~,p)
gamma_f = 0.9;
gamma_k = 0.1;
gamma = 0.9;
zeta = 0;
beta = 0.1;
f = zeros(6,1);
% equation according to the pdf file
dp = [ % velocity and acceleration
    p(2);(p(3)+p(11)-2*p(1))+((2*zeta-beta)*p(4)+(2*zeta-beta)*p(12)-2*(2*zeta-beta)*p(2))+f(1);
    p(4);(p(1)+p(5)-2*p(3))+((2*zeta-beta)*p(2)+(2*zeta-beta)*p(6)-2*(2*zeta-beta)*p(4))+f(2);
    p(6);(p(3)+p(7)-2*p(5))+((2*zeta-beta)*p(4)+(2*zeta-beta)*p(8)-2*(2*zeta-beta)*p(6))+f(3);
    p(8);(p(5)+p(9)-2*p(7))+((2*zeta-beta)*p(6)+(2*zeta-beta)*p(10)-2*(2*zeta-beta)*p(8))+f(4);
    p(10);(p(7)+p(11)-2*p(9))+((2*zeta-beta)*p(8)+(2*zeta-beta)*p(12)-2*(2*zeta-beta)*p(10))+f(5);
    p(2);(p(9)+p(1)-2*p(11))+((2*zeta-beta)*p(10)+(2*zeta-beta)*p(2)-2*(2*zeta-beta)*p(12))+f(6);     
    ];
% the tension is calculated by the difference in the position
% dp = [ % velocity and acceleration
%     p(2);((1-beta)*p(3)+(1-beta)*p(11)-(2-beta)*p(1))+((2*zeta)*p(4)+(2*zeta)*p(12)-2*(2*zeta)*p(2))+f(1);
%     p(4);((1-beta)*p(1)+(1-beta)*p(5)-(2-beta)*p(3))+((2*zeta)*p(2)+(2*zeta)*p(6)-2*(2*zeta)*p(4))+f(2);
%     p(6);((1-beta)*p(3)+(1-beta)*p(7)-(2-beta)*p(5))+((2*zeta)*p(4)+(2*zeta)*p(8)-2*(2*zeta)*p(6))+f(3);
%     p(8);((1-beta)*p(5)+(1-beta)*p(9)-(2-beta)*p(7))+((2*zeta)*p(6)+(2*zeta)*p(10)-2*(2*zeta)*p(8))+f(4);
%     p(10);((1-beta)*p(7)+(1-beta)*p(11)-(2-beta)*p(9))+((2*zeta)*p(8)+(2*zeta)*p(12)-2*(2*zeta)*p(10))+f(5);
%     p(2);((1-beta)*p(9)+(1-beta)*p(1)-(2-beta)*p(11))+((2*zeta)*p(10)+(2*zeta)*p(2)-2*(2*zeta)*p(12))+f(6);     
%     ];
% revised version1 of newton second law
% dp = [ % velocity and acceleration
%     p(2);((1-beta)*p(3)+(1-beta)*p(9)-(1-beta)*p(11)-(1-beta)*p(1))+((2*zeta)*p(4)+(2*zeta)*p(12)-2*(2*zeta)*p(2))+f(1);
%     p(4);((1-beta)*p(1)+(1-beta)*p(5)-(2-beta)*p(3))+((2*zeta)*p(2)+(2*zeta)*p(6)-2*(2*zeta)*p(4))+f(2);
%     p(6);((1-beta)*p(3)+(1-beta)*p(7)-(2-beta)*p(5))+((2*zeta)*p(4)+(2*zeta)*p(8)-2*(2*zeta)*p(6))+f(3);
%     p(8);((1-beta)*p(5)+(1-beta)*p(9)-(2-beta)*p(7))+((2*zeta)*p(6)+(2*zeta)*p(10)-2*(2*zeta)*p(8))+f(4);
%     p(10);((1-beta)*p(7)+(1-beta)*p(11)-(2-beta)*p(9))+((2*zeta)*p(8)+(2*zeta)*p(12)-2*(2*zeta)*p(10))+f(5);
%     p(2);((1-beta)*p(3)+(1-beta)*p(9)-(1-beta)*p(1)-(1-beta)*p(11))+((2*zeta)*p(10)+(2*zeta)*p(2)-2*(2*zeta)*p(12))+f(6);     
%     ];
% % revised version2 of newton second law
% dp = [ % velocity and acceleration
%     p(2);(p(3)+p(9)-p(11)-p(1))+((2*zeta-beta)*p(4)+(2*zeta-beta)*p(12)-(2*zeta-beta)*p(2)-(2*zeta-beta)*p(12))+f(1);
%     p(4);(p(1)+p(5)-2*p(3))+((2*zeta-beta)*p(2)+(2*zeta-beta)*p(6)-2*(2*zeta-beta)*p(4))+f(2);
%     p(6);(p(3)+p(7)-2*p(5))+((2*zeta-beta)*p(4)+(2*zeta-beta)*p(8)-2*(2*zeta-beta)*p(6))+f(3);
%     p(8);(p(5)+p(9)-2*p(7))+((2*zeta-beta)*p(6)+(2*zeta-beta)*p(10)-2*(2*zeta-beta)*p(8))+f(4);
%     p(10);(p(7)+p(11)-2*p(9))+((2*zeta-beta)*p(8)+(2*zeta-beta)*p(12)-2*(2*zeta-beta)*p(10))+f(5);
%     p(2);(p(3)+p(9)-p(1)-p(11))+((2*zeta-beta)*p(10)+(2*zeta-beta)*p(2)-(2*zeta-beta)*p(2)-(2*zeta-beta)*p(12))+f(6);     
%     ];
for i = 1:6
    if  p(2*i)>0
        f(i) = -gamma_f*gamma_k*gamma;
    elseif p(2*i)<0
        f(i) = gamma_k*gamma;
    else
        dp(2*i) = 0;
    end
end
end