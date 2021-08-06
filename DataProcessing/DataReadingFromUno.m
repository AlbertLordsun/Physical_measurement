%% Pressure signal reading from Arduino UNO
% Description: combined with the programs from Arduino UNO, 
% read the single data in real time from arduino Uno board;

%% Serial configuration and initialization
clc
clear all
UNO = serial('COM3','BaudRate',9600,'Terminator','CR');
fopen(UNO);
x=[];
set(0,'defaultfigurecolor','w') ;

%% Pressure signal reading, display and saving
for i = 1:1:200  % Sampling frequency is 100Hz, here is for 100s;
    a = str2num(fgetl(UNO));
    if(size(a,2)~=5)
        continue
    end
    x = [x;a];
    p1=x(:,1); p2=x(:,2); p3=x(:,3); p4=x(:,4); p5=x(:,5);
    subplot(5,1,1); plot(p1);ylabel('P1');set(gca,'ytick',[]);
    subplot(5,1,2); plot(p2);ylabel('P2');set(gca,'ytick',[]);
    subplot(5,1,3); plot(p3);ylabel('P3');set(gca,'ytick',[]);
    subplot(5,1,4); plot(p4);ylabel('P4');set(gca,'ytick',[]);
    subplot(5,1,5); plot(p5);ylabel('P5');set(gca,'ytick',[]);
    xlabel('time');
    grid;
    drawnow;
    hold on;
end

save('pressures','x');

%% closing the serial port;
fclose(UNO);
% when the serial port is still occupied, use the following command;
delete(instrfindall);