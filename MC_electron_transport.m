clc;clear;close all

T = 300; %k temperature
K = 1.38064852e-23;% J/k boltzmann const
me= 9.10938356e-31; % kg mass of electron 
Eme= 0.26* me;% kg effective mass of electron
vth = sqrt(3*K*T/Eme);% m/s thermal volocity

N = 10; %number of electons
dimx =200;%nm
dimy =100;%nm
xpos = randi(dimx+1,1,N)-1;% nm electron position
ypos = randi(dimy+1,1,N)-1;% nm electron position
vx   = zeros(1,N);
vy   = zeros(1,N);

% positions particals untill there are no overlaps
while (length(xpos) ~= length(unique(xpos))) && (length(ypos) ~= length(unique(ypos)))
    xpos = randi(dimx+1,1,N)-1;
    ypos = randi(dimy+1,1,N)-1;
end

vx = vth* cos(rand(1,N)*2*pi);
vy = vth* sin(rand(1,N)*2*pi);


% plots the 
for i =1:N
plot(xpos(:,i),ypos(:,i),'o')
hold on
quiver(xpos(:,i),ypos(:,i),vx(i),vy(i),0.0001)

end
xlim([0,dimx])
ylim([0,dimy])
hold off