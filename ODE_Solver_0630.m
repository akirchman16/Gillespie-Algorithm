clear all;

% This code will plot solutions to the ODEs that describe deterministic
% kinetics. Equations are based on Joseph's kinetics powerpoint. An Euler
% method is used to solve/approximate ODEs.

k_f = 3;    %reaction rate constants
k_r = 1;

nA(1) = 10; %initial conditions
nB(1) = 0;

dt = 0.01;  %time things
t(1) = 0;

for i = 1:1000
    matrix = [1-(k_f*dt) k_r*dt; k_f*dt 1-(k_r*dt)];
    population = matrix*[nA(i);nB(i)];
    nA(i+1) = population(1,1);
    nB(i+1) = population(2,1);
    t(i+1) = t(i)+dt;
end

figure();
plot(t,nA,'r');
hold on;
plot(t,nB,'black');
xlabel('Time');
ylabel('Species Populations');
legend('nA','nB');
