clear all;
close all;

% This code is to test the effect of increasing the number of molecules
% in the Gillespie algorithm and comparing it to the ODE model. This
% reaction first will be solved for A<-->B

t(1) = 0;   %starts time at zero

k_f = 2;    %rate constant of forward reaction (A-->B)
k_r = 1;    %rate constant of reverse reaction (B-->A)
nA(1) = 100; %initial number of A molecules
nB(1) = 0; %initial number of B molecultes
Iterations = 5*(nA(1)+nB(1));   %number of times the algorithm will run

ODEnA(1) = nA(1);
ODEnB(1) = nB(1);

for i = 1:Iterations
    R_tot(i) = (k_f.*nA(i))+(k_r.*nB(i));  %calculates the total reaction rate                   
    dt(i) = random('Exponential',1/R_tot(i)); %selects random value for time interval (stores history in array)
    
    P(i) = (k_f.*nA(i))./(R_tot(i));  %probability of A becoming B
    if rand <= P(i) %tests probability
        nA(i+1) = nA(i)-1;  %decreases count of present A molecules
        nB(i+1) = nB(i)+1;  %increases count of present B molecules
    else           %otherwise a reverse reaction has to occur
        nA(i+1) = nA(i)+1;  %increases count of present A molecules
        nB(i+1) = nB(i)-1;  %decreases count of present B molecules
    end
    t(i+1) = t(i)+dt(i);    %advances time
    
    matrix = [1-(k_f*dt(i)) k_r*dt(i); k_f*dt(i) 1-(k_r*dt(i))];    %matrix to describe the change in species population
    ODEpop = matrix*[ODEnA(i); ODEnB(i)]; %matrix multiplication to describe new populations
    ODEnA(i+1) = ODEpop(1,1);   %stores value for the next population amount
    ODEnB(i+1) = ODEpop(2,1);
end

figure();   %plots number of molecules over time
scatter(t,nA,5,'red','filled');
hold on;
scatter(t,nB,5,'black','filled');
plot(t,ODEnA,'red');
plot(t,ODEnB,'black');
xlabel('Time');
ylabel('Molecule Count');
title('A \leftrightarrow B');
legend('X_A','X_B');

eqA = round(mean(nA)); %average values of each type of molecule
eqB = round(mean(nB));
ODEeqA = round(mean(ODEnA));    %average values of ODE model
ODEeqB = round(mean(ODEnB));