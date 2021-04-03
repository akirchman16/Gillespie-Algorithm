clear all;
close all;

% This code is a simple example of the Gillespie Algorithm based on the
% wikipedia age on the Gillespie Algorithm (Reversible binding of A and B
% to form AB dimers). The reaction described here is A+B <--> AB.

t(1) = 0;   %starts time at zero

k_f = 5;    %rate constant of forward reaction (A+B-->AB)
k_r = 1;    %rate constant of reverse reaction (AB-->A+B)
nA(1) = 10; %initial number of A molecules
nB(1) = 0; %initial number of B molecultes
% nAB(1) = 0; %initial number of AB dimers
Iterations = 50;   %number of times the algorithm will run

for i = 1:Iterations
    R_tot(i) = (k_f.*nA(i))+(k_r.*nB(i));  %calculates the total reaction rate                   
    dt(i) = random('Exponential',1/R_tot(i)); %selects random value for time interval (stores history in array)
    
    P(i) = (k_f.*nA(i))./(R_tot(i));  %probability of A and B molecules binding event
    if rand <= P(i) %tests probability
        nA(i+1) = nA(i)-1;  %decreases count of present A molecules
        nB(i+1) = nB(i)+1;  %decreases count of present B molecules
%         nAB(i+1) = nAB(i)+1;    %increases count of present AB dimers
    else           %otherwise a dissociation event has to occur
        nA(i+1) = nA(i)+1;  %increases count of present A molecules
        nB(i+1) = nB(i)-1;  %increases count of present B molecules
%         nAB(i+1) = nAB(i)-1;    %decreases count of present AB dimers
    end
    t(i+1) = t(i)+dt(i);    %advances time
end

figure();   %plots number of molecules over time
plot(t,nA,'b');
hold on;
plot(t,nB,'black');
% plot(t,nAB,'r');
xlabel('Time');
ylabel('Molecule Count');
legend('nA','nB','nAB');

eqA = round(mean(nA)); %average values of each type of molecule
eqB = round(mean(nB));
% eqAB = round(mean(nAB));