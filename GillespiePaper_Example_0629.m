clear all;
close all;

% Gillespie algorithm simulation of a unimolecular chemical reaction
% (S-->P). Utilizes the Stochastic Simulation Algorithm (SSA) from the
% Gillespie paper (Daniel T. Gillespie, 2007). Copy of the example given in
% the Gillespie paper.

k_f = 1;            %reaction rate constant
Omega = 1;          %volume of system
Iterations = 100;   %number of times the code runs

a_j = zeros(1,Iterations);  %allocates memory
a_0 = zeros(1,Iterations);

c = k_f;      %constant of reaction is rate constant due to unimolecular
t(1) = 0;       %initializes time
x(1) = 100;     %initial value for reactant species (S1)
ODEx(1) = x(1); %initial value for RRE solution

for i = 1:Iterations
    a_j(i) = (c)*x(i);        %calculates propensity function based on reaction constant
    a_0(i) = a_j(i);      %calculates sum of all a_j(x)
    
    r_1(i) = rand; %random values from uniformly distribution over unit interval
    r_2(i) = rand;
    
    tau(i) = (1/a_0(i))*log(1/r_1(i)); %calculates tau value based on Gillespie Eq.10a
    for n = 1:Iterations        %for loop to calculate smallest integer for j (Gillespie Eq. 10b)
        if sum(a_j(1:n)) > (r_2(i))*a_0(i) %finds smallest integer satisfying Eq. 10b
            j(i) = n;      %stores j value
            break;      %if this works we have a j value!
        end
    end
    
    t(i+1) = t(i)+tau(i);               %advances time
    x(i+1) = x(i) - 1;                  %updates molecule populations
    
    if x(i+1) < 0   %ends loop if population of reactant hits zero
        break
    end
end

ODEx = 0:0.1:max(t);
ODEy = x(1)*exp(-c*ODEx);

figure();
scatter(t,x,5,'r','filled');    %plots SSA trajectories
hold on;    
plot(ODEx,ODEy,'b');               %plots RRE solution
xlabel('Time');
ylabel('X(t)');
ylim([0 x(1)]);
legend('Stochastic Simulation','ODE Simulation');
title('Unimolecular Reaction');