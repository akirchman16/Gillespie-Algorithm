clear all;
close all;

% This code will generate populations over time for the reaction: A+B<-->C
% It will compare the results from the Differential Equation model as well
% as the Direct and First Reaction Methods from the Gillespie Algorithm.

k_f = 0.1;    %kinetic rate constants for forward and reverse reactions
k_r = 0.1;

Iterations = 1000;   %number of individual reactions to occur
FinalTime = 5;  %final time for simulation

%Memory Allocation
t_DE = zeros(1,1+Iterations);
A_DE = zeros(1,1+Iterations);
B_DE = zeros(1,1+Iterations);
C_DE = zeros(1,1+Iterations);
diffA = zeros(1,Iterations);
diffB = zeros(1,Iterations);
diffC = zeros(1,Iterations);

%initial values
t_DE(1) = 0;    %initial time values for each method
t_Direct(1) = 0;
t_First(1) = 0;

A_DE(1) = 50;   %initial A-population values
A_Direct(1) = 50;
A_First(1) = 50;
B_DE(1) = 30;   %initial B-population values
B_Direct(1) = 30;
B_First(1) = 30;
C_DE(1) = 0;    %initial C-population values
C_Direct(1) = 0;
C_First(1) = 0;

% DIFFERENTIAL EQUATION MODEL
dt = FinalTime/Iterations;  %Time Step
for a = 2:Iterations+1
    diffA(a-1) = -k_f*A_DE(a-1)*B_DE(a-1)+k_r*C_DE(a-1); %DE for A
    diffB(a-1) = -k_f*A_DE(a-1)*B_DE(a-1)+k_r*C_DE(a-1); %DE for B
    diffC(a-1) = k_f*A_DE(a-1)*B_DE(a-1)-k_r*C_DE(a-1); %DE for C
    
    A_DE(a) = diffA(a-1)*dt+A_DE(a-1);  %Euler method for A
    B_DE(a) = diffB(a-1)*dt+B_DE(a-1);  %Euler method for B
    C_DE(a) = diffC(a-1)*dt+C_DE(a-1);  %Euler method for C
    
    t_DE(a) = t_DE(a-1)+dt;
end

% GILLESPIE - DIRECT METHOD
b = 1;  %counter for events in Direct Method
while t_Direct < FinalTime
    b = b+1;
    
    a_1_Direct(b-1) = k_f*A_Direct(b-1)*B_Direct(b-1); %propensity function for j=1
    a_2_Direct(b-1) = k_r*C_Direct(b-1);   %propensity function for j=2
    a_0(b-1) = a_1_Direct(b-1)+a_2_Direct(b-1);   %total propensity function
    
    r_1 = rand;
    r_2 = rand;
    tau_Direct(b-1) = (1/a_0(b-1))*log(1/rand);    %time interval
    
    if a_1_Direct(b-1) > r_2*a_0(b-1)  %forward reaction;
        A_Direct(b) = A_Direct(b-1)-1;
        B_Direct(b) = B_Direct(b-1)-1;
        C_Direct(b) = C_Direct(b-1)+1;
    else     %reverse reaction
        A_Direct(b) = A_Direct(b-1)+1;
        B_Direct(b) = B_Direct(b-1)+1;
        C_Direct(b) = C_Direct(b-1)-1;
    end
    
    t_Direct(b) = t_Direct(b-1) + tau_Direct(b-1); %advance time
end


% GILLESPIE - FIRST REACTION METHOD
c = 1;
while t_First < FinalTime
    c = c+1;
    
    a_1_First(c-1) = k_f*A_First(c-1)*B_First(c-1);
    a_2_First(c-1) = k_r*C_First(c-1);
    
    T = [(1/a_1_First(c-1))*log(1/rand),(1/a_2_First(c-1))*log(1/rand)];  %time intervals
    
    tau_First(c-1) = min(T);
    j_First(c-1) = find(T == min(T));
    
    if j_First(c-1) == 1    %forward reaction
        A_First(c) = A_First(c-1)-1;
        B_First(c) = B_First(c-1)-1;
        C_First(c) = C_First(c-1)+1;
    elseif j_First(c-1) == 2    %reverse reaction
        A_First(c) = A_First(c-1)+1;
        B_First(c) = B_First(c-1)+1;
        C_First(c) = C_First(c-1)-1;
    end
    
    t_First(c) = t_First(c-1)+tau_First(c-1);   %advance time
end

figure(1);  
subplot(2,2,1); %plot of A population
plot(t_DE,A_DE,'k');
hold on;
scatter(t_Direct,A_Direct,3,'r','filled');
scatter(t_First,A_First,3,'b','filled');
title('A Population: A + B <--> C');
ylim([0 A_DE(1)]);
xlim([0 FinalTime]);
xlabel('Time, t');
ylabel('Population');
legend('Diff. Eq.','Direct','First Reaction');
box on;

subplot(2,2,2); %plot of B population
plot(t_DE,B_DE,'k');
hold on;
scatter(t_Direct,B_Direct,3,'r','filled');
scatter(t_First,B_First,3,'b','filled');
title('B Population: A + B <--> C');
ylim([0 B_DE(1)]);
xlim([0 FinalTime]);
xlabel('Time, t');
ylabel('Population');
legend('Diff. Eq.','Direct','First Reaction');
box on;

subplot(2,2,3); %plot of C population
plot(t_DE,C_DE,'k');
hold on;
scatter(t_Direct,C_Direct,3,'r','filled');
scatter(t_First,C_First,3,'b','filled');
title('C Population: A + B <--> C');
ylim([0 max(C_DE)]);
xlim([0 FinalTime]);
xlabel('Time, t');
ylabel('Population');
legend('Diff. Eq.','Direct','First Reaction');
box on;

figure(2); %everything together
plot(t_DE,A_DE,'r');
hold on;
plot(t_DE,B_DE,'b');
plot(t_DE,C_DE,'g');
scatter(t_Direct,A_Direct,5,'r','+');
scatter(t_Direct,B_Direct,5,'b','+');
scatter(t_Direct,C_Direct,5,'g','+');
scatter(t_First,A_First,5,'r','^');
scatter(t_First,B_First,5,'b','^');
scatter(t_First,C_First,5,'g','^');
title('Everything Together');
xlim([0 FinalTime]);
ylim([0 max([A_DE(1),B_DE(1),max(C_DE)])]);
xlabel('Time, t');
ylabel('Populations');
box on;