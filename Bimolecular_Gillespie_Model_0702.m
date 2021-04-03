clear all;
close all;

% This code will model a bimolecular reaction over time using the Gillespie
% algorithm. The reaction being modeled is the reversible reaction of
% A+B<-->C. It will plot populations of molecules over time and give a
% total time of the reaction.

xA(1) = 1000;   %initial population of A
xB(1) = 1000;   %initial population of B
xC(1) = 0;     %initial population of C
Omega = 1;     %volume of container in which reaction occurs
k_f = 0.1;    %reaction rate constant for the forward reaction (A+B-->C)
k_r = 0.1;    %reaction rate constant for the reverse reaction (C-->A+B)

t(1) = 0;           %initializes time at t=0
c_f = k_f/Omega;    %calculation of c for forward reaction
c_r = k_r;          %calculation of c for reverse reaction

ODExA(1) = xA(1);   %intial populations for ODE approximation
ODExB(1) = xB(1);
ODExC(1) = xC(1);

Equilibrium = 0;
Loops = 0;

ForwardCount = 0;
ReverseCount = 0;

while ~Equilibrium
    a_f(Loops+1) = c_f*xA(Loops+1)*xB(Loops+1);   %propensity functions
    a_r(Loops+1) = c_r*xC(Loops+1);
    a_0(Loops+1) = a_f(Loops+1)+a_r(Loops+1);     %sum of propensity functions
    
    r_1(Loops+1) = rand; %random numbers to be used to calculate tau and j
    r_2(Loops+1) = rand;
    
    tau(Loops+1) = (1/a_0(Loops+1))*log(1/r_1(Loops+1));  %Monte Carlo method to find tau
    if a_f(Loops+1) > r_2(Loops+1)*a_0(Loops+1)  %determines which reaction occurs
        j(Loops+1) = 1;           %forward reaction
    else
        j(Loops+1) = 2;           %reverse reaction
    end
    
    if j(Loops+1) == 1               %updates populations for a forward reaction
        xA((Loops+1)+1) = xA(Loops+1)-1;
        xB((Loops+1)+1) = xB(Loops+1)-1;
        xC((Loops+1)+1) = xC(Loops+1)+1;
        ForwardCount = ForwardCount+1;
    elseif j(Loops+1) == 2           %updates populations for a reverse reaction
        xA((Loops+1)+1) = xA(Loops+1)+1;
        xB((Loops+1)+1) = xB(Loops+1)+1;
        xC((Loops+1)+1) = xC(Loops+1)-1;
        ReverseCount = ReverseCount+1;
    end
    t((Loops+1)+1) = t(Loops+1)+tau(Loops+1);   %advances time by time step, tau
    
    state(Loops+1,1) = xA(Loops+1);   %creates history of state vectors after each
    state(Loops+1,2) = xB(Loops+1);   %reaction
    state(Loops+1,3) = xC(Loops+1);
    
    dx = [(-c_f*tau(Loops+1)),(c_r*tau(Loops+1)); (-c_f*tau(Loops+1)),(c_r*tau(Loops+1)); (c_f*tau(Loops+1)),(-c_r*tau(Loops+1))]*[ODExA(Loops+1)*ODExB(Loops+1);ODExC(Loops+1)];
    ODEold = [ODExA(Loops+1); ODExB(Loops+1); ODExC(Loops+1)];    %Euler method to solve ODEs
    ODEnew = dx + ODEold;
    ODExA((Loops+1)+1) = ODEnew(1,1);
    ODExB((Loops+1)+1) = ODEnew(2,1);
    ODExC((Loops+1)+1) = ODEnew(3,1);
    
    ODEstate(Loops+1,1) = ODExA(Loops+1);
    ODEstate(Loops+1,2) = ODExB(Loops+1);
    ODEstate(Loops+1,3) = ODExC(Loops+1);
    
    if Loops > 1000
        xC_EqTest = [xC(Loops-100) xC(Loops)]; %state of system 50 loops ago and now
        xC_Change = diff(xC_EqTest);    %difference between the two states
        if abs((xC_Change)) <= 1    %if system hasn't changed by more than 1 molecule in 50 iterations it's at equilibrium
            Equilibrium = 1;
        else
            Equilibrium = 0;
        end
    end
    
    Loops = Loops+1;
end

eqxA = round(mean(xA(round(0.75*Loops):Loops)));  %equilibrium pop.
eqxB = round(mean(xB(round(0.75*Loops):Loops)));
eqxC = round(mean(xC(round(0.75*Loops):Loops)));
eqODExA = round(mean(ODExA(round(0.75*Loops):Loops)));   %ODE equilibrium
eqODExB = round(mean(ODExB(round(0.75*Loops):Loops)));
eqODExC = round(mean(ODExC(round(0.75*Loops):Loops)));

figure();                          %plot of molecular populations
scatter(t,xA,5,'red','filled');
hold on;
scatter(t,xB,5,'green','filled');
scatter(t,xC,5,'blue','filled');
plot(t,ODExA,'red');
plot(t,ODExB,'green');
plot(t,ODExC,'blue');
xlabel('Time');
ylabel('Populations');
ylim([0 max(max(state))]);
title('A+B<-->C');
legend('xA','xB','xC');

total_time = max(t) %displays total time of reaction after so many reactions

% figure();
% histogram(tau,100);   %histogram of time steps
% xlabel('Tau, Time Step');
% title('Tau Histogram');