clear all;
close all;

%Attempts to create an exponential distribution using a Monte Carlo
%technique as described in the "Stochastic Processes" powerpoint in the
%google drive folder.

m = 5;
Iterations = 1000;

P_max = m;  %maximum value of PDF

for i = 1:Iterations
    test = 0;   %determines whether x_try value will be kept or rejected
    while test == 0
        x_try(i) = rand;   %random trial value between 0 and 1
        PDF_x(i) = m*exp(-m*x_try(i)); %PDF of random value x_try
        if PDF_x(i)/P_max >= rand
            Dist(i) = x_try(i); %accept x_try
            test = 1;
        else
            test = 0;   %otherwise reject x_try
        end
    end
end
Average_MC = mean(Dist);

figure();
histogram(Dist,100);
hold on;
title('Monte Carlo Exponential');