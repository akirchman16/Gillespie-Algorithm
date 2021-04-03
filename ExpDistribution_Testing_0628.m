clear all;

%testing how to choose a value from a custom exponential distribution.

for i = 1:1000
    r(i) = random('Exponential',3);
end

figure();
histogram(r,100);