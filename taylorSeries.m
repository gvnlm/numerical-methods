close all
clear
clc

syms x;

% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nth order taylor polynomial for f(x) about x0
f(x) = 5*x^3 - 3*x + 4/(2-x);
x0 = 0.5;
n = 5;

% domain of plot
domain = linspace(-1, 1, 100);   

% CALCULATE Pn(x) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% P0(x)
df(x) = f(x);
Pn(x) = df(x0)/factorial(0) * (x - x0)^0; % = f(x0)

% P1(x) - Pn(x)
for k = 1 : n

    df(x) = diff(df(x));
    Pn(x) = Pn(x) + df(x0)/factorial(k) * (x - x0)^k;
    
end

% print Pn(x)
fprintf('Pn(x) = %s\n', char(Pn(x)));

% PLOT f(x) AND Pn(x) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on;

plot(domain, f(domain), 'lineWidth', 1);
plot(domain, Pn(domain), 'o', 'markerSize', 3.5, 'markerFaceColor', 'r');
xline(x0, '--');

xlabel('x');
legend('f(x)', 'Pn(x)', 'x0');