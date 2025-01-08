close all;
clear;
clc;

% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the 2 data points that the linear interpolation will pass through
xData = [0, 1];
yData = [1, 2];

% domain of the linear interpolation and lagrange functions
domain = linspace(min(xData), max(xData), 50);

% CALCULATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lagrange functions
L0 = (domain - xData(2))./(xData(1) - xData(2));
L1 = (domain - xData(1))./(xData(2) - xData(1));

% components of f
L0fx0 = L0*yData(1);
L1fx1 = L1*yData(2);

% f
f = L0fx0 + L1fx1;

% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;
plot(xData, yData, 'ro', 'markerFaceColor', 'r');
plot(domain, L0);
plot(domain, L1);
plot(domain, f, 'r', 'lineWidth', 2);
xlabel('x');
ylabel('y');
legend('data points', 'L0(x)', 'L1(x)', 'f1(x)');

figure;
hold on;
plot(xData, yData, 'ro', 'markerFaceColor', 'r');
plot(domain, L0fx0);
plot(domain, L1fx1);
plot(domain, f, 'r', 'lineWidth', 2);
xlabel('x');
ylabel('y');
legend('data points', 'L0(x)*f(x0)', 'L1(x)*f(x1)', 'f1(x)');