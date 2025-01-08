close all;
clear;
clc;

% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the 3 data points that the quadratic interpolation will pass through
xData = [0, 1, 2];
yData = [2, 3, 6];

% domain of the linear interpolation and lagrange functions
domain = linspace(min(xData), max(xData), 50);

% CALCULATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lagrange functions
x = domain;
x0 = xData(1);
x1 = xData(2);
x2 = xData(3);
L0 = (x - x1).*(x - x2)./((x0 - x1).*(x0 - x2));
L1 = (x - x2).*(x - x0)./((x1 - x2).*(x1 - x0));
L2 = (x - x1).*(x - x0)./((x2 - x1).*(x2 - x0));

% components of f
L0fx0 = L0*yData(1);
L1fx1 = L1*yData(2);
L2fx2 = L2*yData(3);

% f
f = L0fx0 + L1fx1 + L2fx2;

% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;
plot(xData, yData, 'ro', 'markerFaceColor', 'r');
plot(domain, L0);
plot(domain, L1);
plot(domain, L2);
plot(domain, f, 'r', 'lineWidth', 2);
xlabel('x');
ylabel('y');
legend('data points', 'L0(x)', 'L1(x)', 'L2(x)', 'f1(x)');

figure;
hold on;
plot(xData, yData, 'ro', 'markerFaceColor', 'r');
plot(domain, L0fx0);
plot(domain, L1fx1);
plot(domain, L2fx2);
plot(domain, f, 'r', 'lineWidth', 2);
xlabel('x');
ylabel('y');
legend('data points', 'L0(x)*f(x0)', 'L1(x)*f(x1)', 'L2(x)*f(x2)', 'f1(x)');