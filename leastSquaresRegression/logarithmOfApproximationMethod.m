close all;
clear;
clc;

% raw data 
xRaw = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50];
yRaw = [16, 25, 32, 33, 38, 36, 39, 40, 42, 42];

% linearise data
xLinearised = 1./xRaw;
yLinearised = 1./yRaw;

% perform linear least-squares regression on the linearised data
[a, b] = leastSquaresLinear(xLinearised, yLinearised);

% solve for p and r 
r = 1./b;
p = a.*r;

x = linspace(0, 55, 100);
% plug p and r back into original equation
y = r.*(x./(p + x));

% plot
hold on;

plot(x, y, 'lineWidth', 1);
plot(xRaw, yRaw, 'o', 'markerFaceColor', 'r');

xlabel('x');
ylabel('y');