close all;
clear;
clc;

% actual function
x = linspace(0, 5, 50);
fx = 1 - exp(-x);

% input
x0 = 1;
x1 = 2;
fx0 = 1 - exp(-x0);
fx1 = 1 - exp(-x1);

% calculate the linear newton interpolation f1(x)
f1x = fx0 + (fx1 - fx0)/(x1 - x0) * (x - x0);

% plot
hold on;

plot(x, fx, 'lineWidth', 1);
plot(x, f1x, 'lineWidth', 1);

plot(x0, fx0, 'o', 'markerFaceColor', 'r');
plot(x1, fx1, 'o', 'markerFaceColor', 'r');

xlabel('x');