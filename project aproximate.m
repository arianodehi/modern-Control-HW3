clc;
clear;
close all;
% parameters
R = 6371;
k = 1.2;
r = 0.5;
m = 800;
f = 10000;
g = 9.8;

% defining state equation
dx = @(t, x) [x(2); f/m - g*(R/(R + x(1))) - (k/m)*(x(2))^2 * exp(-x(1)/r)];

% sloving diffrentional equation
tspan = [0 100];
x0 = [0; 0];
[t, x] = ode45(dx, tspan, x0);

%ploting
figure;
subplot(2, 1, 1);
plot(t, x(:, 1), 'b', 'LineWidth', 2);
hold on;
plot(t, x(:, 2), 'r', 'LineWidth', 2);
xlabel('time');
ylabel('m&m/s');
legend('distance', 'velosity');
title('non linear ploting');

% aproximating
A = [0 1; -(g*R)/(R+x0(1))^2 -2*(k/m)*x0(2)*exp(-x0(1)/r)];
B = [0; f/m];
C = eye(2);
D = zeros(2, 1);
sys = ss(A, B, C, D);
[y_linear, t_linear] = lsim(sys, f*ones(size(t)), t, x0);


subplot(2, 1, 2);
plot(t, x(:, 1), 'b', 'LineWidth', 2);
hold on;
plot(t, y_linear(:, 1), 'g--', 'LineWidth', 2);
plot(t, x(:, 2), 'r', 'LineWidth', 2);
plot(t, y_linear(:, 2), 'm--', 'LineWidth', 2);
xlabel('time');
ylabel('m&m/s');
legend('distance ', 'aproximate distance', 'speed', 'aproximate speed');
title('linear ploting');