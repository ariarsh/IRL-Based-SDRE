% SDRE Controller for Nonlinear System Regulation
% System: dx/dt = A(x) x + B u
% A(x) = [1-x1^2, 1; 1+x1*x2, -1], B = [1; 1]
% Q = I (2x2), R = 1 (constants)
% Computes total cost using trapezoidal rule and simulation time
% Plots x1, x2, and feedback gains K in subplots

clear; clc; close all;

% Parameters
Q = eye(2);  % State weighting (constant)
R = 1;       % Input weighting (constant)

% Start timing
tic;

% Simulation
x0 = [3; 1];         % Initial condition
tspan = [0, 10];     % Time span
[t, x] = ode45(@(t,x) sdre_dynamics(t, x, Q, R), tspan, x0);

% Compute control input history, cost, and feedback gain
u = zeros(length(t), 1);
K_history = zeros(length(t), 2); % Store K = [K_1, K_2]
cost_integral = zeros(length(t), 1);
for i = 1:length(t)
    xi = x(i, :)';
    A = [1 - xi(1)^2, 1; 1 + xi(1)*xi(2), -1];
    B = [1; 1];
    [P, ~, ~] = care(A, B, Q, R);
    K = (1/R) * B' * P;
    K_history(i, :) = K;
    u(i) = -K * xi;
    cost_integral(i) = xi' * Q * xi + u(i) * R * u(i);
end

% Compute total cost using trapezoidal rule
total_cost = trapz(t, cost_integral);

% Stop timing
sim_time = toc;

% Display results
fprintf('Total cost: %.4f\n', total_cost);
fprintf('Simulation time: %.4f seconds\n', sim_time);

% Plots
figure;
subplot(2,1,1);
plot(t, x(:,1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'x_1'); hold on;
plot(t, x(:,2), 'r--', 'LineWidth', 1.5, 'DisplayName', 'x_2');
ylabel('States');
legend('show');
grid on;

subplot(2,1,2);
plot(t, K_history(:,1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'K_1'); hold on;
plot(t, K_history(:,2), 'r--', 'LineWidth', 1.5, 'DisplayName', 'K_2');
ylabel('Feedback Gains');
xlabel('Time (s)');
legend('show');
grid on;


% SDRE-controlled dynamics function for ode45
function dxdt = sdre_dynamics(t, x, Q, R)
    % Compute A(x)
    A = [1 - x(1)^2, 1; 1 + x(1)*x(2), -1];
    
    % B (constant)
    B = [1; 1];
    
    % Solve ARE: A'P + PA - P B R^{-1} B' P + Q = 0
    [P, ~, ~] = care(A, B, Q, R);
    
    % Feedback gain K = R^{-1} B' P
    K = (1/R) * B' * P;
    
    % Control input u = -K x
    u = -K * x;
    
    % Dynamics: dx/dt = A(x) x + B u
    dxdt = A * x + B * u;
end