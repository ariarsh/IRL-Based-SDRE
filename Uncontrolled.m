% Uncontrolled Simulation of Nonlinear System
% System: dx/dt = A(x) x
% A(x) = [1-x1^2, 1; 1+x1*x2, -1]
% Q = I (2x2) for cost computation
% Computes total cost using trapezoidal rule and simulation time
% Plots x1, x2, and zero feedback gains K in subplots

clear; clc; close all;

% Parameters
Q = eye(2);  % State weighting (constant) for cost computation

% Start timing
tic;

% Simulation
x0 = [3; 1];         % Initial condition
tspan = [0, 10];     % Time span
[t, x] = ode45(@(t,x) uncontrolled_dynamics(t, x), tspan, x0);

% Compute cost (no control input, u = 0)
cost_integral = zeros(length(t), 1);
for i = 1:length(t)
    xi = x(i, :)';
    cost_integral(i) = xi' * Q * xi;
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
plot(t, zeros(size(t)), 'b-', 'LineWidth', 1.5, 'DisplayName', 'K_1'); hold on;
plot(t, zeros(size(t)), 'r--', 'LineWidth', 1.5, 'DisplayName', 'K_2');
ylabel('Feedback Gains');
xlabel('Time (s)');
legend('show');
grid on;


% Uncontrolled dynamics function for ode45
function dxdt = uncontrolled_dynamics(t, x)
    % Compute A(x)
    A = [1 - x(1)^2, 1; 1 + x(1)*x(2), -1];
    
    % Dynamics: dx/dt = A(x) x
    dxdt = A * x;
end