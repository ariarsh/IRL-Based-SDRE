% IRL-based SDRE Controller for Nonlinear System Regulation
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
    P = irl_are(A, B, Q, R);
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
    
    % Solve ARE using IRL
    P = irl_are(A, B, Q, R);
    
    % Feedback gain K = R^{-1} B' P
    K = (1/R) * B' * P;
    
    % Control input u = -K x
    u = -K * x;
    
    % Dynamics: dx/dt = A(x) x + B u
    dxdt = A * x + B * u;
end

% IRL function to solve ARE
function P = irl_are(A, B, Q, R)
    n = size(A, 1);
    nP = 100; % maximum iterations
    M = 15; % number of sampling trajectories (at least n(n+1)/2 = 3 for n=2)
    Ts = 0.001; % sampling time
    T = 0.1; % integral time
    time = 0:Ts:T;
    Nt = numel(time);
    
    K = zeros(nP + 1, n);
    poles = -(1:n); % initial poles for stabilization
    K(1, :) = place(A, B, poles); % initial policy
    
    for j = 1:nP
        PHI = zeros(M, n*(n+1)/2);
        SAI = zeros(M, 1);
        
        for i = 1:M
            x = zeros(n, Nt);
            x(:, 1) = randn(n, 1);
            u = zeros(1, Nt);
            r = zeros(1, Nt);
            
            for k = 1:Nt-1
                u(k) = -K(j, :) * x(:, k) + 0.01 * randn; % probe noise
                x(:, k+1) = x(:, k) + Ts * (A * x(:, k) + B * u(k));
                r(k) = x(:, k)' * Q * x(:, k) + u(k) * R * u(k);
            end
            
            SAI(i) = trapz(time(1:Nt-1), r(1:Nt-1)); % integral using trapezoidal rule
            PHI(i, :) = ComputeXbar(x(:, 1)') - ComputeXbar(x(:, Nt)');
        end
        
        Pbar = PHI \ SAI;
        P = ConvertPbarToP(Pbar, n);
        
        K(j+1, :) = inv(R) * B' * P;
        
        if norm(K(j+1, :) - K(j, :)) < 1e-4
            break;
        end
    end
    
    % Final P from the last iteration
    P = ConvertPbarToP(Pbar, n);
end

% Helper function: Compute Xbar (vectorized upper triangular part of x*x')
function Xbar = ComputeXbar(X)
    X = X(:)';
    Xbar = [];
    for i = 1:numel(X)
        Xbar = [Xbar X(i) * X(i:end)];
    end
end

% Helper function: Convert Pbar to symmetric P matrix
function P = ConvertPbarToP(Pbar, n)
    P = zeros(n, n);
    idx = 1;
    for i = 1:n
        for jj = i:n
            if i == jj
                P(i, i) = Pbar(idx);
            else
                P(i, jj) = Pbar(idx) / 2;
                P(jj, i) = P(i, jj);
            end
            idx = idx + 1;
        end
    end
end