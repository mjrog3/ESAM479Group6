function out = generate_noisy_lorenz(sigma, b, r, IC, tspan, noise_levels)
    % sigma = 10;
    % b = 8/3;
    % r = 28; % Typical value for observing chaos
    % initial_conditions = [1 1 1; 1.01 1 1]; % Initial conditions [x0; y0; z0]
    % tspan = [0 10]; % Extended time span for the simulation
    % noise_levels = [1 1 1]; % x y and z noise levels

    opts = odeset('RelTol',1e-5,'AbsTol',1e-8); % Setting solver options for better accuracy

    [t, sol] = ode45(@(t, y) lorenzODE(t, y, sigma, b, r, noise_levels), tspan, IC, opts);
    x = sol(:, 1);
    y = sol(:, 2);
    z = sol(:, 3);

    out = [t x y z];
end

function dy = lorenzODE(t, y, sigma, b, r, noise_levels)
    dy = zeros(3,1);
    dy(1) = sigma * (y(2) - y(1)) + normrnd(0,noise_levels(1));
    dy(2) = y(1) * (r - y(3)) - y(2) + normrnd(0,noise_levels(2));
    dy(3) = y(1) * y(2) - b * y(3) + normrnd(0,noise_levels(3));
end
