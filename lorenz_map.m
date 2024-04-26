function lorenz_map()
    % Lorenz system parameters
    sigma = 10;
    b = 8/3;
    r = 28; % Typical value for observing chaos
    initial_conditions = [1; 1; 1]; % Initial conditions [x0; y0; z0]
    tspan = [0 1.5]; % Time span to observe one Lyapunov time

    % Solver options
    opts1 = odeset('RelTol',1e-7,'AbsTol',1e-7); % Higher precision
    opts2 = odeset('RelTol',1e-5,'AbsTol',1e-5); % Lower precision

    % Solve with high precision
    [t1, sol1] = ode45(@(t, y) lorenzODE(t, y, sigma, b, r), tspan, initial_conditions, opts1);
    x1 = sol1(:, 1);
    y1 = sol1(:, 2);
    z1 = sol1(:, 3);

    % Solve with lower precision
    [t2, sol2] = ode45(@(t, y) lorenzODE(t, y, sigma, b, r), tspan, initial_conditions, opts2);
    x2 = sol2(:, 1);
    y2 = sol2(:, 2);
    z2 = sol2(:, 3);

    % Plotting
    figure;
    subplot(2, 1, 1);
    plot(t1, x1, 'b', t2, x2, 'r--');
    legend('High precision', 'Low precision');
    title('Comparison of x(t) solutions');
    xlabel('Time');
    ylabel('x');

    subplot(2, 1, 2);
    plot3(x1, y1, z1, 'b', x2, y2, z2, 'r--');
    legend('High precision', 'Low precision');
    title('3D Trajectory of Lorenz Attractor');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    grid on;

end

function dy = lorenzODE(t, y, sigma, b, r)
    dy = zeros(3,1);
    dy(1) = sigma * (y(2) - y(1));
    dy(2) = y(1) * (r - y(3)) - y(2);
    dy(3) = y(1) * y(2) - b * y(3);
end
