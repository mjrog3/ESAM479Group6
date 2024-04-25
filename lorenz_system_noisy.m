function lorenz_system_noisy()
    sigma = 10;
    b = 8/3;
    r = 28; % Typical value for observing chaos
    initial_conditions = [1 1 1; 1.01 1 1]; % Initial conditions [x0; y0; z0]
    tspan = [0 10]; % Extended time span for the simulation
    noise_levels = [1 1 1]; % x y and z noise levels

    opts = odeset('RelTol',1e-5,'AbsTol',1e-8); % Setting solver options for better accuracy

    for i = 1:2
        [t, sol] = ode45(@(t, y) lorenzODE(t, y, sigma, b, r, noise_levels), tspan, initial_conditions(i,:), opts);
        x = sol(:, 1);
        y = sol(:, 2);
        z = sol(:, 3);

        figure; % Create a new figure window
        subplot(3, 1, 1);
        plot(t, x);
        title(['x(t) with r = ', num2str(r)]);
        xlabel('Time');
        ylabel('x');

        subplot(3, 1, 2);
        plot(t, y);
        title(['y(t) with r = ', num2str(r)]);
        xlabel('Time');
        ylabel('y');

        subplot(3, 1, 3);
        plot3(x, y, z);
        title('3D Trajectory of Lorenz Attractor');
        xlabel('x');
        ylabel('y');
        zlabel('z');
        grid on; % Enable grid for better visualization
    end
end

function dy = lorenzODE(t, y, sigma, b, r, noise_levels)
    dy = zeros(3,1);
    dy(1) = sigma * (y(2) - y(1)) + normrnd(0,noise_levels(1));
    dy(2) = y(1) * (r - y(3)) - y(2) + normrnd(0,noise_level(2));
    dy(3) = y(1) * y(2) - b * y(3) + normrnd(0,noise_level(3));
end
