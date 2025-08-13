function fixedBreakpointsDeltaApproximator(f_handle, domain, B, I)
    % FBSD with fixed breakpoints and visualization
    X_min = domain(1);
    X_max = domain(2);
    x = linspace(X_min, X_max, B)';  % fixed breakpoints

    % Interpolation grid for f(x)
    x_data = linspace(X_min, X_max, 1000);
    f_data = f_handle(x_data);

    % YALMIP variables
    s = sdpvar(B,1);
    phi = sdpvar(B,1);
    mu = sdpvar(1);

    Constraints = [];

    % phi_b = f(x_b) + s_b
    for b = 1:B
        f_interp = interp1(x_data, f_data, x(b), 'linear', 'extrap');
        Constraints = [Constraints, phi(b) == f_interp + s(b)];
    end

    % Shift bounds
    Constraints = [Constraints, -mu <= s <= mu];

    % Deviation constraints at test points
    for b = 2:B
        for i = 1:I
            t = i / (I + 1);
            xbi = (1 - t) * x(b-1) + t * x(b);
            fbi = interp1(x_data, f_data, xbi, 'linear', 'extrap');
            lbi = (1 - t) * phi(b-1) + t * phi(b);

            Constraints = [Constraints;
                lbi - fbi <= mu;
                fbi - lbi <= mu;
            ];
        end
    end

    Objective = mu;
    options = sdpsettings('solver','gurobi','verbose',0);
    t_start = tic;
    sol = optimize(Constraints, Objective, options);
    elapsed_time = toc(t_start);  % End timer
    if sol.problem == 0
        mu_val = value(mu);
        phi_val = value(phi);

        disp('FBSD solution found.');
        disp(['Minimized deviation μ = ', num2str(mu_val)]);
        disp(['Computation time: ', num2str(elapsed_time), ' seconds']);
        disp('Fixed breakpoints:');
        disp(x');

        % Plotting
        xx = linspace(X_min, X_max, 1000);
        fx = f_handle(xx);
        fx_upper = fx + mu_val;
        fx_lower = fx - mu_val;

        % Reconstruct piecewise linear approximation
        phi_interp = zeros(size(xx));
        for i = 2:B
            idx = xx >= x(i-1) & xx <= x(i);
            t = (xx(idx) - x(i-1)) / (x(i) - x(i-1));
            phi_interp(idx) = (1 - t) * phi_val(i-1) + t * phi_val(i);
        end

        % Plot
        figure;
        hold on;

        % Gray shaded region: f(x) ± μ
        fill([xx, fliplr(xx)], [fx_upper, fliplr(fx_lower)], ...
             [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);  % shaded band

        % Plot the actual and approximated functions
        plot(xx, fx, 'b-', 'LineWidth', 1);
        plot(xx, phi_interp, 'r-', 'LineWidth', 1.5);
        scatter(x, value(phi), 20, 'ro', 'filled');

        legend('f(x) \pm \delta', 'f(x)', 'Approximated f(x)', 'Location', 'Best');
        title(['FBSD Approximation with \delta = ', num2str(mu_val)]);
        xlabel('x'); ylabel('y');
        grid on; box on;
    else
        disp('Solver failed:');
        disp(sol.info);
    end
end

