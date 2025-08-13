function alphaForwardHeuristc(f_handle, X_minus, X_plus, delta, alpha, D)
% α-Forward Heuristic with Backward Iteration + Plotting
% Uses δ as the only constraint (up to 1e-5 tolerance)
    t_start = tic;
    x_breaks = X_minus;
    s_breaks = 0;
    x_b = X_minus;

    while x_b < X_plus
        x_next = X_plus;
        found = false;

        while ~found
            x_next = x_b + alpha * (x_next - x_b);

            if x_next > X_plus
                x_next = X_plus;
            end

            for d = 1:D
                s_try = ((2*d / (D + 1)) - 1) * delta;

                x_grid = linspace(x_b, x_next, 50);
                f_vals = f_handle(x_grid);
                y0 = f_handle(x_b) + s_breaks(end);
                y1 = f_handle(x_next) + s_try;
                phi_vals = y0 + (y1 - y0) * (x_grid - x_b) / (x_next - x_b);

                % Only check deviation within δ + 1e-5
                if max(abs(phi_vals - f_vals)) <= delta + 1e-5
                    found = true;
                    s_next = s_try;
                    break;
                end
            end
        end

        if found
            x_breaks(end+1) = x_next;
            s_breaks(end+1) = s_next;
            x_b = x_next;
        else
            break;
        end
    end
    elapsed_time = toc(t_start);
    fprintf('α-Forward Heuristic completed in %.6f seconds\n', elapsed_time);
    % Plotting
    x_dense = linspace(X_minus, X_plus, 1000);
    f_dense = f_handle(x_dense);
    phi_dense = zeros(size(x_dense));

    for i = 2:length(x_breaks)
        idx = x_dense >= x_breaks(i-1) & x_dense <= x_breaks(i);
        y0 = f_handle(x_breaks(i-1)) + s_breaks(i-1);
        y1 = f_handle(x_breaks(i)) + s_breaks(i);
        phi_dense(idx) = y0 + (y1 - y0) * ...
            (x_dense(idx) - x_breaks(i-1)) / (x_breaks(i) - x_breaks(i-1));
    end

    figure;
    plot(x_dense, f_dense, 'b', 'LineWidth', 2); hold on;
    plot(x_dense, phi_dense, 'r--', 'LineWidth', 2);
    fill([x_dense fliplr(x_dense)], ...
         [f_dense + delta, fliplr(f_dense - delta)], ...
         [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    scatter(x_breaks, f_handle(x_breaks) + s_breaks, 'ro', 'filled');
    legend('f(x)', 'φ(x)', 'δ-tube', 'Breakpoints','Location','Best');
    title(sprintf('α-Forward Heuristic (δ = %.3f)', delta));
    xlabel('x'); ylabel('y');
    grid on;
end