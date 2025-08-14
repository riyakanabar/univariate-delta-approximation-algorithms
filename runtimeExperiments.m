%% Runtime scaling experiments for FBSD & α-Forward on f(x) = sin(x)/x, x ∈ [-20,20]
clear; clc; close all;

%% Function and domain
f_handle = @(x) sin(x)./x;  
domain   = [-20, 20];

%% ---------------------- FBSD: time vs. number of breakpoints ----------------------
B_list = [2 3 4 5 6 7 8 9 10 11 12 14 16 17 18 20 22 24 25 27 28 30 32];
%B_list  = [4 5 6 7 8 9 10 11 12 14 15 18 20 21 25 28 30 35 40 50];  % adjust as needed
n_reps  = 3;                          % average over reps
time_fbsd = zeros(size(B_list));

% (Optional) warm-up call (un-timed) to avoid first-call overhead:
%fixedBreakpointsDeltaApproximator(f_handle, domain, B_list(1), 100); 

for k = 1:numel(B_list)
    B = B_list(k);
    t = zeros(n_reps,1);
    for r = 1:n_reps
        tic;
        fixedBreakpointsDeltaApproximator(f_handle, domain, B, 100); 
        t(r) = toc;
    end
    time_fbsd(k) = mean(t);
end

% Fit linear and quadratic
[p1_lin, ~]  = polyfit(B_list, time_fbsd, 1);
[p2_quad, ~] = polyfit(B_list, time_fbsd, 2);

% Compute R^2 for each
yhat_lin   = polyval(p1_lin, B_list);
R2_lin     = 1 - sum((time_fbsd - yhat_lin).^2) / sum((time_fbsd - mean(time_fbsd)).^2);
yhat_quad  = polyval(p2_quad, B_list);
R2_quad    = 1 - sum((time_fbsd - yhat_quad).^2) / sum((time_fbsd - mean(time_fbsd)).^2);
disp(["lin",R2_lin])
disp(["quad",R2_quad])
% plot the smooth curve from the winner
if R2_quad > R2_lin
    p_fit = p2_quad; fit_label_fbsd = sprintf('Quadratic fit (R^2=%.3f)', R2_quad);
else
    p_fit = p1_lin;  fit_label_fbsd = sprintf('Linear fit (R^2=%.3f)',   R2_lin);
end

% Create a fine x-grid for smooth fit curve
B_fine   = linspace(min(B_list), max(B_list), 300);
yhat_fine = polyval(p_fit, B_fine);

% Plot: points only for measured data; smooth red curve for the fit
figure('Name','FBSD runtime scaling','Color','w'); 
plot(B_list, time_fbsd, 'ko', 'MarkerFaceColor','k', 'MarkerSize',4); hold on; grid on;
plot(B_fine, yhat_fine, 'r-', 'LineWidth',2);
xlabel('Number of breakpoints (B)'); ylabel('Computation time (s)');
legend('Measured time', fit_label_fbsd, 'Location','northwest');


%% --------------- α-Forward: time vs. precision (δ) ----------------
delta_list = [2 1.8 1.65 1.5 1.45 1.35 1.2 1.1 1.0 0.9 0.8 0.7 0.75 0.6 0.5 0.3 0.20 0.15 0.12 0.10 0.08 0.01]; 
alpha = 0.99; 
D = 100; 
time_af = zeros(size(delta_list)); 
n_reps = 3; 
% Warm-up (optional) 
%alphaForwardHeuristic(f_handle, domain(1), domain(2), delta_list(1), alpha, D); 
for k = 1:numel(delta_list) 
    del = delta_list(k); 
    t = zeros(n_reps,1); 
    for r = 1:n_reps 
        tic; 
        alphaForwardHeuristic(f_handle, domain(1), domain(2), del, alpha, D); 
        t(r) = toc; 
    end 
    time_af(k) = mean(t); 
end
% -------- Robust piecewise fit (auto split with spike detection) ----------
% Sort δ descending so high-δ (flat) comes first
[delta_s, idx_s] = sort(delta_list, 'descend');
time_s = time_af(idx_s);

n = numel(delta_s);

% 1) Estimate baseline from the top-N high-δ points
Nbase = min(max(5, round(0.25*n)), n-4);   
mu0   = mean(time_s(1:Nbase));
sd0   = std(time_s(1:Nbase));
kthr  = 3;                                

% 2) First spike
idx_spike = find(time_s(Nbase+1:end) > mu0 + kthr*sd0, 1, 'first');
if ~isempty(idx_spike)
    idx_spike = idx_spike + Nbase;          % convert to absolute index
    search_max = max(idx_spike-1, 2);       % constant region cannot pass the spike
else
    search_max = n-3;
end

best.SSE = inf; best.idxSplit = NaN; best.model = ''; best.p = []; best.R2 = NaN;

for s = 2 : min(search_max, n-3)   % at least 2 high-δ points, 3 low-δ points
    % High-δ constant segment
    y_hi   = time_s(1:s);
    cval   = mean(y_hi);
    SSE_hi = sum( (y_hi - cval).^2 );

    % Low-δ candidate fits
    x_lo   = delta_s(s+1:end);
    y_lo   = time_s(s+1:end);
    SST_lo = sum( (y_lo - mean(y_lo)).^2 );

    % Linear tail
    p_lin = polyfit(x_lo, y_lo, 1);
    yhat_lin = polyval(p_lin, x_lo);
    SSE_lin = sum((y_lo - yhat_lin).^2);

    % Quadratic tail
    p_quad = polyfit(x_lo, y_lo, 2);
    yhat_quad = polyval(p_quad, x_lo);
    SSE_quad = sum((y_lo - yhat_quad).^2);

    if SSE_lin <= SSE_quad
        SSE_lo = SSE_lin;  p_lo = p_lin;  model_lo = 'linear';
        R2_lo  = 1 - SSE_lin / max(SST_lo, eps);
    else
        SSE_lo = SSE_quad; p_lo = p_quad; model_lo = 'quadratic';
        R2_lo  = 1 - SSE_quad / max(SST_lo, eps);
    end

    SSE_tot = SSE_hi + SSE_lo;
    if SSE_tot < best.SSE
        best.SSE      = SSE_tot;
        best.idxSplit = s;
        best.cval     = cval;
        best.model    = model_lo;
        best.p        = p_lo;
        best.R2       = R2_lo;
    end
end

s = best.idxSplit;
delta_star = 0.5*(delta_s(s) + delta_s(s+1));  

% High-δ constant segment
delta_hi = linspace(delta_s(1), delta_s(s), 200);
y_hi     = best.cval * ones(size(delta_hi));

% Low-δ fitted segment
delta_lo = linspace(delta_s(s+1), delta_s(end), 200);
y_lo     = polyval(best.p, delta_lo);

% Plot (markers only + smooth piecewise curve)
figure('Name','Alpha-forward runtime scaling (piecewise)','Color','w');
plot(delta_list, time_af, 'ks', 'MarkerFaceColor','k', 'MarkerSize',4); hold on; grid on;
plot(delta_hi, y_hi, 'r-', 'LineWidth', 2);
plot(delta_lo, y_lo, 'r-', 'LineWidth', 2);
set(gca,'XDir','reverse');  xlabel('Precision \delta'); ylabel('Computation time (s)');
legend('Measured time', ...
       sprintf('Piecewise: const (high \\delta) + %s (low \\delta, R^2=%.3f)', best.model, best.R2), ...
       'Location','northwest');
fprintf('Estimated split δ* ≈ %.4f | Low-δ model: %s | R^2 = %.3f\n', delta_star, best.model, best.R2);
