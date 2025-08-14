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

tol = 0.02;  

[delta_s, idx_s] = sort(delta_list, 'descend');
time_s = time_af(idx_s);

% Start from largest δ, find where values deviate > tol from initial mean
candidates = [];
for s = 2:numel(delta_s)-2  % leave at least 2 points for the tail
    mean_hi = mean(time_s(1:s));
    if abs(time_s(s+1) - mean_hi) > tol
        break;  % deviation found
    end
    candidates = [candidates, s];
end

if isempty(candidates)
    split_idx = round(numel(delta_s)/2); % fallback
else
    split_idx = candidates(end);  % last index in flat region
end

% Split δ value for reporting
delta_star = delta_s(split_idx);

% --- High-δ constant segment ---
delta_hi = linspace(delta_s(1), delta_star, 200);
cval = mean(time_s(1:split_idx));
y_hi = cval * ones(size(delta_hi));

% --- Low-δ linear fit ---
x_lo = delta_s(split_idx+1:end);
y_lo = time_s(split_idx+1:end);

[x_lo_asc, ord] = sort(x_lo, 'ascend');
y_lo_asc = y_lo(ord);

p_lin = polyfit(x_lo_asc, y_lo_asc, 1);
delta_lo = linspace(min(x_lo_asc), max(x_lo_asc), 200);
y_lo_fit = polyval(p_lin, delta_lo);

% R² for tail
SSE = sum((y_lo_asc - polyval(p_lin, x_lo_asc)).^2);
SST = sum((y_lo_asc - mean(y_lo_asc)).^2);
R2_tail = 1 - SSE / max(SST, eps);

% --- Plot ---
figure('Name','Auto piecewise const + linear','Color','w');
plot(delta_list, time_af, 'ks', 'MarkerFaceColor','k', 'MarkerSize',4); hold on; grid on;
plot(delta_hi, y_hi, 'r-', 'LineWidth', 2);
plot(delta_lo, y_lo_fit, 'r-', 'LineWidth', 2);
set(gca,'XDir','reverse');
xlabel('Precision \delta'); ylabel('Computation time (s)');
legend('Measured time', ...
    sprintf('Piecewise Constant Linear'), ...
    'Location','northwest');

fprintf('Auto δ* ≈ %.4f | Const=%.4f | Tail slope=%.4f | R^2=%.3f\n', ...
    delta_star, cval, p_lin(1), R2_tail);
