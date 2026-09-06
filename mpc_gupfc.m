
% Step 7: Load-Aware MPC for GUPFC
% Optimises r1, gamma1, r2, gamma2, Qsh at each time step
% using predicted load trajectory

% Load factor profile
LF = [0.67, 0.65, 0.65, 0.66, 0.68, 0.72, 0.80, 0.88, ...
      0.92, 0.95, 0.97, 0.99, 1.00, 0.99, 0.98, 0.97, ...
      1.00, 1.05, 1.10, 1.15, 1.18, 1.20, 1.15, 1.05];

% GUPFC placement
s_bus  = 3;
r1_bus = 7;
r2_bus = 17;

% Series impedance
R_se = 0.01; X_se = 0.15;
Z_se = R_se + 1j*X_se;
y_se = 1/Z_se;

% MPC weights (cost function)
w_v = 10.0;   % voltage deviation weight
w_p = 1.0;    % power loss weight
w_u = 0.1;    % control effort weight

% Prediction horizon
Np = 5;

% GUPFC variable bounds
lb = [0,    0,   0,    0,   -0.5];
ub = [0.15, 360, 0.15, 360,  0.5];

% Build contingency Y-bus (branch 5 removed)
linedata_cont = linedata;
linedata_cont(5,:) = [];
nb = size(busdata,1);
Ybus_cont = zeros(nb,nb);
for k = 1:size(linedata_cont,1)
    from=linedata_cont(k,1); to=linedata_cont(k,2);
    R=linedata_cont(k,3); X=linedata_cont(k,4); B=linedata_cont(k,5);
    if R==0; y=1/(1j*X); else; y=1/(R+1j*X); end
    yc=1j*B/2;
    Ybus_cont(from,from)=Ybus_cont(from,from)+y+yc;
    Ybus_cont(to,to)=Ybus_cont(to,to)+y+yc;
    Ybus_cont(from,to)=Ybus_cont(from,to)-y;
    Ybus_cont(to,from)=Ybus_cont(to,from)-y;
end

% Add GUPFC virtual branches
Ybus_cont(s_bus,s_bus)   = Ybus_cont(s_bus,s_bus)   + y_se*2;
Ybus_cont(r1_bus,r1_bus) = Ybus_cont(r1_bus,r1_bus) + y_se;
Ybus_cont(r2_bus,r2_bus) = Ybus_cont(r2_bus,r2_bus) + y_se;
Ybus_cont(s_bus,r1_bus)  = Ybus_cont(s_bus,r1_bus)  - y_se;
Ybus_cont(r1_bus,s_bus)  = Ybus_cont(r1_bus,s_bus)  - y_se;
Ybus_cont(s_bus,r2_bus)  = Ybus_cont(s_bus,r2_bus)  - y_se;
Ybus_cont(r2_bus,s_bus)  = Ybus_cont(r2_bus,s_bus)  - y_se;

% Storage arrays
P_loss_mpc = zeros(1,24);
V_dev_mpc  = zeros(1,24);
u_history  = zeros(24,5);

% Initial GUPFC settings
u_prev = [0.05, 30, 0.05, 30, 0.10];

% Nominal loads
P_nom = busdata(:,3);
Q_nom = busdata(:,4);

fprintf('Running Load-Aware MPC over 24 hours...\n\n');

for t = 1:24
    % Get predicted load trajectory for next Np steps
    t_pred = t:min(t+Np-1, 24);
    if length(t_pred) < Np
        t_pred = [t_pred, repmat(24, 1, Np-length(t_pred))];
    end
    LF_pred = LF(t_pred);

    % Optimise GUPFC settings using fmincon
    options = optimoptions('fmincon', 'Display', 'off', ...
                           'MaxIterations', 100, 'TolFun', 1e-6);

    u_opt = fmincon(@(u) mpc_cost_fn(u, LF_pred, P_nom, Q_nom, ...
                    Ybus_cont, busdata, linedata_cont, s_bus, r1_bus, ...
                    r2_bus, y_se, w_v, w_p, w_u, u_prev, Np), ...
                    u_prev, [], [], [], [], lb, ub, [], options);

    % Apply first optimal action
    [Ploss_t, Vdev_t] = run_gupfc_lf(u_opt, LF(t), P_nom, Q_nom, ...
                                      Ybus_cont, busdata, linedata_cont, ...
                                      s_bus, r1_bus, r2_bus, y_se);

    P_loss_mpc(t)  = Ploss_t;
    V_dev_mpc(t)   = Vdev_t;
    u_history(t,:) = u_opt;
    u_prev = u_opt;

    fprintf('Hour %2d | LF=%.2f | Loss=%.4f MW | Vdev=%.4f | r1=%.3f g1=%.1f r2=%.3f g2=%.1f Qsh=%.3f\n', ...
        t-1, LF(t), Ploss_t, Vdev_t, u_opt(1), u_opt(2), u_opt(3), u_opt(4), u_opt(5));
end

% Summary
% ===== Fair comparison: fixed GUPFC evaluated at the SAME load profile =====
u_fixed = [0.05, 30, 0.05, 30, 0.10];   % the fixed setting used as baseline
P_loss_fixed = zeros(1,24);
for t = 1:24
    P_loss_fixed(t) = run_gupfc_lf(u_fixed, LF(t), P_nom, Q_nom, ...
                        Ybus_cont, busdata, linedata_cont, ...
                        s_bus, r1_bus, r2_bus, y_se);
end

% Nominal-load index (LF = 1.00) for clean head-to-head
[~, idx_nom] = min(abs(LF - 1.00));

fprintf('\n===== MPC RESULTS SUMMARY =====\n');
fprintf('Average losses, fixed GUPFC:   %.4f MW\n', mean(P_loss_fixed));
fprintf('Average losses, load-aware MPC: %.4f MW\n', mean(P_loss_mpc));
fprintf('Average improvement:            %.2f %%\n', ...
        100*(mean(P_loss_fixed)-mean(P_loss_mpc))/mean(P_loss_fixed));

fprintf('\n--- At nominal load (LF=1.00, hour %d) ---\n', idx_nom-1);
fprintf('Fixed GUPFC:      %.4f MW\n', P_loss_fixed(idx_nom));
fprintf('Load-aware MPC:   %.4f MW\n', P_loss_mpc(idx_nom));
fprintf('Improvement:      %.2f %%\n', ...
        100*(P_loss_fixed(idx_nom)-P_loss_mpc(idx_nom))/P_loss_fixed(idx_nom));

fprintf('\n--- Voltage performance (worst-case load bus) ---\n');
fprintf('MPC peak Vdev over 24h:  %.4f p.u.\n', max(V_dev_mpc));
fprintf('MPC mean Vdev over 24h:  %.4f p.u.\n', mean(V_dev_mpc));

fprintf('\n--- Reference (single-load snapshots) ---\n');
fprintf('Contingency, no GUPFC:   32.95 MW\n');
fprintf('Peak demand (LF=1.20):   MPC = %.4f MW, Fixed = %.4f MW\n', ...
        P_loss_mpc(22), P_loss_fixed(22));

% Save results
save('mpc_results.mat', 'P_loss_mpc', 'V_dev_mpc', 'u_history', 'LF');

% Plot
hours = 0:23;
figure(2);
subplot(2,1,1);
plot(hours, P_loss_mpc, 'b-o', 'LineWidth', 2);
xlabel('Hour'); ylabel('Losses (MW)');
title('Active Losses with Load-Aware MPC GUPFC');
grid on;

subplot(2,1,2);
plot(hours, V_dev_mpc, 'r-o', 'LineWidth', 2);
xlabel('Hour'); ylabel('Voltage Deviation (p.u.)');
title('Voltage Deviation with Load-Aware MPC');
grid on;
