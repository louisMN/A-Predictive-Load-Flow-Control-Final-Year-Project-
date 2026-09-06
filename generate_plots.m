% ============================================================
% generate_plots.m
% Produces Chapter 4 figures from the load-aware MPC simulation.
% Run AFTER mpc_gupfc.m (needs its variables in the workspace),
% or it will reload from mpc_results.mat and rebuild what it needs.
% ============================================================

% ---- Make sure we have the data ----
if ~exist('P_loss_mpc','var')
    load('mpc_results.mat');   % brings P_loss_mpc, V_dev_mpc, u_history, LF
end

% ---- Rebuild the fixed-GUPFC baseline curve if not in workspace ----
if ~exist('P_loss_fixed','var')
    % Rebuild contingency Ybus + GUPFC (same as mpc_gupfc.m)
    linedata_cont = linedata; linedata_cont(5,:) = [];
    nb = size(busdata,1);
    Ybus_cont = zeros(nb,nb);
    for k = 1:size(linedata_cont,1)
        f=linedata_cont(k,1); t=linedata_cont(k,2);
        R=linedata_cont(k,3); X=linedata_cont(k,4); B=linedata_cont(k,5);
        if R==0; y=1/(1j*X); else; y=1/(R+1j*X); end
        yc=1j*B/2;
        Ybus_cont(f,f)=Ybus_cont(f,f)+y+yc; Ybus_cont(t,t)=Ybus_cont(t,t)+y+yc;
        Ybus_cont(f,t)=Ybus_cont(f,t)-y;    Ybus_cont(t,f)=Ybus_cont(t,f)-y;
    end
    R_se=0.01; X_se=0.15; y_se=1/(R_se+1j*X_se);
    s_bus=3; r1_bus=7; r2_bus=17;
    Ybus_cont(s_bus,s_bus)=Ybus_cont(s_bus,s_bus)+y_se*2;
    Ybus_cont(r1_bus,r1_bus)=Ybus_cont(r1_bus,r1_bus)+y_se;
    Ybus_cont(r2_bus,r2_bus)=Ybus_cont(r2_bus,r2_bus)+y_se;
    Ybus_cont(s_bus,r1_bus)=Ybus_cont(s_bus,r1_bus)-y_se; Ybus_cont(r1_bus,s_bus)=Ybus_cont(r1_bus,s_bus)-y_se;
    Ybus_cont(s_bus,r2_bus)=Ybus_cont(s_bus,r2_bus)-y_se; Ybus_cont(r2_bus,s_bus)=Ybus_cont(r2_bus,s_bus)-y_se;

    P_nom = busdata(:,3); Q_nom = busdata(:,4);
    u_fixed = [0.05, 30, 0.05, 30, 0.10];
    P_loss_fixed = zeros(1,24);
    V_dev_fixed  = zeros(1,24);
    for t = 1:24
        [P_loss_fixed(t), V_dev_fixed(t)] = run_gupfc_lf(u_fixed, LF(t), ...
            P_nom, Q_nom, Ybus_cont, busdata, linedata_cont, ...
            s_bus, r1_bus, r2_bus, y_se);
    end
end

hours = 0:23;

% ============================================================
% FIGURE 1 — Active power losses: MPC vs fixed GUPFC over 24h
% ============================================================
f1 = figure('Color','w','Position',[100 100 800 450]);
plot(hours, P_loss_fixed, 's--', 'Color',[0.85 0.33 0.10], ...
     'LineWidth',1.8,'MarkerSize',6,'MarkerFaceColor',[0.85 0.33 0.10]); hold on;
plot(hours, P_loss_mpc, 'o-', 'Color',[0.00 0.45 0.74], ...
     'LineWidth',2.2,'MarkerSize',6,'MarkerFaceColor',[0.00 0.45 0.74]);
hold off; grid on;
xlabel('Hour of Day','FontSize',12);
ylabel('Active Power Loss (MW)','FontSize',12);
title('Active Power Losses over 24-Hour Load Cycle','FontSize',13,'FontWeight','bold');
lg1 = legend('Fixed GUPFC','Load-Aware MPC','Location','northwest');
set(lg1,'FontSize',11);
xlim([0 23]); set(gca,'FontSize',11,'XTick',0:2:23);
saveas(f1,'fig4_1_losses.png');

% ============================================================
% FIGURE 2 — Worst-case voltage deviation over 24h
% ============================================================
f2 = figure('Color','w','Position',[100 100 800 450]);
plot(hours, V_dev_mpc, 'o-', 'Color',[0.47 0.67 0.19], ...
     'LineWidth',2.2,'MarkerSize',6,'MarkerFaceColor',[0.47 0.67 0.19]); hold on;
plot([0 23],[0.060 0.060],'--','Color',[0.6 0.6 0.6],'LineWidth',1.5);
text(0.3, 0.063, 'Paper target (0.060)','Color',[0.6 0.6 0.6],'FontSize',10);
hold off; grid on;
xlabel('Hour of Day','FontSize',12);
ylabel('Worst-Case Voltage Deviation (p.u.)','FontSize',12);
title('Voltage Deviation under Load-Aware MPC','FontSize',13,'FontWeight','bold');
xlim([0 23]); ylim([0 0.08]); set(gca,'FontSize',11,'XTick',0:2:23);
saveas(f2,'fig4_2_voltage.png');

% ============================================================
% FIGURE 3 — Daily load factor profile (the novel input)
% ============================================================
f3 = figure('Color','w','Position',[100 100 800 450]);
area(hours, LF, 'FaceColor',[0.30 0.30 0.35],'FaceAlpha',0.25, ...
     'EdgeColor',[0.20 0.20 0.25],'LineWidth',2); hold on;
plot(hours, LF, 'o-','Color',[0.20 0.20 0.25],'LineWidth',2,'MarkerSize',5,...
     'MarkerFaceColor',[0.20 0.20 0.25]);
plot([0 23],[1.0 1.0],':','Color',[0.7 0.3 0.3],'LineWidth',1.3);
text(0.3, 1.02, 'Nominal (1.0)','Color',[0.7 0.3 0.3],'FontSize',10);
hold off; grid on;
xlabel('Hour of Day','FontSize',12);
ylabel('Load Factor (p.u.)','FontSize',12);
title('Time-Varying Daily Load Profile','FontSize',13,'FontWeight','bold');
xlim([0 23]); ylim([0.55 1.30]); set(gca,'FontSize',11,'XTick',0:2:23);
saveas(f3,'fig4_3_loadprofile.png');

% ============================================================
% FIGURE 4 — Combined 2x1: losses + voltage (compact for doc)
% ============================================================
f4 = figure('Color','w','Position',[100 100 800 700]);
subplot(2,1,1);
plot(hours, P_loss_fixed,'s--','Color',[0.85 0.33 0.10],'LineWidth',1.8,'MarkerSize',5); hold on;
plot(hours, P_loss_mpc,'o-','Color',[0.00 0.45 0.74],'LineWidth',2.2,'MarkerSize',5);
hold off; grid on;
ylabel('Loss (MW)','FontSize',11);
title('(a) Active Power Losses: Fixed vs Load-Aware MPC','FontSize',12,'FontWeight','bold');
lg4 = legend('Fixed GUPFC','Load-Aware MPC','Location','northwest');
set(lg4,'FontSize',10);
xlim([0 23]); set(gca,'FontSize',10,'XTick',0:2:23);

subplot(2,1,2);
plot(hours, V_dev_mpc,'o-','Color',[0.47 0.67 0.19],'LineWidth',2.2,'MarkerSize',5); hold on;
plot([0 23],[0.060 0.060],'--','Color',[0.6 0.6 0.6],'LineWidth',1.3);
hold off; grid on;
xlabel('Hour of Day','FontSize',11); ylabel('V dev (p.u.)','FontSize',11);
title('(b) Worst-Case Voltage Deviation under MPC','FontSize',12,'FontWeight','bold');
xlim([0 23]); ylim([0 0.08]); set(gca,'FontSize',10,'XTick',0:2:23);
saveas(f4,'fig4_4_combined.png');

% ---- Console summary of what was saved ----
fprintf('\n===== PLOTS GENERATED =====\n');
fprintf('fig4_1_losses.png       - MPC vs fixed losses\n');
fprintf('fig4_2_voltage.png      - voltage deviation\n');
fprintf('fig4_3_loadprofile.png  - daily load profile\n');
fprintf('fig4_4_combined.png     - combined losses+voltage\n');
fprintf('\nAvg MPC loss:   %.4f MW\n', mean(P_loss_mpc));
fprintf('Avg fixed loss: %.4f MW\n', mean(P_loss_fixed));
fprintf('Improvement:    %.2f %%\n', 100*(mean(P_loss_fixed)-mean(P_loss_mpc))/mean(P_loss_fixed));