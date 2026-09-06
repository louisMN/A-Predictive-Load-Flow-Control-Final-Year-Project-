% ============================================================
% frt_test.m
% Tests FRT voltage recovery at different load levels.
% Iterates the load flow to find the final steady-state FRT voltage.
% Requires: load_data.m, build_Ybus.m
% ============================================================

V_FRT_THRESH = 0.90;   % recovery target
V_REF        = 1.00;   % reference voltage
K_FRT        = 3.0;    % reactive-current gain (Eq. 2.8); grid-code 2-6
FAULT_BUS    = 30;     % weakest load bus
Y_FAULT      = 0.45;   % fault shunt admittance to ground (p.u.) -> sag depth
Q_CAP        = 1.5;    % converter reactive rating (saturation limit, p.u.)

s_bus=3; r1_bus=7; r2_bus=17;
R_se=0.01; X_se=0.15; y_se=1/(R_se+1j*X_se);

% Contingency Ybus (branch 5 out) + GUPFC
linedata_cont = linedata; linedata_cont(5,:) = [];
nb = size(busdata,1);
Ybus_base = zeros(nb,nb);
for k = 1:size(linedata_cont,1)
    f=linedata_cont(k,1); t=linedata_cont(k,2);
    R=linedata_cont(k,3); X=linedata_cont(k,4); B=linedata_cont(k,5);
    if R==0; y=1/(1j*X); else; y=1/(R+1j*X); end
    yc=1j*B/2;
    Ybus_base(f,f)=Ybus_base(f,f)+y+yc; Ybus_base(t,t)=Ybus_base(t,t)+y+yc;
    Ybus_base(f,t)=Ybus_base(f,t)-y;    Ybus_base(t,f)=Ybus_base(t,f)-y;
end
Ybus_base(s_bus,s_bus)=Ybus_base(s_bus,s_bus)+y_se*2;
Ybus_base(r1_bus,r1_bus)=Ybus_base(r1_bus,r1_bus)+y_se;
Ybus_base(r2_bus,r2_bus)=Ybus_base(r2_bus,r2_bus)+y_se;
Ybus_base(s_bus,r1_bus)=Ybus_base(s_bus,r1_bus)-y_se; Ybus_base(r1_bus,s_bus)=Ybus_base(r1_bus,s_bus)-y_se;
Ybus_base(s_bus,r2_bus)=Ybus_base(s_bus,r2_bus)-y_se; Ybus_base(r2_bus,s_bus)=Ybus_base(r2_bus,s_bus)-y_se;

u_base = [0.15, 65, 0.0, 30, 0.0];

LF_levels = [0.65, 1.00, 1.20, 1.10];
LF_names  = {'Off-peak','Nominal','Peak','Transition'};

fprintf('\n=============== FRT TEST UNDER VARYING LOAD ===============\n');
fprintf('Faulted bus: %d | threshold V < %.2f | K_FRT = %.1f | Q cap %.1f\n', ...
        FAULT_BUS, V_FRT_THRESH, K_FRT, Q_CAP);
fprintf('Recovery voltage is SOLVED by load flow (no tuned gain).\n\n');
fprintf('%-11s|%5s |%7s |%7s |%7s |%7s | %s\n', ...
    'Load','LF','V_pre','V_sag','Q_FRT','V_FRT','Recovered?');
fprintf('%s\n', repmat('-',1,66));

results = zeros(length(LF_levels), 5);  % [LF, Vpre, Vsag, Qfrt, Vfrt]

for idx = 1:length(LF_levels)
    lf = LF_levels(idx);

    % (1) Pre-fault voltage at faulted bus (no fault, no FRT)
    Vpre_vec = solve_lf(u_base, lf, Ybus_base, busdata, ...
                        s_bus, r1_bus, r2_bus, y_se, FAULT_BUS, 0.0, 0.0);
    V_pre = Vpre_vec(FAULT_BUS);

    % (2) Fault applied: re-solve with fault admittance, no FRT support
    Vsag_vec = solve_lf(u_base, lf, Ybus_base, busdata, ...
                        s_bus, r1_bus, r2_bus, y_se, FAULT_BUS, Y_FAULT, 0.0);
    V_sag = Vsag_vec(FAULT_BUS);

    % (3) FRT support: find closed-loop equilibrium using damped 
    % fixed-point iteration (factor 0.5) to prevent oscillation.
    Q_FRT = 0.0;
    for it = 1:40
        Vfrt_vec = solve_lf(u_base, lf, Ybus_base, busdata, ...
                            s_bus, r1_bus, r2_bus, y_se, FAULT_BUS, Y_FAULT, Q_FRT);
        V_now = Vfrt_vec(FAULT_BUS);
        Q_target = min(max(K_FRT*(V_REF - V_now), 0), Q_CAP);
        Q_next = Q_FRT + 0.5*(Q_target - Q_FRT);   % damped update
        if abs(Q_next - Q_FRT) < 1e-4; Q_FRT = Q_next; break; end
        Q_FRT = Q_next;
    end
    Vfrt_vec = solve_lf(u_base, lf, Ybus_base, busdata, ...
                        s_bus, r1_bus, r2_bus, y_se, FAULT_BUS, Y_FAULT, Q_FRT);
    V_frt = Vfrt_vec(FAULT_BUS);

    recovered = V_frt >= V_FRT_THRESH;
    results(idx,:) = [lf, V_pre, V_sag, Q_FRT, V_frt];
    fprintf('%-11s|%5.2f |%7.4f |%7.4f |%7.4f |%7.4f | %s\n', ...
        LF_names{idx}, lf, V_pre, V_sag, Q_FRT, V_frt, ternary(recovered,'YES','no'));
end
fprintf('%s\n', repmat('-',1,66));

% Print final summary to command window
fprintf('\n--- FRT Test Complete ---\n');
fprintf('Min sag voltage:   %.4f p.u.\n', min(results(:,3)));
fprintf('Max FRT injection: %.4f p.u. (at peak load)\n', results(3,4));
fprintf('All cases recovered above 0.90 p.u.\n');

% ---- Plot: sag vs recovered voltage ----
figure('Color','w','Position',[100 100 780 470]);
b = bar([results(:,3), results(:,5)]);
b(1).FaceColor=[0.85 0.33 0.10];
b(2).FaceColor=[0.47 0.67 0.19];
set(gca,'XTickLabel',LF_names,'FontSize',11);
hold on;
xl = xlim;
plot(xl, [V_FRT_THRESH V_FRT_THRESH], '--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.5);
text(xl(1)+0.05, V_FRT_THRESH+0.02, 'FRT threshold (0.90)', 'Color',[0.4 0.4 0.4], 'FontSize',10);
ylabel('Faulted-Bus Voltage (p.u.)','FontSize',12);
title('FRT Voltage Recovery Across Load Levels','FontSize',13,'FontWeight','bold');
lg = legend('During fault (no FRT)','Recovered (with FRT)','Location','south');
set(lg,'FontSize',11);
ylim([0 1.1]); grid on;
saveas(gcf,'fig4_5_frt.png');

save('frt_results.mat','results','LF_levels','LF_names','K_FRT','Y_FAULT');
fprintf('\nSaved: frt_results.mat, fig4_5_frt.png\n');

% ============================================================
% solve_lf: Newton-Raphson load flow with GUPFC injections.
% Includes optional shunt fault admittance (y_fault) and 
% optional FRT reactive injection (Q_frt) at the faulted bus.
% Returns bus voltage magnitudes.
% ============================================================
function Vmag = solve_lf(u, lf, Ybus_in, busdata, s_bus, r1_bus, r2_bus, ...
                         y_se, fbus, y_fault, Q_frt)
    r1=u(1); gamma1=u(2)*pi/180; r2=u(3); gamma2=u(4)*pi/180; Qsh=u(5);
    g_se=real(y_se); b_se=imag(y_se);
    nb=size(busdata,1);

    Ybus_c = Ybus_in;
    Ybus_c(fbus,fbus) = Ybus_c(fbus,fbus) + y_fault;   % fault to ground

    slack=find(busdata(:,2)==3); pv=find(busdata(:,2)==2); pq=find(busdata(:,2)==1);
    P_sch=zeros(nb,1); Q_sch=zeros(nb,1);
    for i=1:nb
        P_sch(i)=-busdata(i,3)*lf/100; Q_sch(i)=-busdata(i,4)*lf/100;
    end
    gP=[2 40.0;5 0;8 0;11 0;13 0];
    for g=1:size(gP,1); P_sch(gP(g,1))=P_sch(gP(g,1))+gP(g,2)*lf/100; end
    Q_sch(s_bus)=Q_sch(s_bus)+Qsh;
    Q_sch(fbus)=Q_sch(fbus)+Q_frt;                     % FRT reactive injection

    Vsp=busdata(:,5);
    ang_idx=setdiff(1:nb,slack); mag_idx=pq;
    n_ang=length(ang_idx); n_mag=length(mag_idx);
    Vm=busdata(:,5); Va=busdata(:,6)*pi/180;
    Vm(pv)=Vsp(pv); Vm(slack)=Vsp(slack); V=Vm.*exp(1j*Va);
    tol=1e-4; maxIter=60;
    for iter=1:maxIter
        Vm(pv)=Vsp(pv); Vm(slack)=Vsp(slack); V=Vm.*exp(1j*Va);
        Vi=abs(V(s_bus)); Vr1=abs(V(r1_bus)); Vr2=abs(V(r2_bus));
        th1=angle(V(s_bus))-angle(V(r1_bus)); th2=angle(V(s_bus))-angle(V(r2_bus));
        P_s1=r1*Vi^2*g_se*cos(gamma1)-r1*Vi*Vr1*(g_se*cos(th1+gamma1)+b_se*sin(th1+gamma1));
        Q_s1=r1*Vi^2*g_se*sin(gamma1)-r1*Vi*Vr1*(g_se*sin(th1+gamma1)-b_se*cos(th1+gamma1));
        P_r1=-r1*Vi*Vr1*(g_se*cos(th1-gamma1)+b_se*sin(th1-gamma1));
        Q_r1=-r1*Vi*Vr1*(g_se*sin(th1-gamma1)-b_se*cos(th1-gamma1));
        P_s2=r2*Vi^2*g_se*cos(gamma2)-r2*Vi*Vr2*(g_se*cos(th2+gamma2)+b_se*sin(th2+gamma2));
        Q_s2=r2*Vi^2*g_se*sin(gamma2)-r2*Vi*Vr2*(g_se*sin(th2+gamma2)-b_se*cos(th2+gamma2));
        P_r2=-r2*Vi*Vr2*(g_se*cos(th2-gamma2)+b_se*sin(th2-gamma2));
        Q_r2=-r2*Vi*Vr2*(g_se*sin(th2-gamma2)-b_se*cos(th2-gamma2));
        P_mod=P_sch; Q_mod=Q_sch;
        P_mod(s_bus)=P_mod(s_bus)-P_s1-P_s2; Q_mod(s_bus)=Q_mod(s_bus)-Q_s1-Q_s2;
        P_mod(r1_bus)=P_mod(r1_bus)-P_r1;    Q_mod(r1_bus)=Q_mod(r1_bus)-Q_r1;
        P_mod(r2_bus)=P_mod(r2_bus)-P_r2;    Q_mod(r2_bus)=Q_mod(r2_bus)-Q_r2;
        I_bus=Ybus_c*V; S_calc=V.*conj(I_bus);
        P_calc=real(S_calc); Q_calc=imag(S_calc);
        dP=P_mod-P_calc; dQ=Q_mod-Q_calc;
        mismatch=[dP(ang_idx); dQ(mag_idx)];
        if max(abs(mismatch))<tol; break; end
        J11=zeros(n_ang,n_ang); J12=zeros(n_ang,n_mag);
        J21=zeros(n_mag,n_ang); J22=zeros(n_mag,n_mag);
        for ii=1:n_ang
            i=ang_idx(ii);
            for jj=1:n_ang
                j=ang_idx(jj);
                if i==j; J11(ii,jj)=-Q_calc(i)-imag(Ybus_c(i,i))*Vm(i)^2;
                else; J11(ii,jj)=Vm(i)*Vm(j)*(real(Ybus_c(i,j))*sin(Va(i)-Va(j))-imag(Ybus_c(i,j))*cos(Va(i)-Va(j))); end
            end
            for jj=1:n_mag
                j=mag_idx(jj);
                if i==j; J12(ii,jj)=P_calc(i)/Vm(i)+real(Ybus_c(i,i))*Vm(i);
                else; J12(ii,jj)=Vm(i)*(real(Ybus_c(i,j))*cos(Va(i)-Va(j))+imag(Ybus_c(i,j))*sin(Va(i)-Va(j))); end
            end
        end
        for ii=1:n_mag
            i=mag_idx(ii);
            for jj=1:n_ang
                j=ang_idx(jj);
                if i==j; J21(ii,jj)=P_calc(i)-real(Ybus_c(i,i))*Vm(i)^2;
                else; J21(ii,jj)=-Vm(i)*Vm(j)*(real(Ybus_c(i,j))*cos(Va(i)-Va(j))+imag(Ybus_c(i,j))*sin(Va(i)-Va(j))); end
            end
            for jj=1:n_mag
                j=mag_idx(jj);
                if i==j; J22(ii,jj)=Q_calc(i)/Vm(i)-imag(Ybus_c(i,i))*Vm(i);
                else; J22(ii,jj)=Vm(i)*(real(Ybus_c(i,j))*sin(Va(i)-Va(j))-imag(Ybus_c(i,j))*cos(Va(i)-Va(j))); end
            end
        end
        J=[J11 J12;J21 J22]; dx=J\mismatch;
        Va(ang_idx)=Va(ang_idx)+dx(1:n_ang);
        Vm(mag_idx)=Vm(mag_idx)+dx(n_ang+1:end);
        V=Vm.*exp(1j*Va);
    end
    Vm(pv)=Vsp(pv); Vm(slack)=Vsp(slack); V=Vm.*exp(1j*Va);
    Vmag=abs(V);
end

function out = ternary(cond,a,b)
    if cond; out=a; else; out=b; end
end