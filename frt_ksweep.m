% ============================================================
% frt_ksweep.m
% Tests FRT voltage recovery for different K_FRT gains (2, 3, 4, 5).
% Uses the same damped fixed-point iteration as frt_test.m.
% Requires: load_data.m, build_Ybus.m
% ============================================================

V_FRT_THRESH = 0.90;
V_REF        = 1.00;
FAULT_BUS    = 30;
Y_FAULT      = 0.45;        % fault shunt admittance (same as frt_test.m)
Q_CAP        = 1.5;         % converter reactive rating (saturation)

K_values  = [2.0, 3.0, 4.0, 5.0];
LF_levels = [0.65, 1.00, 1.20, 1.10];
LF_names  = {'Off-peak','Nominal','Peak','Transition'};

s_bus=3; r1_bus=7; r2_bus=17;
R_se=0.01; X_se=0.15; y_se=1/(R_se+1j*X_se);

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

% Pre-fault and faulted (no-FRT) voltages per load level -- solved once
V_pre = zeros(1,length(LF_levels));
V_sag = zeros(1,length(LF_levels));
for idx = 1:length(LF_levels)
    Vp = solve_lf(u_base, LF_levels(idx), Ybus_base, busdata, s_bus, r1_bus, r2_bus, y_se, FAULT_BUS, 0.0, 0.0);
    V_pre(idx) = Vp(FAULT_BUS);
    Vs = solve_lf(u_base, LF_levels(idx), Ybus_base, busdata, s_bus, r1_bus, r2_bus, y_se, FAULT_BUS, Y_FAULT, 0.0);
    V_sag(idx) = Vs(FAULT_BUS);
end

fprintf('\n========= FRT K-FACTOR SENSITIVITY SWEEP =========\n');
fprintf('Faulted bus: %d | Y_fault: %.2f | Q cap: %.1f p.u.\n\n', ...
        FAULT_BUS, Y_FAULT, Q_CAP);

Vrec = zeros(length(K_values), length(LF_levels));
Qfrt = zeros(length(K_values), length(LF_levels));

for ki = 1:length(K_values)
    K = K_values(ki);
    fprintf('--- K_FRT = %.1f ---\n', K);
    fprintf('%-11s | V_sag  | Q_FRT | V_FRT  | Recovered?\n','Load');
    fprintf('%s\n', repmat('-',1,50));
    for idx = 1:length(LF_levels)
        lf_now = LF_levels(idx);
        Q_FRT = 0.0;
        for it = 1:40
            Vf = solve_lf(u_base, lf_now, Ybus_base, busdata, ...
                          s_bus, r1_bus, r2_bus, y_se, FAULT_BUS, Y_FAULT, Q_FRT);
            V_now = Vf(FAULT_BUS);
            Q_target = min(max(K*(V_REF - V_now), 0), Q_CAP);
            Q_next = Q_FRT + 0.5*(Q_target - Q_FRT);   % damped update
            if abs(Q_next - Q_FRT) < 1e-4; Q_FRT = Q_next; break; end
            Q_FRT = Q_next;
        end
        Vf = solve_lf(u_base, lf_now, Ybus_base, busdata, ...
                      s_bus, r1_bus, r2_bus, y_se, FAULT_BUS, Y_FAULT, Q_FRT);
        Vfrt = Vf(FAULT_BUS);
        Vrec(ki,idx) = Vfrt; Qfrt(ki,idx) = Q_FRT;
        rec = ternary(Vfrt>=V_FRT_THRESH,'YES','no');
        fprintf('%-11s | %.4f | %.3f | %.4f | %s\n', ...
            LF_names{idx}, V_sag(idx), Q_FRT, Vfrt, rec);
    end
    fprintf('\n');
end

% ---- Plot: recovered voltage vs K ----
figure('Color','w','Position',[100 100 780 470]);
markers = {'o-','s-','^-','d-'};
cols = [0.00 0.45 0.74; 0.47 0.67 0.19; 0.85 0.33 0.10; 0.49 0.18 0.56];
hold on;
for idx = 1:length(LF_levels)
    plot(K_values, Vrec(:,idx), markers{idx}, 'Color',cols(idx,:), ...
        'LineWidth',2,'MarkerSize',7,'MarkerFaceColor',cols(idx,:));
end
plot([1.8 5.2], [V_FRT_THRESH V_FRT_THRESH], '--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.5);
text(1.85, V_FRT_THRESH+0.005, 'Recovery threshold (0.90)', 'Color',[0.4 0.4 0.4], 'FontSize',10);
plot([3.0 3.0], ylim, ':', 'Color',[0.6 0.2 0.2], 'LineWidth',1.5);
text(3.05, min(ylim)+0.01, 'Chosen K = 3.0', 'Color',[0.6 0.2 0.2], 'FontSize',10);
hold off; grid on;
xlabel('FRT Gain K_{FRT}','FontSize',12);
ylabel('Recovered Faulted-Bus Voltage (p.u.)','FontSize',12);
title('FRT Recovery vs K-Factor Across Load Levels','FontSize',13,'FontWeight','bold');
lg = legend(LF_names,'Location','southeast');
set(lg,'FontSize',10);
xlim([1.8 5.2]); set(gca,'XTick',K_values,'FontSize',11);
saveas(gcf,'fig4_6_ksweep.png');

save('frt_ksweep_results.mat','K_values','LF_levels','LF_names','Vrec','Qfrt','V_pre','V_sag');
fprintf('Saved: frt_ksweep_results.mat, fig4_6_ksweep.png\n');

% ---- solve_lf: identical to frt_test.m (NR + GUPFC + fault + FRT Q) ----
function Vmag = solve_lf(u, lf, Ybus_in, busdata, s_bus, r1_bus, r2_bus, ...
                         y_se, fbus, y_fault, Q_frt)
    r1=u(1); gamma1=u(2)*pi/180; r2=u(3); gamma2=u(4)*pi/180; Qsh=u(5);
    g_se=real(y_se); b_se=imag(y_se);
    nb=size(busdata,1);
    Ybus_c = Ybus_in;
    Ybus_c(fbus,fbus) = Ybus_c(fbus,fbus) + y_fault;
    slack=find(busdata(:,2)==3); pv=find(busdata(:,2)==2); pq=find(busdata(:,2)==1);
    P_sch=zeros(nb,1); Q_sch=zeros(nb,1);
    for i=1:nb
        P_sch(i)=-busdata(i,3)*lf/100; Q_sch(i)=-busdata(i,4)*lf/100;
    end
    gP=[2 40.0;5 0;8 0;11 0;13 0];
    for g=1:size(gP,1); P_sch(gP(g,1))=P_sch(gP(g,1))+gP(g,2)*lf/100; end
    Q_sch(s_bus)=Q_sch(s_bus)+Qsh;
    Q_sch(fbus)=Q_sch(fbus)+Q_frt;
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