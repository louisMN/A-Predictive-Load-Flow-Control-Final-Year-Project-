function [Ploss, Vdev] = run_gupfc_lf(u, lf, P_nom, Q_nom, Ybus_c, ...
                                        busdata, linedata_c, s_bus, ...
                                        r1_bus, r2_bus, y_se)
% Runs NR load flow with GUPFC and returns losses and voltage deviation.
% Generators (buses 2,5,8,11,13) treated as PROPER PV buses:
%   voltage magnitude HELD at setpoint, reactive power Q solved as free var.
% Real generation scales with load factor so slack does not absorb all growth.
%
% u = [r1, gamma1_deg, r2, gamma2_deg, Qsh]

r1     = u(1);
gamma1 = u(2)*pi/180;
r2     = u(3);
gamma2 = u(4)*pi/180;
Qsh    = u(5);

g_se = real(y_se);
b_se = imag(y_se);

nb = size(busdata,1);

% ---- Bus type sets (from busdata column 2: 3=slack, 2=PV, 1=PQ) ----
slack = find(busdata(:,2)==3);
pv    = find(busdata(:,2)==2);
pq    = find(busdata(:,2)==1);

% ---- Scheduled power: loads scale with load factor lf ----
P_sch = zeros(nb,1);
Q_sch = zeros(nb,1);
for i = 1:nb
    P_sch(i) = -busdata(i,3)*lf/100;   % load P scaled by lf
    Q_sch(i) = -busdata(i,4)*lf/100;   % load Q scaled by lf
end

% ---- Generator REAL power: scale with lf so total gen tracks total load ----
% (Q is NOT scheduled here - PV buses solve for their own Q.)
gendata_local = [2 40.0; 5 0; 8 0; 11 0; 13 0];   % [bus, P_MW]
for g = 1:size(gendata_local,1)
    bg = gendata_local(g,1);
    P_sch(bg) = P_sch(bg) + gendata_local(g,2)*lf/100;   % gen P scales with lf
end

% ---- GUPFC shunt reactive injection at sending bus ----
Q_sch(s_bus) = Q_sch(s_bus) + Qsh;

% ---- PV voltage setpoints: HELD constant throughout iterations ----
% Use the scheduled magnitudes from busdata as the PV setpoints.
Vsp = busdata(:,5);

% ---- Index sets ----
ang_idx = setdiff(1:nb, slack);   % angle unknown at all non-slack buses
mag_idx = pq;                     % magnitude unknown only at PQ buses
n_ang   = length(ang_idx);
n_mag   = length(mag_idx);

% ---- Initial state; PV + slack magnitudes pinned to setpoint ----
Vm = busdata(:,5);
Va = busdata(:,6)*pi/180;
Vm(pv)    = Vsp(pv);      % enforce PV setpoint
Vm(slack) = Vsp(slack);  % enforce slack setpoint
V = Vm.*exp(1j*Va);

tol=1e-4; maxIter=30;

for iter = 1:maxIter
    % Re-pin PV/slack magnitudes every iteration (they never drift)
    Vm(pv)    = Vsp(pv);
    Vm(slack) = Vsp(slack);
    V = Vm.*exp(1j*Va);

    Vi  = abs(V(s_bus));
    Vr1 = abs(V(r1_bus));
    Vr2 = abs(V(r2_bus));
    th1 = angle(V(s_bus)) - angle(V(r1_bus));
    th2 = angle(V(s_bus)) - angle(V(r2_bus));

    % GUPFC series-converter power injections
    P_s1 = r1*Vi^2*g_se*cos(gamma1) - r1*Vi*Vr1*(g_se*cos(th1+gamma1)+b_se*sin(th1+gamma1));
    Q_s1 = r1*Vi^2*g_se*sin(gamma1) - r1*Vi*Vr1*(g_se*sin(th1+gamma1)-b_se*cos(th1+gamma1));
    P_r1 = -r1*Vi*Vr1*(g_se*cos(th1-gamma1)+b_se*sin(th1-gamma1));
    Q_r1 = -r1*Vi*Vr1*(g_se*sin(th1-gamma1)-b_se*cos(th1-gamma1));
    P_s2 = r2*Vi^2*g_se*cos(gamma2) - r2*Vi*Vr2*(g_se*cos(th2+gamma2)+b_se*sin(th2+gamma2));
    Q_s2 = r2*Vi^2*g_se*sin(gamma2) - r2*Vi*Vr2*(g_se*sin(th2+gamma2)-b_se*cos(th2+gamma2));
    P_r2 = -r2*Vi*Vr2*(g_se*cos(th2-gamma2)+b_se*sin(th2-gamma2));
    Q_r2 = -r2*Vi*Vr2*(g_se*sin(th2-gamma2)-b_se*cos(th2-gamma2));

    P_mod = P_sch; Q_mod = Q_sch;
    P_mod(s_bus)  = P_mod(s_bus)  - P_s1 - P_s2;
    Q_mod(s_bus)  = Q_mod(s_bus)  - Q_s1 - Q_s2;
    P_mod(r1_bus) = P_mod(r1_bus) - P_r1;
    Q_mod(r1_bus) = Q_mod(r1_bus) - Q_r1;
    P_mod(r2_bus) = P_mod(r2_bus) - P_r2;
    Q_mod(r2_bus) = Q_mod(r2_bus) - Q_r2;

    I_bus  = Ybus_c*V;
    S_calc = V.*conj(I_bus);
    P_calc = real(S_calc);
    Q_calc = imag(S_calc);
    dP = P_mod - P_calc;
    dQ = Q_mod - Q_calc;
    mismatch = [dP(ang_idx); dQ(mag_idx)];
    if max(abs(mismatch)) < tol; break; end

    J11=zeros(n_ang,n_ang); J12=zeros(n_ang,n_mag);
    J21=zeros(n_mag,n_ang); J22=zeros(n_mag,n_mag);
    for ii=1:n_ang
        i=ang_idx(ii);
        for jj=1:n_ang
            j=ang_idx(jj);
            if i==j; J11(ii,jj)=-Q_calc(i)-imag(Ybus_c(i,i))*Vm(i)^2;
            else; J11(ii,jj)=Vm(i)*Vm(j)*(real(Ybus_c(i,j))*sin(Va(i)-Va(j))-imag(Ybus_c(i,j))*cos(Va(i)-Va(j)));
            end
        end
        for jj=1:n_mag
            j=mag_idx(jj);
            if i==j; J12(ii,jj)=P_calc(i)/Vm(i)+real(Ybus_c(i,i))*Vm(i);
            else; J12(ii,jj)=Vm(i)*(real(Ybus_c(i,j))*cos(Va(i)-Va(j))+imag(Ybus_c(i,j))*sin(Va(i)-Va(j)));
            end
        end
    end
    for ii=1:n_mag
        i=mag_idx(ii);
        for jj=1:n_ang
            j=ang_idx(jj);
            if i==j; J21(ii,jj)=P_calc(i)-real(Ybus_c(i,i))*Vm(i)^2;
            else; J21(ii,jj)=-Vm(i)*Vm(j)*(real(Ybus_c(i,j))*cos(Va(i)-Va(j))+imag(Ybus_c(i,j))*sin(Va(i)-Va(j)));
            end
        end
        for jj=1:n_mag
            j=mag_idx(jj);
            if i==j; J22(ii,jj)=Q_calc(i)/Vm(i)-imag(Ybus_c(i,i))*Vm(i);
            else; J22(ii,jj)=Vm(i)*(real(Ybus_c(i,j))*sin(Va(i)-Va(j))-imag(Ybus_c(i,j))*cos(Va(i)-Va(j)));
            end
        end
    end
    J=[J11 J12;J21 J22];
    dx=J\mismatch;
    Va(ang_idx)=Va(ang_idx)+dx(1:n_ang);
    Vm(mag_idx)=Vm(mag_idx)+dx(n_ang+1:end);
    V=Vm.*exp(1j*Va);
end

% Final re-pin (safety) before computing metrics
Vm(pv)    = Vsp(pv);
Vm(slack) = Vsp(slack);
V = Vm.*exp(1j*Va);

% ---- Line losses ----
Ploss=0;
for k=1:size(linedata_c,1)
    from=linedata_c(k,1); to=linedata_c(k,2);
    R=linedata_c(k,3); X=linedata_c(k,4);
    if R==0&&X==0; continue; end
    z=R+1j*X;
    I_line=(V(from)-V(to))/z;
    S_loss=I_line*conj(I_line)*z;
    Ploss=Ploss+real(S_loss);
end
Ploss = Ploss*100;

% ---- Voltage deviation (paper Eq. 21 / Eq. 31) ----
% Worst-case deviation among CONTROLLABLE buses (PQ load buses).
% Generator (PV) and slack bus magnitudes are fixed by setpoint and
% lie outside GUPFC control, so they are excluded from the metric.
Vmag = abs(V);
dev_all = abs(Vmag - 1.0);
dev_all(pv)    = 0;   % exclude generator setpoints
dev_all(slack) = 0;   % exclude slack setpoint
Vdev = max(dev_all);
end