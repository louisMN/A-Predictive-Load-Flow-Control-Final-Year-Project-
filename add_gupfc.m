% Step 5: GUPFC at bus triplet (3, 7, 17)
% Paper creates NEW virtual branches — no existing line needed
% Z_se = series converter impedance

s_bus  = 3;    % sending bus
r1_bus = 7;    % receiving bus 1
r2_bus = 17;   % receiving bus 2

% GUPFC control variables
r1     = 0.05;
gamma1 = 30 * pi/180;
r2     = 0.05;
gamma2 = 30 * pi/180;
Qsh    = 0.10;

% Series converter impedance (GUPFC creates new branches with this impedance)
R_se = 0.01;
X_se = 0.15;
Z_se = R_se + 1j*X_se;
y_se = 1/Z_se;
g_se = real(y_se);
b_se = imag(y_se);

fprintf('GUPFC series admittance: g_se=%.4f, b_se=%.4f\n', g_se, b_se);

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

% ADD new GUPFC virtual branches to Y-bus
% Branch 1: bus 3 to bus 7
Ybus_cont(s_bus,s_bus)   = Ybus_cont(s_bus,s_bus)   + y_se;
Ybus_cont(r1_bus,r1_bus) = Ybus_cont(r1_bus,r1_bus) + y_se;
Ybus_cont(s_bus,r1_bus)  = Ybus_cont(s_bus,r1_bus)  - y_se;
Ybus_cont(r1_bus,s_bus)  = Ybus_cont(r1_bus,s_bus)  - y_se;

% Branch 2: bus 3 to bus 17
Ybus_cont(s_bus,s_bus)   = Ybus_cont(s_bus,s_bus)   + y_se;
Ybus_cont(r2_bus,r2_bus) = Ybus_cont(r2_bus,r2_bus) + y_se;
Ybus_cont(s_bus,r2_bus)  = Ybus_cont(s_bus,r2_bus)  - y_se;
Ybus_cont(r2_bus,s_bus)  = Ybus_cont(r2_bus,s_bus)  - y_se;

% Scheduled power
P_sch=zeros(nb,1); Q_sch=zeros(nb,1);
for i=1:nb
    P_sch(i)=-busdata(i,3)/100;
    Q_sch(i)=-busdata(i,4)/100;
end
gendata_local=[2 40.0;5 0;8 0;11 0;13 0];
for g=1:size(gendata_local,1)
    P_sch(gendata_local(g,1))=P_sch(gendata_local(g,1))+gendata_local(g,2)/100;
end

% Add shunt injection at sending bus
Q_sch(s_bus) = Q_sch(s_bus) + Qsh;

% Initial voltages
Vm=busdata(:,5); Va=busdata(:,6)*pi/180;
V=Vm.*exp(1j*Va);
slack=find(busdata(:,2)==3);
pq=find(busdata(:,2)==1);
ang_idx=setdiff(1:nb,slack);
mag_idx=pq;
n_ang=length(ang_idx); n_mag=length(mag_idx);

tol=1e-4; maxIter=50;

for iter=1:maxIter
    % Current voltages
    Vi  = abs(V(s_bus));
    Vr1 = abs(V(r1_bus));
    Vr2 = abs(V(r2_bus));
    th1 = angle(V(s_bus)) - angle(V(r1_bus));  % theta_s - theta_r1
    th2 = angle(V(s_bus)) - angle(V(r2_bus));  % theta_s - theta_r2

    % GUPFC power injections at sending bus (line 1)
    P_s1 = r1*Vi^2*g_se*cos(gamma1) ...
         - r1*Vi*Vr1*(g_se*cos(th1+gamma1) + b_se*sin(th1+gamma1));
    Q_s1 = r1*Vi^2*g_se*sin(gamma1) ...
         - r1*Vi*Vr1*(g_se*sin(th1+gamma1) - b_se*cos(th1+gamma1));

    % GUPFC power injections at receiving bus 1
    P_r1 = -r1*Vi*Vr1*(g_se*cos(th1-gamma1) + b_se*sin(th1-gamma1));
    Q_r1 = -r1*Vi*Vr1*(g_se*sin(th1-gamma1) - b_se*cos(th1-gamma1));

    % GUPFC power injections at sending bus (line 2)
    P_s2 = r2*Vi^2*g_se*cos(gamma2) ...
         - r2*Vi*Vr2*(g_se*cos(th2+gamma2) + b_se*sin(th2+gamma2));
    Q_s2 = r2*Vi^2*g_se*sin(gamma2) ...
         - r2*Vi*Vr2*(g_se*sin(th2+gamma2) - b_se*cos(th2+gamma2));

    % GUPFC power injections at receiving bus 2
    P_r2 = -r2*Vi*Vr2*(g_se*cos(th2-gamma2) + b_se*sin(th2-gamma2));
    Q_r2 = -r2*Vi*Vr2*(g_se*sin(th2-gamma2) - b_se*cos(th2-gamma2));

    % Modified scheduled power
    P_mod = P_sch; Q_mod = Q_sch;
    P_mod(s_bus)  = P_mod(s_bus)  - P_s1 - P_s2;
    Q_mod(s_bus)  = Q_mod(s_bus)  - Q_s1 - Q_s2;
    P_mod(r1_bus) = P_mod(r1_bus) - P_r1;
    Q_mod(r1_bus) = Q_mod(r1_bus) - Q_r1;
    P_mod(r2_bus) = P_mod(r2_bus) - P_r2;
    Q_mod(r2_bus) = Q_mod(r2_bus) - Q_r2;

    % NR iteration
    I_bus=Ybus_cont*V; S_calc=V.*conj(I_bus);
    P_calc=real(S_calc); Q_calc=imag(S_calc);
    dP=P_mod-P_calc; dQ=Q_mod-Q_calc;
    mismatch=[dP(ang_idx); dQ(mag_idx)];
    if max(abs(mismatch))<tol
        fprintf('GUPFC NR converged in %d iterations\n',iter);
        break;
    end

    J11=zeros(n_ang,n_ang); J12=zeros(n_ang,n_mag);
    J21=zeros(n_mag,n_ang); J22=zeros(n_mag,n_mag);
    for ii=1:n_ang
        i=ang_idx(ii);
        for jj=1:n_ang
            j=ang_idx(jj);
            if i==j
                J11(ii,jj)=-Q_calc(i)-imag(Ybus_cont(i,i))*Vm(i)^2;
            else
                J11(ii,jj)=Vm(i)*Vm(j)*(real(Ybus_cont(i,j))*sin(Va(i)-Va(j))...
                          -imag(Ybus_cont(i,j))*cos(Va(i)-Va(j)));
            end
        end
        for jj=1:n_mag
            j=mag_idx(jj);
            if i==j
                J12(ii,jj)=P_calc(i)/Vm(i)+real(Ybus_cont(i,i))*Vm(i);
            else
                J12(ii,jj)=Vm(i)*(real(Ybus_cont(i,j))*cos(Va(i)-Va(j))...
                          +imag(Ybus_cont(i,j))*sin(Va(i)-Va(j)));
            end
        end
    end
    for ii=1:n_mag
        i=mag_idx(ii);
        for jj=1:n_ang
            j=ang_idx(jj);
            if i==j
                J21(ii,jj)=P_calc(i)-real(Ybus_cont(i,i))*Vm(i)^2;
            else
                J21(ii,jj)=-Vm(i)*Vm(j)*(real(Ybus_cont(i,j))*cos(Va(i)-Va(j))...
                           +imag(Ybus_cont(i,j))*sin(Va(i)-Va(j)));
            end
        end
        for jj=1:n_mag
            j=mag_idx(jj);
            if i==j
                J22(ii,jj)=Q_calc(i)/Vm(i)-imag(Ybus_cont(i,i))*Vm(i);
            else
                J22(ii,jj)=Vm(i)*(real(Ybus_cont(i,j))*sin(Va(i)-Va(j))...
                          -imag(Ybus_cont(i,j))*cos(Va(i)-Va(j)));
            end
        end
    end
    J=[J11 J12;J21 J22];
    dx=J\mismatch;
    Va(ang_idx)=Va(ang_idx)+dx(1:n_ang);
    Vm(mag_idx)=Vm(mag_idx)+dx(n_ang+1:end);
    V=Vm.*exp(1j*Va);
end

% Calculate losses
P_loss_g=0; Q_loss_g=0;
for k=1:size(linedata_cont,1)
    from=linedata_cont(k,1); to=linedata_cont(k,2);
    R=linedata_cont(k,3); X=linedata_cont(k,4);
    if R==0&&X==0; continue; end
    z=R+1j*X; I_line=(V(from)-V(to))/z;
    S_loss=I_line*conj(I_line)*z;
    P_loss_g=P_loss_g+real(S_loss);
    Q_loss_g=Q_loss_g+imag(S_loss);
end
P_loss_g=P_loss_g*100;
V_dev_g=sum(abs(abs(V)-1.0));

fprintf('\n===== RESULTS WITH GUPFC AT (3, 7, 17) =====\n');
fprintf('Active Loss with GUPFC:    %.4f MW\n', P_loss_g);
fprintf('Voltage Deviation:         %.4f p.u.\n', V_dev_g);
fprintf('\n--- Comparison ---\n');
fprintf('Base case (no outage):     17.73 MW\n');
fprintf('Outage, no GUPFC:          32.95 MW\n');
fprintf('Outage + GUPFC:            %.4f MW\n', P_loss_g);
fprintf('Loss reduction from GUPFC: %.2f%%\n', (32.9533-P_loss_g)/32.9533*100);