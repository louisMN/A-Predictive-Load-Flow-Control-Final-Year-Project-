% Step 4: N-1 Contingency — Remove Branch 5 (Bus 2 to Bus 5)

% Branch 5 is line from bus 2 to bus 5 (row 5 in linedata)
critical_branch = 5;

fprintf('Removing branch %d (Bus %d - Bus %d)\n', ...
    critical_branch, linedata(critical_branch,1), linedata(critical_branch,2));

% Create modified line data with branch 5 removed
linedata_contingency = linedata;
linedata_contingency(critical_branch, :) = [];

% Rebuild Y-bus with branch removed
nb = size(busdata,1);
Ybus_cont = zeros(nb,nb);

for k = 1:size(linedata_contingency,1)
    from = linedata_contingency(k,1);
    to   = linedata_contingency(k,2);
    R    = linedata_contingency(k,3);
    X    = linedata_contingency(k,4);
    B    = linedata_contingency(k,5);
    
    if R == 0
        y = 1/(1j*X);
    else
        y = 1/(R + 1j*X);
    end
    yc = 1j*B/2;
    
    Ybus_cont(from,from) = Ybus_cont(from,from) + y + yc;
    Ybus_cont(to,to)     = Ybus_cont(to,to)     + y + yc;
    Ybus_cont(from,to)   = Ybus_cont(from,to)   - y;
    Ybus_cont(to,from)   = Ybus_cont(to,from)   - y;
end

% Run NR load flow with contingency Y-bus
Ybus_saved = Ybus;
Ybus = Ybus_cont;

tol=1e-4; maxIter=50; nb=size(busdata,1);
slack=find(busdata(:,2)==3);
pv=find(busdata(:,2)==2);
pq=find(busdata(:,2)==1);

P_sch=zeros(nb,1); Q_sch=zeros(nb,1);
for i=1:nb
    P_sch(i)=-busdata(i,3)/100;
    Q_sch(i)=-busdata(i,4)/100;
end
gendata=[2 40.0;5 0;8 0;11 0;13 0];
for g=1:size(gendata,1)
    P_sch(gendata(g,1))=P_sch(gendata(g,1))+gendata(g,2)/100;
end

Vm=busdata(:,5); Va=busdata(:,6)*pi/180;
V=Vm.*exp(1j*Va);
ang_idx=setdiff(1:nb,slack);
mag_idx=pq;
n_ang=length(ang_idx); n_mag=length(mag_idx);

converged=false;
for iter=1:maxIter
    I_bus=Ybus*V; S_calc=V.*conj(I_bus);
    P_calc=real(S_calc); Q_calc=imag(S_calc);
    dP=P_sch-P_calc; dQ=Q_sch-Q_calc;
    dP_red=dP(ang_idx); dQ_red=dQ(mag_idx);
    mismatch=[dP_red;dQ_red];
    if max(abs(mismatch))<tol
        fprintf('Contingency NR converged in %d iterations\n',iter);
        converged=true; break;
    end
    J11=zeros(n_ang,n_ang); J12=zeros(n_ang,n_mag);
    J21=zeros(n_mag,n_ang); J22=zeros(n_mag,n_mag);
    for ii=1:n_ang
        i=ang_idx(ii);
        for jj=1:n_ang
            j=ang_idx(jj);
            if i==j; J11(ii,jj)=-Q_calc(i)-imag(Ybus(i,i))*Vm(i)^2;
            else; J11(ii,jj)=Vm(i)*Vm(j)*(real(Ybus(i,j))*sin(Va(i)-Va(j))-imag(Ybus(i,j))*cos(Va(i)-Va(j)));
            end
        end
        for jj=1:n_mag
            j=mag_idx(jj);
            if i==j; J12(ii,jj)=P_calc(i)/Vm(i)+real(Ybus(i,i))*Vm(i);
            else; J12(ii,jj)=Vm(i)*(real(Ybus(i,j))*cos(Va(i)-Va(j))+imag(Ybus(i,j))*sin(Va(i)-Va(j)));
            end
        end
    end
    for ii=1:n_mag
        i=mag_idx(ii);
        for jj=1:n_ang
            j=ang_idx(jj);
            if i==j; J21(ii,jj)=P_calc(i)-real(Ybus(i,i))*Vm(i)^2;
            else; J21(ii,jj)=-Vm(i)*Vm(j)*(real(Ybus(i,j))*cos(Va(i)-Va(j))+imag(Ybus(i,j))*sin(Va(i)-Va(j)));
            end
        end
        for jj=1:n_mag
            j=mag_idx(jj);
            if i==j; J22(ii,jj)=Q_calc(i)/Vm(i)-imag(Ybus(i,i))*Vm(i);
            else; J22(ii,jj)=Vm(i)*(real(Ybus(i,j))*sin(Va(i)-Va(j))-imag(Ybus(i,j))*cos(Va(i)-Va(j)));
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
P_loss_cont=0; Q_loss_cont=0;
for k=1:size(linedata_contingency,1)
    from=linedata_contingency(k,1); to=linedata_contingency(k,2);
    R=linedata_contingency(k,3); X=linedata_contingency(k,4);
    if R==0&&X==0; continue; end
    z=R+1j*X; I_line=(V(from)-V(to))/z;
    S_loss=I_line*conj(I_line)*z;
    P_loss_cont=P_loss_cont+real(S_loss);
    Q_loss_cont=Q_loss_cont+imag(S_loss);
end
P_loss_cont=P_loss_cont*100;
Q_loss_cont=Q_loss_cont*100;

% Voltage deviations
V_base = busdata(:,5);
V_cont = abs(V);
V_dev = sum(abs(V_cont - 1.0));

fprintf('\n===== CONTINGENCY RESULTS (Branch 5 removed) =====\n');
fprintf('Total Active Loss:     %.4f MW\n', P_loss_cont);
fprintf('Total Reactive Loss:   %.4f MVAR\n', Q_loss_cont);
fprintf('Total Voltage Deviation: %.4f p.u.\n', V_dev);

% Violations
fprintf('\nVoltage violations (outside 0.95-1.05 pu):\n');
violations = 0;
for i=1:nb
    if abs(V(i))<0.95 || abs(V(i))>1.05
        fprintf('  Bus %2d: %.4f pu\n', i, abs(V(i)));
        violations = violations+1;
    end
end
if violations==0; fprintf('  None\n'); end

% Restore original Ybus
Ybus = Ybus_saved;