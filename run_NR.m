% Step 3: Newton-Raphson Load Flow

tol     = 1e-4;  % convergence tolerance
maxIter = 50;    % max iterations
nb      = size(busdata, 1);

% Extract bus types
slack = find(busdata(:,2) == 3);   % slack bus (bus 1)
pv    = find(busdata(:,2) == 2);   % generator buses
pq    = find(busdata(:,2) == 1);   % load buses

% Scheduled power (convert MW to per unit)
P_sch = zeros(nb,1);
Q_sch = zeros(nb,1);

for i = 1:nb
    P_sch(i) = -busdata(i,3) / 100;
    Q_sch(i) = -busdata(i,4) / 100;
end

% Add generator injections for PV buses
gendata = [
2  40.0   50.0
5  0      37.0
8  0      37.3
11 0      16.2
13 0      10.6
];

for g = 1:size(gendata,1)
    bus_g = gendata(g,1);
    P_sch(bus_g) = P_sch(bus_g) + gendata(g,2)/100;
    Q_sch(bus_g) = Q_sch(bus_g) + gendata(g,3)/100;
end

% Initial voltage guess from bus data
Vm = busdata(:,5);         % voltage magnitude
Va = busdata(:,6)*pi/180;  % voltage angle in radians
V  = Vm .* exp(1j*Va);

% Buses where angle is unknown (all except slack)
ang_idx = setdiff(1:nb, slack);
% Buses where magnitude is unknown (PQ buses only)
mag_idx = pq;

n_ang = length(ang_idx);
n_mag = length(mag_idx);

fprintf('Starting NR iterations...\n');

for iter = 1:maxIter

    % Calculate power from voltages
    I_bus  = Ybus * V;
    S_calc = V .* conj(I_bus);
    P_calc = real(S_calc);
    Q_calc = imag(S_calc);

    % Mismatches
    dP = P_sch - P_calc;
    dQ = Q_sch - Q_calc;

    % Remove slack from dP, remove slack+PV from dQ
    dP_red = dP(ang_idx);
    dQ_red = dQ(mag_idx);
    mismatch = [dP_red; dQ_red];

    % Check convergence
    if max(abs(mismatch)) < tol
        fprintf('Converged in %d iterations\n', iter);
        converged = true;
        break;
    end

    % Build Jacobian
    J11 = zeros(n_ang, n_ang);
    J12 = zeros(n_ang, n_mag);
    J21 = zeros(n_mag, n_ang);
    J22 = zeros(n_mag, n_mag);

    for ii = 1:n_ang
        i = ang_idx(ii);
        for jj = 1:n_ang
            j = ang_idx(jj);
            if i == j
                J11(ii,jj) = -Q_calc(i) - imag(Ybus(i,i))*Vm(i)^2;
            else
                J11(ii,jj) = Vm(i)*Vm(j)*(real(Ybus(i,j))*sin(Va(i)-Va(j)) ...
                             - imag(Ybus(i,j))*cos(Va(i)-Va(j)));
            end
        end
        for jj = 1:n_mag
            j = mag_idx(jj);
            if i == j
                J12(ii,jj) = P_calc(i)/Vm(i) + real(Ybus(i,i))*Vm(i);
            else
                J12(ii,jj) = Vm(i)*(real(Ybus(i,j))*cos(Va(i)-Va(j)) ...
                             + imag(Ybus(i,j))*sin(Va(i)-Va(j)));
            end
        end
    end

    for ii = 1:n_mag
        i = mag_idx(ii);
        for jj = 1:n_ang
            j = ang_idx(jj);
            if i == j
                J21(ii,jj) = P_calc(i) - real(Ybus(i,i))*Vm(i)^2;
            else
                J21(ii,jj) = -Vm(i)*Vm(j)*(real(Ybus(i,j))*cos(Va(i)-Va(j)) ...
                              + imag(Ybus(i,j))*sin(Va(i)-Va(j)));
            end
        end
        for jj = 1:n_mag
            j = mag_idx(jj);
            if i == j
                J22(ii,jj) = Q_calc(i)/Vm(i) - imag(Ybus(i,i))*Vm(i);
            else
                J22(ii,jj) = Vm(i)*(real(Ybus(i,j))*sin(Va(i)-Va(j)) ...
                             - imag(Ybus(i,j))*cos(Va(i)-Va(j)));
            end
        end
    end

    J  = [J11 J12; J21 J22];
    dx = J \ mismatch;

    % Update voltages
    Va(ang_idx) = Va(ang_idx) + dx(1:n_ang);
    Vm(mag_idx) = Vm(mag_idx) + dx(n_ang+1:end);
    V = Vm .* exp(1j*Va);
end

% Calculate total losses
P_loss = 0;
Q_loss = 0;

for k = 1:size(linedata,1)
    from = linedata(k,1);
    to   = linedata(k,2);
    R    = linedata(k,3);
    X    = linedata(k,4);
    if R == 0 && X == 0, continue; end
    z      = R + 1j*X;
    I_line = (V(from) - V(to)) / z;
    S_loss = I_line * conj(I_line) * z;
    P_loss = P_loss + real(S_loss);
    Q_loss = Q_loss + imag(S_loss);
end

P_loss = P_loss * 100;
Q_loss = Q_loss * 100;

% Display results
fprintf('\n===== LOAD FLOW RESULTS =====\n');
fprintf('Bus | Voltage (pu) | Angle (deg)\n');
fprintf('----+--------------+------------\n');
for i = 1:nb
    fprintf(' %2d |    %.4f    |   %7.4f\n', ...
        i, abs(V(i)), angle(V(i))*180/pi);
end
fprintf('\nTotal Active Loss:   %.4f MW\n', P_loss);
fprintf('Total Reactive Loss: %.4f MVAR\n', Q_loss);