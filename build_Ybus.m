% Step 2: Build Y-bus admittance matrix

nb = size(busdata, 1);   % 30 buses
nl = size(linedata, 1);  % 40 lines

Ybus = zeros(nb, nb);    % start with empty matrix

for k = 1:nl
    from = linedata(k,1);
    to   = linedata(k,2);
    R    = linedata(k,3);
    X    = linedata(k,4);
    B    = linedata(k,5);

    % Line admittance
    if R == 0
        y = 1/(1j*X);    % pure transformer
    else
        y = 1/(R + 1j*X);
    end

    % Shunt charging
    yc = 1j*B/2;

    % Fill Y-bus
    Ybus(from,from) = Ybus(from,from) + y + yc;
    Ybus(to,to)     = Ybus(to,to)     + y + yc;
    Ybus(from,to)   = Ybus(from,to)   - y;
    Ybus(to,from)   = Ybus(to,from)   - y;
end

disp('Y-bus built successfully');
disp(['Y-bus size: ' num2str(size(Ybus,1)) 'x' num2str(size(Ybus,2))]);