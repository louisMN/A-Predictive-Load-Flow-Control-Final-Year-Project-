function cost = mpc_cost_fn(u, LF_pred, P_nom, Q_nom, Ybus_c, ...
                             busdata, linedata_c, s_bus, r1_bus, r2_bus, ...
                             y_se, w_v, w_p, w_u, u_prev, Np)
% MPC cost function — evaluates cost over prediction horizon
cost = 0;
for k = 1:Np
    lf = LF_pred(k);
    [Ploss, Vdev] = run_gupfc_lf(u, lf, P_nom, Q_nom, Ybus_c, ...
                                  busdata, linedata_c, s_bus, r1_bus, ...
                                  r2_bus, y_se);
    du   = u - u_prev;
    cost = cost + w_v*Vdev^2 + w_p*Ploss^2 + w_u*norm(du)^2;
end
end
