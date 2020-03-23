function [v_TM, v_TE, i_TM, i_TE] = trxline_PEC(k0, kro)
kz0 = -1j.*sqrt(-(k0.^2-kro.^2));
zeta0 = 120*pi;
Z0_TM = zeta0.*kz0./k0;
Z0_TE = zeta0.*k0./kz0;

v_TE = 1;
v_TM = 1;
i_TE = v_TE./Z0_TE;
i_TM = v_TM./Z0_TM;
end

