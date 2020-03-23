function [v_TM, v_TE, i_TM, i_TE] = trxline_PPW(k0, er, h, kro)
kz0 = -1j.*sqrt(-(k0.^2-kro.^2));
ks = k0.*sqrt(er);
kzs = -1j.*sqrt(-(ks.^2-kro.^2));
zeta0 = 120*pi;
zeta_s = zeta0/sqrt(er);
Zs_TM = zeta_s.*kzs./ks;
Z0_TM = zeta0.*kz0./k0;
Zs_TE = zeta_s.*ks./kzs;
Z0_TE = zeta0.*k0./kz0;

zup_TM = 0; 
zup_TE = 0; 
zdn_TM = 1j.*Zs_TM.*tan(kzs.*h);
zdn_TE = 1j.*Zs_TE.*tan(kzs.*h);

i_TE = 1./zdn_TE;
i_TM = 1./zdn_TM;
v_TE = 1;
v_TM = 1;
% Gamma = -1;
% Zpar_TM = zup_TM.*zdn_TM./(zup_TM+zdn_TM);
% Zpar_TE = zup_TE.*zdn_TE./(zup_TE+zdn_TE);
% 
% if (z>h)
%     V0_pl_TM = Zpar_TM.*exp(1j.*kz0*h);
%     V0_pl_TE = Zpar_TE.*exp(1j.*kz0*h);
%     v_TM = V0_pl_TM.*exp(-1j.*kz0.*z);
%     v_TE = V0_pl_TE.*exp(-1j.*kz0.*z);
%     i_TM = V0_pl_TM./Z0_TM.*exp(-1j.*kz0.*z);
%     i_TE = V0_pl_TE./Z0_TE.*exp(-1j.*kz0.*z);
% elseif (z>0) && (z<h)
%     Vs_pl_TM = Zpar_TM.*1./(exp(1j.*kzs*h)+Gamma.*exp(-1j.*kzs*h)); 
%     Vs_pl_TE = Zpar_TE.*1./(exp(1j.*kzs*h)+Gamma.*exp(-1j.*kzs*h)); 
%     v_TM = Vs_pl_TM.*(exp(1j.*kzs*z)+Gamma.*exp(-1j.*kzs.*z));
%     v_TE = Vs_pl_TE.*(exp(1j.*kzs*z)+Gamma.*exp(-1j.*kzs.*z));
%     i_TM = Vs_pl_TM./Zs_TM.*(exp(1j.*kzs*z)-Gamma.*exp(-1j.*kzs.*z));
%     i_TE = Vs_pl_TE./Zs_TE.*(exp(1j.*kzs*z)-Gamma.*exp(-1j.*kzs.*z));
% else 
%     v_TM = 0;
%     v_TE = 0;
%     i_TM = 0;
%     i_TE = 0;
% end
end

