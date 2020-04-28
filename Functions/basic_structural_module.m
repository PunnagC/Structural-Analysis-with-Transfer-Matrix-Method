function [ structural,w_heaveTMM, w_pitchTMM, fHz_imbalance ] = basic_structural_module(segc,panelc,sim,Mpoint, FTM_size,Bend_root_guess,Tor_root_guess,idx_start,idx_end )%,matxRR
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% All parameters are to be evaluated at r = rmode

%% Error checking at INPUT
IN_min = 7;  %min inputs
IN_max = 9; %max inputs
narginchk(IN_min,IN_max);
if nargin <= IN_min
    idx_start = 1;
    idx_end   = numel(panelc.rhoA);
end

geom.type = segc.type;
geom.R    = segc.l;
% geom.idx  = panelc.idx;
geom.np   = segc.Npanel;

l      = segc.l;
type   = segc.type;
rhoA   = segc.rhoA;
rhoJ   = segc.rhoJ;
EIlag  = segc.EIlag;
EIbend = segc.EIbending;
BC     = sim.BC;
rmode  = segc.rmode;
Acs    = segc.Acs;
Jp     = segc.Jp;
T      = 1*sim.T;
% keyboard
%% Calculating mode shapes in lead-lag DOF 
syms x
if sim.lag == 'Y'
    cprintf('*string','***************[Lead/Lag]***************\n'); 
    method = 3;
    [PHI_v_x, w_leadTMM, panel_lag] = TMM_bending(l, type, rhoA, EIlag, Mpoint,...
        T*1, BC, FTM_size, method, rmode,1.1.*[12, 37, 83, 135], sim.root_range_lag ,sim.nlag,5);% 1.1.*[12, 37, 83, 135]
    
    w_leadTMM = reshape(w_leadTMM, [1,numel(w_leadTMM)]);
    strf = sprintf('% 4.4f\t',w_leadTMM*0.5/pi);
    strw = sprintf('% 4.4f\t',w_leadTMM);
    fprintf('Lag natural frequency(Hz), %s\n',strf);
    fprintf('                  (rad/s), %s\n',strw);
    disp(['Segmental Lead/lag Mode Shape, ',char(966),'(x) : ']);
    disp(vpa(PHI_v_x,5));
    structural.lag = panel_lag;
%     if strcmp(sim.plot_mode_shape,'Y') == 1
%         fig1 = mode_shape_plot_panel(1,panel_lag,'Lag',panelc.rmid,sim.BC);
%     end
    structural.PHI_lagx  = PHI_v_x;
else
    panel_lag = 0;
end
%% Calculating mode shapes in transverse Heave/bending
% keyboard
if sim.bending == 'Y'
    cprintf('*string','***************[Bending/Heave]***************\n');
    method = 3;
    [PHI_h_x, w_heaveTMM, panel_bending] = TMM_bending(l, type, rhoA, ...
        EIbend ,Mpoint,T*1,BC,FTM_size, method,rmode, Bend_root_guess,sim.root_range_bending ,sim.nbending,1);
    w_heaveTMM = reshape(w_heaveTMM,[1,numel(w_heaveTMM)]);
    disp('----------------------------------------------------------------');
    strf = sprintf('% 4.4f\t',w_heaveTMM*0.5/pi);
    strw = sprintf('% 4.4f\t',w_heaveTMM);
    fprintf('Heave natural frequency(Hz), %s\n',strf);
    fprintf('                    (rad/s), %s\n',strw);
    disp(['Segmental Heave Mode Shape, ',char(966),'(x) : ']);
    disp(vpa(PHI_h_x,5));
    structural.bending   = panel_bending;
    structural.PHI_bendx = PHI_h_x;
    cprintf('*string','*****************[Completed]*****************\n');
else
    panel_bending = 0;
end
%% Calculating mode shapes in torsion/pitching
if sim.torsion == 'Y'
    cprintf('*string','***************[Torsion/Pitch]***************\n');   
    method = 3;
    [PSI_a_x, w_pitchTMM, panel_torsion] = TMM_Torsion(l,type,geom,rhoJ,...
        segc.Ggamma,Mpoint,T*1,Jp,Acs,BC,method, rmode, sim.root_range_torsion, sim.ntorsion,1);
    % keyboard
    w_pitchTMM = reshape(w_pitchTMM, [1,numel(w_pitchTMM)]);
    disp('----------------------------------------------------------------');
    strf = sprintf('% 4.4f\t',w_pitchTMM*0.5/pi);
    strw = sprintf('% 4.4f\t',w_pitchTMM);
    fprintf('Pitch natural frequency(Hz), %s\n',strf);
    fprintf('                    (rad/s), %s\n',strw);
    disp(['Segmental Pitch Mode Shape, ',char(968),'(x) : ']);
    disp(vpa(PSI_a_x,5));
    structural.torsion = panel_torsion;
    structural.PSI_torx = PSI_a_x;
        
    
    cprintf('*string','*****************[Completed]*****************\n');
else
    panel_torsion = 0;
end
%% Plotting the structural mode shape functions

if strcmp(sim.plot_mode_shape,'Y') == 1
    if sim.lag == 'Y'
        fig1 = mode_shape_plot_panel(1,panel_lag,'Lag',panelc.rmid,sim.BC, w_leadTMM*0.5/pi,T);
    end
    if sim.bending == 'Y'
        fig2 = mode_shape_plot_panel(2,panel_bending,'Heave',panelc.rmid,sim.BC, w_heaveTMM*0.5/pi,T); %fig2 = mode_shape_plot_panel(2,panel_bending,'Heave',panelc.rmid,sim.BC,w_heaveTMM*0.5/pi,T);
    end
    if sim.torsion == 'Y'
        fig3 = mode_shape_plot_panel(3,panel_torsion,'Pitch',panelc.rmid,sim.BC, w_pitchTMM*0.5/pi,T);
    end
end

%% Rayleigh-Ritz modal mass and linear stiffness matrices
matxRR = Rayleigh_Ritz_structural_modal_matrices(panelc,panel_lag,panel_bending,panel_torsion,0);%,idx_start,idx_end
% keyboard

%% Calculating Inertial Coupling MATRIXES
Icoup  = inertial_coupling_matrices(panel_bending,panel_torsion,panel_lag,sim,panelc,Mpoint,segc);

%% structured output
T      = sim.T;

% For bending
M_heave   = matxRR.M_bend;
K_heave   = matxRR.K_bend_EI + 1*matxRR.K_bend_T(T);
fHz_bend           = diag((0.5/pi)*sqrt(M_heave^-0.5*K_heave*M_heave^-0.5));
structural.f_heave = fHz_bend;%w_heaveTMM*0.5/pi;
strf_bend          = sprintf('% 4.4f\t', fHz_bend);
fprintf('Bending natural frequencies   (Hz) : %s \n', strf_bend);

% For torsion
M_torsional = matxRR.M_torsion;
K_torsional = matxRR.K_tor_Gg + 1*matxRR.K_tor_T(T);
fHz_pitch   = diag((0.5/pi)*sqrt(M_torsional^-0.5*K_torsional*M_torsional^-0.5));

structural.f_pitch = fHz_pitch;%w_pitchTMM*0.5/pi;
strf_tor           = sprintf('% 4.4f\t', fHz_pitch);
fprintf('Torsional natural frequencies (Hz) : %s \n', strf_tor);

cprintf('*Strings','*****Structural natural frequencies (Rayleigh Ritz)*****\n');
fprintf('Compare with frequencies obtained using TMM. \n');

M_struct = [    M_heave   , Icoup.Tor_Bend';
            Icoup.Tor_Bend,  M_torsional];
        
K_struct        = blkdiag(K_heave, K_torsional);
fHz_imbalance   = diag((0.5/pi)*sqrt(M_struct^-0.5*K_struct*M_struct^-0.5));

imbalance_check = segc.e ~= 0 + segc.a ~= 0;

if imbalance_check ~= 0 %If there exists mass imbalance (CM axis does not lie on pivot axis)
    strf_imbalance = sprintf('% 4.4f\t', fHz_imbalance');
    fprintf('Natural freq.(with imbalaces) (Hz) : %s \n', strf_imbalance);
end
structural.Icoup   = Icoup;
structural.matxRR  = matxRR; 

%% Structural damping

if strcmp(sim.include_damping,'N') == 1
    C_h(1:sim.nbending,1:sim.nbending) = zeros;
    C_a(1:sim.ntorsion,1:sim.ntorsion) = zeros;
else
    keyboard
    [C_param_h, C_param_a] = PTFE_AR18_Damping_Data(0.5);
%     M_aero_h = M_aero(1:sim.nbending,1:sim.nbending);
    [C_h,alphah,betah] = Damping_struct_exp(M_heave,K_heave,0*1,C_param_h);
%     M_aero_a = M_aero(sim.nbending+1:sim.nbending+sim.ntorsion,sim.nbending+1:sim.nbending+sim.ntorsion);
    [C_a,alphaa,betaa] = Damping_struct_exp(M_torsional,K_torsional,0*1,C_param_a);
end
C_struct = 1.*double(blkdiag(C_h,C_a));
% structural.matxRR.C_struct = C_struct;

alphaM = 0*0.002;
betaK  = 0;
structural.matxRR.C_struct = alphaM.*M_struct + betaK.*K_struct;

end

