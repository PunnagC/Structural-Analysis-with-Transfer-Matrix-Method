function [Icoup] = inertial_coupling_matrices(panelh,panela,panel_lag,sim,panel,Mpoint,seg)
% This evaluates all the terms other than mass terms which has
% state_dot_dot coefficients
% Bending-Torsion, Lag-Torsion couplings are examples of such
nMp = nnz(Mpoint.mass); %number of point masses

Torsion_Heave_Coup_Mp     = zeros(sim.ntorsion,sim.nbending);
Torsion_Lag_Coup_Mp       = zeros(sim.ntorsion,sim.nlag);
Heave_Coup_Mp_dampinglike = zeros(sim.ntorsion,sim.nbending);
Lag_Coup_Mp_dampinglike   = zeros(sim.ntorsion,sim.nlag);
Heave_modal_Mp            = zeros(1,sim.nbending);
MOI_modal_Mp              = zeros(1,sim.ntorsion);
% keyboard
if nMp ~= 0
    for i = 1:nMp
        
        Mp_panel = seg.Mp_panel(i);
        b        = panel.b(Mp_panel);
        a        = panel.a(Mp_panel);
        e        = panel.e(Mp_panel);
        rhoA     = panel.rhoA(Mp_panel);
        dr       = panel.dr(Mp_panel);
        p        = Mpoint.p(i);
        Mp       = Mpoint.mass(i);%geom.c(2) %kg point mass
        m        = rhoA*dr; %mass of airfoil strip
        cg_dash  = (Mp*(1 + p)*b + m*(1 + e)*b)/(m + Mp);
        e_dash   = cg_dash/b - 1;
        xth_dash = e_dash - a;
        
        PSI         = panela.PSI(Mp_panel,:);
        PHI_bending = panelh.PHI(Mp_panel,:);
        mbxth       = Mp*b*xth_dash;
        
        % for bending-torsion DOF (+ve sign when evaluated at LHS)
        Torsion_Heave_Coup_Mp = Torsion_Heave_Coup_Mp + mbxth*PSI'*PHI_bending;
        
        % Inertial mass due to point mass - useful for base excitations 
        Heave_modal_Mp        = Heave_modal_Mp + Mp*PHI_bending(Mp_panel,:);
        MOI_modal_Mp          = MOI_modal_Mp + mbxth*PSI(Mp_panel,:);
        
        % Bending DOF coupled - damping like terms (+ve sign when evaluated at RHS)
        Heave_Coup_Mp_dampinglike = Heave_Coup_Mp_dampinglike + mbxth*PSI.^3'*PHI_bending; %evaluated at RHS of EOM
        
        if sim.nlag ~=0 || isempty(sim.nlag) == 0 % for lag-torsion DOF
            PHI_lag = panel_lag.PHI(Mp_panel,:) ;
            % Lag-Torsional Inertial coupling (+ve sign when evaluated at LHS)
            Torsion_Lag_Coup_Mp = Torsion_Lag_Coup_Mp + mbxth*PSI.^2'*PHI_lag;
            
            % Lag DOF coupled - damping like terms  (-ve sign when evaluated at RHS)
            Lag_Coup_Mp_dampinglike = - mbxth*PSI.^2'*PHI_lag + Lag_Coup_Mp_dampinglike ;  %evaluated at RHS of EOM
        end
        
        panel.xth(Mp_panel) = xth_dash;
        
    end
else
    Heave_Coup_Mp_dampinglike = 0;
    Torsion_Heave_Coup_Mp     = 0;
    
    Torsion_Lag_Coup_Mp       = 0;
    Lag_Coup_Mp_dampinglike   = 0;
end
% keyboard

%% For distributed system - shortening term names to make compact equations

rhoA        = panel.rhoA;
b           = panel.b;
xth         = panel.xth;
dr          = panel.dr;

PSI = panela.PSI;

PHI_bending = panelh.PHI;
% keyboard
% if isempty(panel_lag) == 0 || panel_lag ~= 0
%     PHI_lag = panel_lag.PHI;%panel_lag.PHI;
% else
    [np,~]  = size(PSI);
    PHI_lag = zeros(np,1);%panel_lag.PHI;
% end

%% Modal inertial mass calculations - useful for base excitations
% keyboard
Heave_modal_distr_mass = (rhoA.*dr)'*PHI_bending; %generalized force
Heave_modal_mass       = -(Heave_modal_Mp + Heave_modal_distr_mass);%Note the negative sign at front
% keyboard
MOI_modal_distr        = (rhoA.*b.*xth.*dr)'*PSI;
% if MOI_modal_Mp ~= []
    MOI_modal          = -(MOI_modal_distr + MOI_modal_Mp);
% else
%     MOI_modal              = -(MOI_modal_distr);
% end

%% Terms for distributed mass imbalance couplings
% Bending-Torsional Inertial coupling (+ve sign when evaluated at LHS)
Torsion_Heave_Coup_distributed = (rhoA.*b.*xth.*dr.*PSI)'*PHI_bending; % int(rhoA*b*xth*[PSI]*[PHI]*dx)


% Lag-Torsional Inertial coupling (+ve sign when evaluated at LHS)
Torsion_Lag_Coup_distributed   = (rhoA.*b.*xth.*dr.*PSI.^2)'*PHI_lag; % this needs to be multiples by .*{ra}'

Torsion_Heave_Coup  = (Torsion_Heave_Coup_distributed + Torsion_Heave_Coup_Mp); %Check the -ve sign at the front when h is measured positive from origin

Torsion_Lag_Coup    = Torsion_Lag_Coup_distributed   + Torsion_Lag_Coup_Mp;

%% Terms for damping like couplings
% Bending DOF coupled - damping like terms  (+ve sign when evaluated at RHS)
Heave_Coup_distributed_dampinglike = (rhoA.*b.*xth.*dr.*PSI.^3)'*PHI_bending;

% Lag DOF coupled - damping like terms (-ve sign when evaluated at RHS)
Lag_Coup_distributed_dampinglike   = -(rhoA.*b.*xth.*dr.*PSI.^2)'*PHI_lag;

Heave_Coup_dampinglike = Heave_Coup_distributed_dampinglike + Heave_Coup_Mp_dampinglike;

Lag_Coup_dampinglike   = Lag_Coup_distributed_dampinglike   + Lag_Coup_Mp_dampinglike;

%% Displaying the matrices
disp('----------------------------------------------------------------');
disp('Inertial Torsion-Bending coupling Matrix, [M_coup] : ');
(double(Torsion_Heave_Coup))
disp('----------------------------------------------------------------');
disp('Inertial Torsion-Lag coupling Matrix, [M_coup] : ');
(double(Torsion_Lag_Coup))

Icoup.Tor_Bend         = Torsion_Heave_Coup;
Icoup.Tor_lag          = Torsion_Lag_Coup;
Icoup.Bend_Clike       = Heave_Coup_dampinglike;
Icoup.Lag_Clike        = Lag_Coup_dampinglike;
Icoup.Heave_modal_mass = Heave_modal_mass;
Icoup.MOI_modal        = MOI_modal;
% keyboard
end
%% Inertial Force for bending DOF

% Heave_modal_Mp         = Mp*panelh.PHI(Mp_panelNo,:); %generalized force
% Heave_modal_distr_mass = (panelh.rhoA_vec.*panelh.dl)'*panelh.PHI; %generalized force
% Heave_modal_mass       = -(Heave_modal_Mp + Heave_modal_distr_mass);%Note the negative sign at front
% 
% %% Inertial Moment for torsional DOF -  small angles only !
% % valid for small angle of theta such that cos(theta) = 1
% % For large angle assumption, it has to be evaluated within ode
% 
% MOI_modal_distr = (panelh.rhoA_vec.*panelh.b.*panelh.xth.*panelh.dl)'*panela.PSI;
% MOI_modal_Mp    = Mp*panelh.b(Mp_panelNo)*panelh.xth(Mp_panelNo)*panela.PSI(Mp_panelNo,:); %generalized force
% MOI_modal       = -(MOI_modal_distr + MOI_modal_Mp);

