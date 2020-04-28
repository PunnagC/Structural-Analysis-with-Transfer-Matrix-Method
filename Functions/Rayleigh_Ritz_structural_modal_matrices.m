function  matxRR = Rayleigh_Ritz_structural_modal_matrices(panel,lag,bending,torsion,T_type)
%% This finds the modal matricess for the system
% They include:
%              1. modal mass matrix
%              2. modal stiffness matrix
%              3. modal stiffness matrix for axial tension

%---------------------------------INPUT-----------------------------------
% panel   - panelwise discretization properties containing (ATLEAST):
%           dr        - panel/strip width
%           Jp        - polar MOI
%           Acs       - cross sectional area
%           Ggamma    - torsional rigidity (G x gamma)
%           rhoA      - effective mass per unit length
%           rhoJ      - effective rho x Jp
%           EIbending - bending rigidity
%           EIlag     - lead/lag rigidity
% lag     - 'MATLAB' strcuture containing (ATLEAST) panelwise :
%           (1) mode shape(PHI)
%           (2) slope slope (dPHI_dx)
%           (3) 2nd derivative of mode shape (d2PHI_dx2) of lead/lag DOF
% bending - same as above for bending DOF
% torsion - same as above for torsional DOF
% T_type  - whether the structure is rotating (helicopter blades) or not
%           T_type = 0 (if not rotating)
%           T_type = 1 (if rotating)

%---------------------------------OUTPUT-----------------------------------
% matxRR -  Rayleigh Ritz matrices organized in a MATLAB 'structure'

%***********************************NOTE***********************************
%       If some DOF's are missing, it is fine CODE will still work
%***********************************NOTE***********************************

%% Start of actual code

% keyboard
dr        = panel.dr;
Jp        = panel.Jp;
Acs       = panel.Acs;
Ggamma    = panel.Ggamma;
rhoA      = panel.rhoA;
rhoJ      = panel.rhoJ;
EIbending = panel.EIbending;
EIlag     = panel.EIlag;

if T_type == 1 %depends on axial distance
    r = rhoA.*panel.rmid;
else %does NOT depend on axial distance (like axial tension)
    r = 1;
end


%% For bending DOF
if isfield(bending,'PHI') == 1
    d2PHI_dx2 = bending.d2PHI_dx2;
    dPHI_dx   = bending.dPHI_dx;
    PHI       = bending.PHI;
    
    M_bending     = (PHI.*dr)'*(PHI.*rhoA); %Mass matrix contribution from distributed system
    K_bending_EI  = (EIbending.*dr.*d2PHI_dx2)'*d2PHI_dx2;
    K_bending_T   = @(T) (T.*r.*dr.*dPHI_dx)'*dPHI_dx;
    
    matxRR.M_bend    = M_bending;
    matxRR.K_bend_EI = K_bending_EI;
    matxRR.K_bend_T  = K_bending_T;
end

%% For lead-lag DOF
if isfield(lag,'PHI') == 1
    M_lag     = (lag.PHI.*dr)'*(lag.PHI.*rhoA); %Mass matrix contribution from distributed system
    K_lag_EI  = (EIlag.*dr.*lag.d2PHI_dx2)'*lag.d2PHI_dx2;
    K_lag_T   = @(T) (T.*r.*dr.*lag.dPHI_dx)'*lag.dPHI_dx;
    
    matxRR.M_lag     = M_lag;
    matxRR.K_lag_EI  = K_lag_EI;
    matxRR.K_lag_T   = K_lag_T;
end

%% For torsion DOF
if isfield(torsion,'PSI') == 1
    PSI     = torsion.PSI;
    dPSI_dx = torsion.dPSI_dx;
    
    M_torsion        = (rhoJ.*dr.*PSI)'*PSI; %Mass matrix contribution from distributed system
    K_torsion_Ggamma = (Ggamma.*dr.*dPSI_dx)'*dPSI_dx;
    K_torsion_T      = @(T) ((T.*Jp./Acs).*r.*dr.*dPSI_dx)'*dPSI_dx;
    
    matxRR.M_torsion = M_torsion;
    matxRR.K_tor_Gg  = K_torsion_Ggamma;
    matxRR.K_tor_T   = K_torsion_T;
end



end