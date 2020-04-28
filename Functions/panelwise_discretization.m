function [ panel ] = panelwise_discretization(segc,Mp)
%This is to discretize the rotor blade/structure panelwise
ns       = numel(segc.chord); %number of structural segments of the blade

% Initializing vectors of panelwise properties
rhoA      = []; % material density
EIbending = []; % Young's modulus
EIlag     = []; % Young's modulus
Em        = [];
dr        = [];
Ggamma    = []; 
rhoJ      = [];
G        = []; % Rigidity modulus
chord    = []; % chord
rmid     = [];
Ibending = []; % MOI in bending
Ilag     = []; % MOI in lead-lag
EIp      = [];
EA       = [];
Jp       = []; % MOI in bending
Acs      = []; % cross sectional area
gamma    = []; % torsion factor/constant
a        = []; % non-dimensional location of pivot point P, as per airfoil schematic
e        = []; % non-dimensional location of CG point C, as per airfoil schematic 
Et       = []; 
thk_tot  = []; % total segmental thickness

% simple formulas for segmental property evaluations
Icm_function   = @(chord,thk) (1/12)*chord*thk^3;               %m^4, bending MOI about COM
Acs_function   = @(chord,thk) chord*thk;    %m^2, cross sectional area
Jcm_function   = @(chord,thk) (chord*thk/12)*(chord^2 + thk^2) ; %m^4, polar MOI about COM + chord*thk*r_arm^2
gamma_function = @(chord,thk) 0.333*chord*thk^3;                %m^4, torsion constant for rectangle


for i = 1:ns
    sgn     = 1;
    Npaneli = segc.Npanel(i) + 1; % number of panels in ith segment
    vec_1s  = ones(segc.Npanel(i),1); %vector of 1's
   
    r_s        = linspace(0,segc.l(i),Npaneli).*1e0; %panels within a segment
    dr_s       = (r_s(2) - r_s(1));
    rmid_s     = 0.5*dr_s + r_s(1:end - 1);
    
    if i > 1
       rmid_s = rmid_s +  sum(segc.l(1:i - 1))*1e0;
    end
  
    %vectors of panelwise properties
    % material property vectors
    EIbending    = [EIbending; segc.EIbending(i)*vec_1s];
    rhoA         = [rhoA; segc.rhoA(i)*vec_1s];
    EIlag        = [EIlag; segc.EIlag(i)*vec_1s];
    Ggamma       = [Ggamma; segc.Ggamma(i)*vec_1s];
    Acs          = [Acs; segc.Acs(i)*vec_1s];
    rhoJ         = [rhoJ; segc.rhoJ(i)*vec_1s];
    rmid         = [rmid; rmid_s']; %global r
    chord        = [chord; segc.chord(i)*vec_1s];
    EA           = [EA; segc.EA(i)*vec_1s]; %MOI for transverse bending
    EIp          = [EIp; segc.EJp(i)*vec_1s]; %MOI for transverse bending
    Em           = [Em; segc.Em(i)*vec_1s];
    thk_tot      = [thk_tot; segc.thk_s(i)*vec_1s];
    
    eP_T{i}      = segc.eP_T{i}.*vec_1s;
    eS  {i}      = segc.eS{i}.*vec_1s;
    Eshi{i}      = segc.Eshi{i}.*vec_1s;
    
    if i == 2
        sgn = -1;
    end
        
    yu  {i}      = sgn.*fliplr(segc.z{i}(2:end)).*vec_1s;
    yl  {i}      = sgn.*fliplr(segc.z{i}(1:end-1)).*vec_1s;
  
    Et       = [Et; segc.Et(i)*vec_1s];
    Ibending = [Ibending; segc.Em(i)*segc.Ibending(i)*vec_1s]; %MOI for transverse bending
    
    Ilag     = [Ilag; segc.Ilag(i)*vec_1s]; %MOI for lead/lag bending
    Jp       = [Jp; segc.Jp(i)*vec_1s]; %polar area MOI for torsion
    
    gamma    = [gamma; segc.gamma(i)*vec_1s];
    
    a        = [a; segc.a(i)*vec_1s];
    e        = [e; segc.e(i)*vec_1s];
    dr       = [dr, dr_s*vec_1s'];
    
    
% keyboard
end


%% Creating Matlab structure 'panel' for panelwise properties
% All these values are evaluated at the center of each panel
% Dividing the blade into n parts will give n-1 panels

panel.Acs        = Acs;
panel.rhoA       = rhoA;

panel.EIbending  = EIbending;
panel.Em         = Em;
panel.EIlag      = EIlag;
panel.EA         = EA;

panel.rhoJ       = rhoJ;
panel.Jp         = Jp;
panel.Ggamma     = Ggamma;
panel.EIp        = EIp;
panel.Et         = Et;
panel.thk_tot    = thk_tot;
panel.rmid       = rmid';
panel.chord      = chord;
panel.dr         = dr';
panel.b          = 0.5.*panel.chord;
panel.a          = a;
panel.e          = e;
panel.xth        = panel.e - panel.a; %aeroelastic mass imbalance CG to P distance

panel.eP_T       = eP_T;
panel.eS         = eS;
panel.Eshi       = Eshi;
panel.yu         = yu;
panel.yl         = yl;

%% Finding location indix(es) of point mass(es) within rmid vector
nMp = numel(Mp.mass); %number of point masses
for i = 1:nMp
    segc.Mp_panel(i) = find(panel.rmid >= Mp.span(i)*sum(segc.l), 1) - 1;
end


end

