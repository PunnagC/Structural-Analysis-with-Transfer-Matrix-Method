clc 
clearvars
close all

addpath('.\Functions\')
addpath('.\Materials\')

%% Geometrical Parameters
L             = 236e-3;
ts            = 0.26e-3;
AR            = 18.1538;
PZT_thk       = 0.1e-3;
PZT_coverage  = 0.17;
e             = 0.0; 
a             = 0.0;

%% Materials used
% [Mat1] = read_material_data('Powerfilm_solar1.matl');
% [Mat2] = read_material_data('QP16N_isotropic_ANSYS.matl');
[Mat1] = read_material_data('PTFE_ribbonFRF.matl');
[Mat2] = read_material_data('PTFE_ribbonFRF.matl');
Mats   = {Mat1, Mat2};

[seg, Mpoint] = geometry_layup(Mats, L,ts,AR,PZT_thk,PZT_coverage,e,a); %Mats,L,ts,AR,PZT_thk,PZT_coverage,e,a
segc          = segmental_discretization(seg, Mats);     %(Mats,L,ts,AR,PZT_thk_ratio,PZT_coverage,e)
[panelc]      = panelwise_discretization(segc,Mpoint);
Pcr           = pi^2*segc.EIbending(2)/(0.25*sum(segc.l(2))^2); %critical buckling load

%% Over riding some simulation parameters
Tfactor      = 0; %applied axial tension factor(non-dimensional)
sim          = get_structural_simulation_options; 
sim.nbending = 3; %No. of bending/heave modes to include
sim.ntorsion = 2; %No. of torsion/pitch modes to include
sim.theta0   = 0;
sim.T        = 0;%Pcr*Tfactor;

%% Structural MODULE
FTM_size  = 1;
idx_start = 1;
Bend_root_guess = 2*pi.*[8.9954   25.4648   49.4972   81.1690].*0.25;
Tor_root_guess  = Bend_root_guess;
idx_end         = numel(panelc.rhoA);
[structural]    = basic_structural_module(segc,panelc,sim,Mpoint,FTM_size,Bend_root_guess,Tor_root_guess,idx_start,idx_end);  %

