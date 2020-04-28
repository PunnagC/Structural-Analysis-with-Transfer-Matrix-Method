function sim = get_structural_simulation_options

% structural modes to be considered
sim.nbending = 4;
sim.ntorsion = 3;
sim.nlag     = 0;

sim.U_speed    = linspace(1,10,10);%3,12,19

% geometric boundary condition
sim.BC                 = 'clamped-clamped';%'clamped-clamped' 'pinned-pinned' 'clamped-free'
sim.root_range_lag     = [1,500];
sim.root_range_bending = [1,1400];
sim.root_range_torsion = [1,800];

sim.bending = 'Y';
sim.torsion = 'Y';
sim.lag     = 'N';

sim.include_damping = 'N';
sim.plot_mode_shape = 'Y';

sim.Pcr             = 0.0044;%0.0467
sim.T               = 0*sim.Pcr;%sim.Pcr*1;

sim.v  = 0.5; %non-dim, chordwise location of interest 
sim.lv = 0.5; %non-dim, spanwise location of interest

if sim.v > 1 || sim.v < -1
   error(' -1 <= sim.v <= 1') 
end
if sim.lv > 1 || sim.v <= 0
   error(' 0 < sim.lv <= 1') 
end
end

