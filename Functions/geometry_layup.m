function [ seg, Mpoint] = geometry_layup_V3(Mats,L,ts,AR,PZT_thk,PZT_coverage,e,a)
%%%%%% FIxed PZT thickness
%% Use this as data input file for aero-elastic analysis
%% Aerodynamic properties/data
% keyboard
PZT_l     = PZT_coverage*L;
%% Structural segmental properties/data
seg.theta0_deg = [0, 0, 0]; %static segmental pre-twist

% geometrical properties
seg.chord = (L/AR).*[1,   1,   1]; %m, chord length
seg.l     = [PZT_l, (L - 2*PZT_l), PZT_l]; %m, span/length
thk{1}    = ts.*    [1, 1, 1]; %m, thickness, assuming flat plate
thk{2}    = [1, 1, 1].*PZT_thk; %m, thickness, assuming flat plate CHANGE THIS [[[[{{{{{{to e-3}}}}}}]]]]


seg.layer_arrangement = [2 0 2 ;
                         1 1 1 ;
                         2 0 2];
[ar, ac] = size(seg.layer_arrangement); %arrangements row (ar) and arrangement column (ac)
seg.thk  = zeros(ar, ac);
% mat_loc  = seg.layer_arrangement;                     
% mat_loc(mat_loc > 0) = 1;
% keyboard
for r = 1:ar
    for c = 1:ac
        MatNo = seg.layer_arrangement(r,c);
        if MatNo~= 0
           seg.thk(r,c) = thk{MatNo}(c);
        end
    end
end
        
% keyboard                     



% airfoil properties
seg.a      = [a, a, a]; %non-dimensional location of pivot point P, as per airfoil schematic
seg.e      = [e, e, e]; %non-dimensional location of CG point C, as per airfoil schematic
seg.Npanel = round(seg.l.*1e3); %segmental structural panels
seg.type   = ['BBB']; % B- beam type, P - point mass type
% seg.type   = ['B','B','B']; % B- beam type, P - point mass type
% point/lumped mass location along span
Mpoint.span  = [0, 0];         % non-dim location along span
Mpoint.p     = [0, 0];         % non-dim location along chord
Mpoint.mass  = [0, 0].*1e-3; % Kg, lumped mass
Mpoint.I     = [0, 0];              % MOI
Mpoint.theta = [0, 0];

%% checking for segmental consistency

nUNQ = unique(seg.layer_arrangement);

nseg = [numel(seg.theta0_deg),ac,...
numel(seg.chord),numel(seg.l),numel(seg.a),numel(seg.e),numel(seg.Npanel)];

nunique = unique(nseg);
if numel(nunique) > 1
    error('Please check you have entered all segmental properties !')
end

if nnz(nUNQ) > numel(Mats)
   error('Total number of unique materials used in "Arrangement" exceeds the total materials supplied')
   disp('Please pass in more materials')
end
    

%% checking for Mp consistency
% nMp = [numel(Mpoint.span),numel(Mpoint.p),numel(Mpoint.mass),numel(Mpoint.I)];
nMp = [nnz(Mpoint.span),nnz(Mpoint.p),nnz(Mpoint.mass),nnz(Mpoint.I)];

nunique = unique(nMp);
if numel(nunique) > 1
    error('Please check you have entered all Point mass (Mp) properties !')
end

%% Checking combination type

strlen = numel(seg.type);
if strlen ~= nnz(Mpoint.span) + ar %if strlen ~= numel(Mpoint.span) + numel(seg.E) 
    error('seg.type ~= total(beam segs + point masses)')
end

end

