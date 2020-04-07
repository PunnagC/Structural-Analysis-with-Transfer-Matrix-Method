function [Char_Matrix_w,BC_sym,root_cols] = TMM_characteristic_equation(BC,GTM_w,DOF_type,FTM_type)
% Finds the symbolic characteristic equation for TMM formulation 
%----------------------------------INPUT-----------------------------------
% BC       - string, geometric BC on either side ('clamped','pinned','free')
%           example: 'clamped-clamped' OR 'clamped-free'
% GTM_w    - global transfer matrix of system (usually a 4x4 matrix)
% symbolic in 'w'
% DOF_type -  type of elastic DOF either 'bend'ing or 'torsion'
%---------------------------------OUTPUT-----------------------------------
% Char_Matrix_w - characteristic matrix symbolic in 'w'
% BC_sym        - [4x2] matrix of symbolic B/C for root (column 1) and tip
%                 (column 2)
%                 Sequence is: [    PHIx  ; %mode shape
%                               d(PHIx)/dx; %slope
%                                    M    ; %moment
%                                    V    ] %shear force
% root_cols     - Non-zero states at the root
%--------------------------------------------------------------------------

%% Error checking at INPUT
IN_min = 3;  %min inputs
IN_max = 4; %max inputs
narginchk(IN_min,IN_max);
if strcmp(DOF_type,'bend') == 1 && nargin < IN_max
   error('If type = "bend", please enter "FTM_size"');  
end

%% Actual evaluations start from here

idx        = strfind(BC,'-');
BC_tip     = BC(idx + 1:end);
BC_root    = BC(1 : idx-1);
BC_vec{1}  = BC_tip;
BC_vec{2}  = BC_root;

if strcmp(DOF_type,'bend') == 1 %% for Bending DOF
    syms Shi N Phi PHI_x M V
    if FTM_type == 1
        BC_sym = sym(zeros(4,2));
        states_sym = [Phi;PHI_x;M;V];
        for i = 1:2
            BCi = BC_vec{i};
            if strcmp(BCi,'clamped') == 1
                states = [0;0;1;1];
            elseif strcmp(BCi,'free') == 1
                states = [1;1;0;0];
            elseif strcmp(BCi,'pinned') == 1
                states = [0;1;0;1];
            else
                error('Either BC is pinned/clamped/free. Pls input only these')
            end
            BC_sym(:,i) = states.*states_sym;
        end
        tip_rows       = find(BC_sym(:,1) == 0); % Zero states at the tip
        root_cols      = find(BC_sym(:,2) ~= 0); % Non-zero states at the root
        Char_Matrix_w  = GTM_w(tip_rows,root_cols);
    else
        BC_sym = sym(zeros(6,2));
        states_sym = [Shi;N;Phi;PHI_x;M;V];
        for i = 1:2
            BCi = BC_vec{i};
            if strcmp(BCi,'clamped') == 1
                states = [0;1;0;0;1;1];
            elseif strcmp(BCi,'free') == 1
                states = [1;0;1;1;0;0];
            elseif strcmp(BCi,'pinned') == 1
                states = [0;1;0;1;0;1];
            else
                error('Either BC is pinned/clamped/free. Pls input only these')
            end
            BC_sym(:,i) = states.*states_sym;
        end
        tip_rows       = find(BC_sym(:,1) == 0); % Zero states at the tip
        root_cols      = find(BC_sym(:,2) ~= 0); % Non-zero states at the root
        Char_Matrix_w  = GTM_w(tip_rows,root_cols);
    end
%     keyboard
    
elseif strcmp(DOF_type,'torsion') == 1 %% for Torsion DOF
    syms Psi Tau
    BC_sym     = sym(zeros(2,2));
    states_sym = [Psi;Tau];
%     keyboard
    for i = 1:2
        BCi = BC_vec{i};
        if strcmp(BCi,'clamped') == 1
            states = [0;1];
        elseif strcmp(BCi,'free') == 1
            states = [1;0];
        elseif strcmp(BCi,'pinned') == 1
            states = [0;0];
        else
            error('Either BC is pinned/clamped/free. Pls input only these')
        end
        BC_sym(:,i) = states.*states_sym;
    end
    tip_rows       = find(BC_sym(:,1) == 0); % Zero states at the tip
    root_cols      = find(BC_sym(:,2) ~= 0); % Non-zero states at the root
    Char_Matrix_w  = GTM_w(tip_rows,root_cols);
end


end

