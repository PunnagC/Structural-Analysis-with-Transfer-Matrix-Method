function [PHI_shape, w_h, strip] = TMM_bending(l_s,type_s,rhoA,EI,Mpoint,T,BC,...
                                  FTM_size,method,r,w_guess111,w_guess,n_roots,w_inc)
                         
%% Evaluates mode shapes and modal frequencies of a segmented beam 
%  structure using TMM formulation

%----------------------------------INPUT-----------------------------------
% l_s      - [1xns] length of each segment
% type_s   - [1xns] element type for each segment (Beam 'B'/Point Mass 'P')
% rhoA     - [1xns] rhoA (mass per unit length) for each segment
% EI       - [1xns] EI (bending rigidity) for each segment
% T        - scalar, applied axial load (+ve = tension, -ve = compression)
% BC       - string, geometric BC on either side ('clamped','pinned','free')
%            example: 'clamped-clamped' OR 'clamped-free'
% FTM_size - whether to perform TMM in 4x4 or more general 6x6 format
% w_guess  - [1x2] guessed min and max root range to search for roots
% r        - cell of size {1xns} where each cell is a vector of size [Npanelx1] 
%           contains numerical global points where the mode shape is evaluated
% n_roots  - scalar, number of modes to evaluate
% w_inc    - used in root finding algorithm, checks for roots in intervals
%            of this value
%---------------------------------OUTPUT-----------------------------------
% PHI_shape - [nsxn_roots] symbolic mode shape expression
% w_h       - [1xn_roots] modal natural frequencies in rad/s
% strip     - structure containing
%             [sum(Npanel)xn_roots] Mode shape
%             [sum(Npanel)xn_roots] Mode shape slope
%             [sum(Npanel)xn_roots] 2nd derivative of mode shape
%             [sum(Npanel)xn_roots] 3rd derivative of mode shape
%             [sum(Npanel)xn_roots] 4th derivative of mode shape
%--------------------------------------------------------------------------

%% Error checking at INPUT
IN_min = 12;  %min inputs
IN_max = 14; %max inputs
narginchk(IN_min,IN_max);
if numel(w_guess) > 1 && nargin < IN_max
   error('If Root Range is provided, # Roots must be provided!');  
elseif numel(w_guess) == 1 && nargin > IN_min
    error((message('If just one Root guess value is provided, remove # Roots and search increment(last 2 arguments)!'))); 
    
end

%% PROBLEM SETUP
syms w x
assume(w > 0)
assume(x >= 0)
% keyboard
nb  = numel(EI);          % number of TMM beam segments
npm = nnz(Mpoint.span); % number of TMM point mass segments
ne  = nb + npm;
% keyboard
if FTM_size == 1
    matrix_size = 4;
elseif FTM_size == 2
    matrix_size = 6;
else
    error('FTM_type can be either 1 or 2 [1 = 4x4 matrix output, 2 = 6x6 matrix output]');
end

%% Finding FTM for every segment
%FTM is a 2D Field Transfer matrix for a segment
%FTM_x_f - 3D matrix where the 3rd dimension is used for segments
%          This matrix is a function of 'x' anf 'f' 
%FTM_f   - 3D matrix like FTM_x_f but only a function of 'f' as segmental
%          lengths are inserted here

BP_idx_global = strfind(char(type_s),'BP');
P_idx_global  = strfind(char(type_s),'P');
B_idx_global  = strfind(char(type_s),'B');

FTM_x_w   = sym(NaN(matrix_size,matrix_size,nb));
FTM_w     = FTM_x_w;
FTM_x_w_B = FTM_x_w;

FTM_x_w_Pm = sym(NaN(matrix_size,matrix_size,npm));

% for beam segments
cb = 0;
cp = 1;
% keyboard
for seg = 1:ne  %sweeping through beam elements 
    if type_s(seg) == 'B'
        cb = cb + 1;
        FTM_x_w_B(:,:,cb)       = FTM_bending(w,x,EI(cb),rhoA(cb),T,FTM_size); %; FTM_Heave(w,x,EI(seg),rhoA(seg))
    elseif type_s(seg) == 'P'
        mass                    = Mpoint.mass(cp); 
        if mass ~= 0
            I                   = Mpoint.I(cp);
            theta               = Mpoint.theta(cp);
            FTM_x_w_Pm(:,:,cp)  = FTM_bending_pointmass(w,I,mass,theta,FTM_size);
            cp                  = cp + 1;
        else
            FTM_x_w_Pm(:,:,cp)  = eye(matrix_size);
        end
    end
end

%% associate/merge FTM(PM=particle mass) with subsequent beam segment

pm_binary = zeros(1,nb); %if PM is present in that beam segment or not
for i = 1:npm
    for j = 1:nb
        check = P_idx_global(i) - B_idx_global(j);
        if check == 1
            beam_idx = j;
            pm_binary(beam_idx) = 1;
        end
    end
end
% pm_binary
cp = 0; %just a counter
for i = 1:nb
    if pm_binary(i) == 1
       cp = cp + 1;
       FTM_x_w(:,:,i) = FTM_x_w_Pm(:,:,cp)*FTM_x_w_B(:,:,i);
    else
        FTM_x_w(:,:,i) = FTM_x_w_B(:,:,i);
    end
    FTM_w(:,:,i) = subs(FTM_x_w(:,:,i),x,l_s(i));
end
FTM_w   = vpa(FTM_w,8);
FTM_x_w = vpa(FTM_x_w,8);  
disp('FTMs evaluated')

%% Finding GTM for entire structure
%GTM_f is a 2D Global Transfer matrix for the entire structure as a
%function of 'f'
%[  SHIx(L)    ;              [ SHIx(0)   ;
%    N(L)      ;                 N(0)     ;
% ---------------------------------------------
%   PHIx(L)    ;                PHIx(0)   ;
% d(PHIx(L))/dx; = [GTM_f]* d(PHIx(0))/dx ;
%    M(L)      ;                 M(0)     ;
%    V(L)]                       V(0) ]

GTM_w = eye(matrix_size);
for seg = 1:nb
    GTM_w = FTM_w(:,:,seg)*GTM_w;    
end

disp('GTM evaluated')  
%% Applying Boundary condition to GTM to get Characteristic equation
% Char_Matrix_w is the matrix whose determinant is the characteristic equation
% Note that Char_Matrix_w is a function of 'w' only.

BC = char(BC);
Char_Matrix_w = TMM_characteristic_equation(BC, GTM_w, 'bend',FTM_size);
Char_Matrix_w = vpa(Char_Matrix_w,5);
Char_Matrix_w = check_characteristic_matrix(Char_Matrix_w,w);
% max_val = max(max(double(subs(Char_Matrix_w,w,1))));
disp('[Bending]Characteristic matrix generated') 

%% Finding Natural frequencies within specified Range
[w_h, nroots] = TMM_root_finding(w_guess, w, Char_Matrix_w,n_roots,w_inc,method,w_guess111);

%% Mode shape 
% Finding the maximum value of a plotted mode shape and then normalizing
% with the max value

for i = 1:nroots
    [PHI_shape_temp] = get_sym_bending_modeshape(w_h(i),w,FTM_x_w,FTM_x_w_B,x,l_s,BC,FTM_size);%FTM_w,(w_a(i),w,FTM_x_w,x,FTM_w,l_s);
    PHI_shape(:,i)   = transpose(simplify(expand(PHI_shape_temp)));
    PHI_shape(:,i)   = vpa(PHI_shape(:,i),8);
end
disp('Finished evaluating symbolic mode shape.') 

%% Mode shape normalization
PHI_value = [];
for i = 1:nb
    PHI_value_s = double(subs(PHI_shape(i,:), x, r{i})); %segmental mode shape
    PHI_value   = vertcat(PHI_value, PHI_value_s);
end

for i = 1:nroots
    max_mode_value(i) = max(abs(PHI_value(:,i)));
    PHI_shape(:,i)    = PHI_shape(:,i)/max_mode_value(i);
end
disp('Finished normalizing symbolic mode shape.') 

%% Discretizing mode shape into strips
PHI              = [];
D_PHI_Dx_strip   = [];
D2_PHI_Dx2_strip = [];
D3_PHI_Dx3_strip = [];
D4_PHI_Dx4_strip = [];

for i = 1:nb
    D_PHI_Dx   = diff(PHI_shape(i,:),x,1);
    D2_PHI_Dx2 = simplify(expand(diff(D_PHI_Dx,x,1)));
    D3_PHI_Dx3 = simplify(expand(diff(D2_PHI_Dx2,x,1)));
    D4_PHI_Dx4 = simplify(expand(diff(D3_PHI_Dx3,x,1)));
    
    PHI              = vertcat(PHI,              double(subs(PHI_shape(i,:),x,r{i})));
    D_PHI_Dx_strip   = vertcat(D_PHI_Dx_strip,   double(subs(D_PHI_Dx,x,r{i})));
    D2_PHI_Dx2_strip = vertcat(D2_PHI_Dx2_strip, double(subs(D2_PHI_Dx2,x,r{i})));
    D3_PHI_Dx3_strip = vertcat(D3_PHI_Dx3_strip, double(subs(D3_PHI_Dx3,x,r{i})));
    D4_PHI_Dx4_strip = vertcat(D4_PHI_Dx4_strip, double(subs(D4_PHI_Dx4,x,r{i})));
    
end

strip.PHI       = PHI;             
strip.dPHI_dx   = D_PHI_Dx_strip;  
strip.d2PHI_dx2 = D2_PHI_Dx2_strip;
strip.d3PHI_dx3 = D3_PHI_Dx3_strip;
strip.d4PHI_dx4 = D4_PHI_Dx4_strip;
disp('Finished generating bending numeric mode shape.')

end


