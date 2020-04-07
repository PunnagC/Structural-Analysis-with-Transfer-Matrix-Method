function [PSI_shape, w_a, strip] = TMM_TorsionV2(l_s,type_s,geom,Ip,GJ,Mpoint,T,J,...
                                              A,BC,method,r,w_guess,n_roots,w_inc)
%----------------------------------INPUT-----------------------------------
% l_s     - [1xns] length of each segment
% type_s  - [1xns] element type for each segment (Beam 'B'/Point Mass 'P')
% rhoA    - [1xns] rhoA (mass per unit length) for each segment
% EI      - [1xns] EI (bending rigidity) for each segment
% T       - scalar, applied axial load (+ve = tension, -ve = compression)
% BC      - string, geometric BC on either side ('clamped','pinned','free')
%           example: 'clamped-clamped' OR 'clamped-free'
% w_guess - [1x2] guessed min and max root range to search for roots
% r       - cell of size {1xns} where each cell is a vector of size [Npanelx1] 
%           contains numerical global points where the mode shape is evaluated
% n_roots -  scalar, number of modes to evaluate
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
IN_min = 12;
IN_max = 15;
narginchk(IN_min,IN_max);
if numel(w_guess) > 1 && nargin < IN_max
   error('If Root Range is provided, # Roots must be provided!');  
elseif numel(w_guess) == 1 && nargin > IN_min
    error('If just one Root guess value is provided, remove # Roots (last argument)!'); 
end

ns = numel(GJ);
% l_s = geom.R;

%% PROBLEM SETUP
syms w x
assume(w > 0);
assume(x >= 0)

nb          = numel(GJ);          % number of TMM beam segments
npm         = nnz(Mpoint.span); % number of TMM point mass segments
ne          = nb + npm;
matrix_size = 2;

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
FTM_w     = sym(NaN(matrix_size,matrix_size,nb));
FTM_x_w_B = FTM_x_w;

% for beam segments
cb = 0;
cp = 1;
% keyboard
for seg = 1:ne  %sweeping through beam elements 
%     keyboard
    if type_s(seg) == 'B'
        cb = cb + 1;
        FTM_x_w_B(:,:,cb)       = FTM_torsion(w,x,GJ(cb),Ip(cb),T,J,A); %; FTM_Heave(w,x,EI(seg),rhoA(seg))
        elseif type_s(seg) == 'P'
        mass                    = Mpoint.mass(cp); 
        if mass ~= 0
            I                   = Mpoint.I(cp);
            theta               = Mpoint.theta(cp);
            FTM_x_w_Pm(:,:,cp)  = FTM_torsion_pointmass(w,I);
            cp                  = cp + 1;
        else
            FTM_x_w_Pm(:,:,cp)  = eye(matrix_size);
        end
    end
end

%% associate/merge FTM(PM) with subsequent beam segment
P_idx_global  = strfind(char(type_s),'P');
B_idx_global  = strfind(char(type_s),'B');
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
%%
% BP_idx = strfind(char(type_s),'BP');
% % keyboard
% cp = 0; %just a counter
% for i = 1:nb
%     PM = find(BP_idx == i,1);
%     if ~isempty(PM) %if point mass is present
%         cp = cp + 1;
%         FTM_x_w(:,:,i) = FTM_x_w_Pm(:,:,cp)*FTM_x_w_B(:,:,i);
%        
%     else
%         FTM_x_w(:,:,i) = FTM_x_w_B(:,:,i);
%     end
%     FTM_w(:,:,i) = subs(FTM_x_w(:,:,i),x,l_s(i));
% end
% FTM_w   = vpa(FTM_w,8);
% FTM_x_w = vpa(FTM_x_w,8);

%% Finding GTM for entire structure
%GTM_f is a 2D Global Transfer matrix for the entire structure as a
%function of 'f'
%[PSIx(L);           [PSIx(0);
% Tau(L) ] = [GTM_f]* Tau(L) ]
GTM_w = eye(2);
for seg = 1:ns
    GTM_w = FTM_w(:,:,seg)*GTM_w;    
end
% GTM_w
%% Applying Boundary condition to GTM to get Characteristic equation
%Char_Matrix_f is the matrix whose determinant is the characteristic equation
% Note that Char_Matrix_f is a function of 'f' only.
BC            = char(BC);
Char_Matrix_w = TMM_characteristic_equation(BC, GTM_w, 'torsion');
disp('[Torsion]Characteristic matrix generated') 


%% Finding Natural frequencies within specified Range
% method = 1;
w_guess111 = [66, 111, 175]; % T = 1*0.0044N
[w_a,nroots] = TMM_root_finding(w_guess, w, Char_Matrix_w,n_roots,w_inc, method, w_guess111);
for i = 1:nroots
    [PSI_shape_temp] = get_sym_torsional_modeshape(w_a(i),w,FTM_x_w,FTM_x_w_B,x,l_s); %(w_eval,w,FTM_x_w,x,l_s)
    PSI_shape(:,i)   = transpose(PSI_shape_temp);
    PSI_shape(:,i)   = vpa(PSI_shape(:,i),5);
end
% keyboard
%% Mode shape normalization
% Finding the maximum value of a plotted mode shape and then normalizing
% with the max value

PSI_value = [];

for i = 1:nb
    PSI_value = vertcat(PSI_value, double(subs(PSI_shape(i,:),x,r{i})));
    
end

for i = 1:nroots
    max_mode_value(i) = max(abs(PSI_value(:,i)));
    PSI_shape(:,i)    = PSI_shape(:,i)/max_mode_value(i);
end
disp('Finished normalizing symbolic mode shape.') 
%% Discretizing mode shape into strips
PSI = [];
D_PSI_Dx_strip = [];
D2_PSI_Dx2_strip = [];
D3_PSI_Dx3_strip = [];

for i = 1:nb
    D_PSI_Dx   = diff(PSI_shape(i,:),x,1);
    D2_PSI_Dx2 = simplify(expand(diff(D_PSI_Dx,x,1)));
    D3_PSI_Dx3 = simplify(expand(diff(D2_PSI_Dx2,x,1)));
    
    PSI              = vertcat(PSI,              double(subs(PSI_shape(i,:),x,r{i})));
    D_PSI_Dx_strip   = vertcat(D_PSI_Dx_strip,   double(subs(D_PSI_Dx,x,r{i})));
    D2_PSI_Dx2_strip = vertcat(D2_PSI_Dx2_strip, double(subs(D2_PSI_Dx2,x,r{i})));
    D3_PSI_Dx3_strip = vertcat(D3_PSI_Dx3_strip, double(subs(D3_PSI_Dx3,x,r{i})));
    
end
strip.PSI       = PSI;%(geom.idx,:);
strip.dPSI_dx   = D_PSI_Dx_strip;%(geom.idx,:);
strip.d2PSI_dx2 = D2_PSI_Dx2_strip;%(geom.idx,:);
strip.d3PSI_dx3 = D3_PSI_Dx3_strip;%(geom.idx,:);
disp('Finished generating torsional numeric mode shape.')

end



