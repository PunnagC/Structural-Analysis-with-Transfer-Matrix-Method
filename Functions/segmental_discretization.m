function segc = segmental_discretization(segc, Mats)
%% This function geometric and material using segmental information
% to evaluate effective segmental properties
% keyboard
[~, ns]     = size(segc.layer_arrangement); %number of structural segments of the blade
e0          = 8.854e-12; %F/m, permittivitty of free space       

Icm_function = @(w,thk) (1/12)*w*thk^3;                %m^4, bending MOI about COM
Acs_function = @(w,thk) w*thk;                         %m^2, cross sectional area
Jcm_function = @(w,thk) (w*thk/12)*(w^2 + thk^2) ; %m^4, polar MOI about COM + chord*thk*r_arm^2

for i = 1:ns %sweeping through segments
    nlunique   = unique(segc.layer_arrangement(:,i));
    Npaneli    = segc.Npanel(i) + 1;
    r_s        = linspace(0,segc.l(i),Npaneli).*1e0; %panels within a segment
    dr_s       = (r_s(2) - r_s(1));
    rmid_s     = 0.5*dr_s + r_s(1:end - 1);
    rmode      = rmid_s;
    width      = segc.chord(i);

    if nnz(nlunique) > 1
        out   = composite_layer_seg_prop(segc.layer_arrangement(:,i),width, segc.thk(:,i), Mats, Icm_function,Acs_function,Jcm_function,e0 ); 
    else
       idx   = find(segc.layer_arrangement(:,i)~=0); %Finding the index number for the material which is non-zero
       mat_i = Mats{segc.layer_arrangement(idx,i)};
       
       if strcmp(mat_i.type, 'Isotropic') == 1
          [S, ~] = isotropic_compliance_matrix(mat_i.M);
       end
%        keyboard
       thk  = segc.thk(idx,i);
       cE   = mat_i.M(1);
       G    = 1/S(6,6);
       z    = [0.5*thk, -0.5*thk];
       d31  = mat_i.d31;
       eP_T = d31*cE;
       
       if i == 1
          Eshi = -1/thk; 
       else
          Eshi = 1/thk; 
       end
       Eshi(~isfinite(Eshi)) = 0;
       
       out     = single_layer_seg_prop(width, thk, cE, G, mat_i.rho,  Icm_function,Acs_function,Jcm_function); 

       out{15} = eP_T;
       out{16} = e0*mat_i.er;
       out{17} = Eshi;
       out{18} = z;
       out{19} = thk;
    end
    
    
    segc.Ibending(i)  = double(out{1});
    segc.Ilag(i)      = double(out{2});
    segc.Jp(i)        = double(out{3});
    segc.Acs(i)       = double(out{4});
    segc.gamma(i)     = double(out{5});
    segc.rhoA(i)      = double(out{6});
    segc.EIbending(i) = double(out{7});
    segc.EIlag(i)     = double(out{8});
    segc.rhoJ(i)      = double(out{9});
    segc.Ggamma(i)    = double(out{10});
    segc.EA(i)        = double(out{11}); 
    segc.EJp(i)       = double(out{12});
    segc.Em(i)        = double(out{13});
    segc.Et(i)        = double(out{14});
    
    segc.eP_T{i}      = double(out{15});
    segc.eS{i}        = double(out{16});
    segc.Eshi{i}      = double(out{17});
    segc.z{i}         = double(out{18});
    
    segc.rmode{i}     = double(rmode');
    segc.thk_s(i)     = double(out{19});
end

end

function [out] = single_layer_seg_prop(chord, thk, E, G, rho,Icm_function,Acs_function,Jcm_function)
% to evaluate segmental properties if the segment is of a single
% layer
Ibending_s = Icm_function(chord,thk);
Ilag_s     = Icm_function(thk,chord);
Jp_s       = Jcm_function(chord,thk) ;
Area_cs_s  = Acs_function(chord,thk);
EA_s       = E.*Area_cs_s;
gamma_s    = Torsion_gamma_function(chord,thk,100,1);
EIp_s      = E.*Jp_s;%seg.E(i).*gamma_s
Et_s       = E.*thk;

         
out = {Ibending_s, Ilag_s, Jp_s, Area_cs_s, gamma_s, rho.*Area_cs_s...
       E.*Ibending_s,E.*Ilag_s, rho.*Jp_s,G.*gamma_s,EA_s,EIp_s, E, Et_s} ;

end




function [out] = composite_layer_seg_prop(layer_arrangement,chord, thk_layers, Mat,Icm_function,Acs_function,Jcm_function,e0)
% to evaluate "effective segmental properties" if the segment is of
% multiple layers
nl        = numel(thk_layers);

rhoA_s    = 0;
Acs       = 0;
thk       = 0;
z         = [];
thk_tot   = sum(thk_layers); %total thickness of the segment
mid_layer = ceil(nl/2); %
z(1)      = -thk_tot/2;
c         = 0;
rhoJ_s    = 0;
for j = 1:nl
    if layer_arrangement(j) ~= 0
        c = c + 1;
        mat_no   = layer_arrangement(j); %material number
        cE       = Mat{mat_no}.cE;
        d31(1,j) = Mat{mat_no}.d31;
        if strcmp(Mat{mat_no}.type, 'Orthotropic') == 1
            
            [~, Q_bar{c},~] = orthotropic_compliance_matrix(Mat{mat_no}.M,0);
            
        elseif strcmp(Mat{mat_no}.type, 'Isotropic') == 1
            
            [~, Q_bar{c}] = isotropic_compliance_matrix(Mat{mat_no}.M);
            
        end
        thk_j  = thk_layers(j); %thickness of jth layer
        Eshi_s = 1/thk_j;        Eshi_s(~isfinite(Eshi_s)) = 0;
        
        if j < mid_layer % layer lies above mid layer
            z(c + 1) = z(c) + thk_layers(j);
            r        = z(c + 1) - 0.5*thk_layers(j);
            Eshi_multiplier = -1;
            
            
        elseif j == mid_layer
            z(c + 1) = z(c) + thk_layers(j)/2;
            r        = 0;
            Eshi_multiplier = 1;
        else % layer lies below mid layer
            r        = -(z(c - mid_layer) + 0.5*thk_layers(j));
            Eshi_multiplier = 1;
        end
        Eshi(1,j) = Eshi_multiplier*Eshi_s;
        Area_cs_s   = Acs_function(chord,thk_j);
        
        Acs    = Acs + Area_cs_s;
        rhoA_s = rhoA_s + Mat{mat_no}.rho*Area_cs_s;
        thk    = thk + thk_j;
        
        J      = Jcm_function(chord,thk_j) + Acs_function(chord,thk_j)*r^2;
        rhoJ_j = Mat{mat_no}.rho*J;
        
        rhoJ_s     = rhoJ_s + rhoJ_j;
        eP_T(1,j)  = d31(1,j)*cE;
        eS(1,j)    = Mat{mat_no}.er*e0;
        
    end
end
    
z = z(1:end - 1);
z = [z, -fliplr(z)];

[A, B, D]            = laminate_matrices(Q_bar,z);
[E1m, E1b, E2b, G_s] = laminate_properties(A, B, D, thk_tot);



Ilag_s      = Icm_function(thk_tot,chord);
Ibending_s  = Icm_function(chord,thk_tot);
J_s         = Jcm_function(chord,thk_tot);
EA_s        = E1m*Acs;
EIbending_s = E1b*Ibending_s;
EIlag_s     = E2b*Icm_function(thk_tot,chord);
Et_s        = E1m.*thk_tot;
gamma_s     = Torsion_gamma_function(chord,thk_tot,100,1);
Ggamma_s    = G_s*gamma_s;

rho_eff = rhoA_s/Acs;
rhoJ_s  = rho_eff*J_s;

out = {Ibending_s, Ilag_s, J_s, Acs, gamma_s, rhoA_s, EIbending_s,...
       EIlag_s,rhoJ_s, Ggamma_s, EA_s, E1m.*J_s, E1m, Et_s, eP_T,eS,Eshi,z,thk};


end

function [E1m, E1b, E2b, G_s] = laminate_properties(A, B, D, thk_tot)
% Composite laminate effective modululii based on CLI theory
A_inv = inv(A);
delta = inv(D - B*A_inv*B);
alpha = A_inv + A_inv*B*delta*B*A_inv;

E1m =  1/(thk_tot  *alpha(1,1)); %membrane modulus along 1 direction
E1b = 12/(thk_tot^3*delta(1,1));
E2b = 12/(thk_tot^3*delta(2,2));
G_s = 12/(thk_tot^3*delta(3,3));

end