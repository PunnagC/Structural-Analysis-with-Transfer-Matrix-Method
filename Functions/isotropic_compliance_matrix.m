function [S, Q] = isotropic_compliance_matrix(Matl)
%% Finds isotropic material properties
%----------------------------------INPUT-----------------------------------
% Mat1 - (compliance matrix)elastic material properties using a MATLAB vector
%---------------------------------OUTPUT-----------------------------------
% S     - compliance matrix
% Q     - reduced stiffness matrix
%--------------------------------------------------------------------------

E = Matl(1);
v   = Matl(2);

G = 0.5*E/(1 + v);

S = [ 1/E -v/E -v/E  0     0     0 ;
     -v/E  1/E -v/E  0     0     0 ;
     -v/E -v/E  1/E  0     0     0 ;
        0     0      0    1/G    0     0 ;
        0     0      0     0    1/G    0 ;
        0     0      0     0     0   1/G];
           
%            keyboard
S_reduced = [S(1,1)  S(1,2)           0;
             S(1,2)  S(1,1)           0;
                0       0   2*(S(1,1) - S(1,2))];
C = inv(S_reduced); %stiffness matrix

Q11 = E/(1 - v^2);
Q12 = abs(v*E/(1 - v^2));
Q66 = G;

Q = [Q11 Q12  0  ;
     Q12 Q11  0  ;
     0    0  Q66];

end