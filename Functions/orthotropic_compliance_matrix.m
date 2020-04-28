function [S, Q_bar, Q] = orthotropic_compliance_matrix(Matl,th)
%% Finds orthotropic material properties
%----------------------------------INPUT-----------------------------------
% Mat1 - (compliance matrix)elastic material properties using a MATLAB vector
% th   - fiber oriantation angle in [deg]
%---------------------------------OUTPUT-----------------------------------
% S     - compliance matrix
% Q     - reduced stiffness matrix
% Q_bar - reduced stiffness matrix including fiber rotation effects
%--------------------------------------------------------------------------

%% Assigning the elastic properties to S matrix 
s33 = Matl(1);
s11 = Matl(2);
s12 = Matl(3);
s13 = Matl(4);
s44 = Matl(5);
s55 = Matl(6);
s66 = Matl(7);
E1  = Matl(8);
Y33 = Matl(9);
s22 = s11;
Y22 = E1;


v23 = 1 - s22*0.5/s33;
s23 = -v23/Y22;
s32 = -v23/Y33;

v12   = -s12*Y22;
v13   = -s13*E1;

S = [s11 s12 s13  0   0   0;
     s12 s22 s23  0   0   0;
     s13 s32 s33  0   0   0;
      0   0   0  s44  0   0;
      0   0   0   0  s55  0;
      0   0   0   0   0  s66];

S_reduced = [S(1,1)  S(1,2)           0;
             S(1,2)  S(1,1)           0;
                0       0   2*(S(1,1) - S(1,2))];  
  
C = inv(S); %stiffness matrix

Q = [C(1,1) C(1,2)      0   ;
     C(2,1) C(2,2)      0   ;
        0     0      C(6,6)];

cth = cosd(th);
sth = sind(th);

T = [  cth^2    sth^2     2*sth*cth   ;
       sth^2    cth^2    -2*sth*cth   ;
     -sth*cth  sth*cth  cth^2 - sth^2];
 
T_inv = inv(T);
Q_bar = T_inv*Q*T_inv'; 

end