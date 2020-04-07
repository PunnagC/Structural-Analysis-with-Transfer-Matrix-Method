function [FTM_T]  = FTM_torsion(w,x,GGamma,rhoJ,T,J,A)
%CALCULATES THE INDIVIDUAL FIELD TRANSFER MATRIXES OF SECTIONS
    %Elements of FTM for Torsion with axial load 'T'
    %'T' is positive in compression
%w=2*pi*f;
GGamma = GGamma + T*J/A;
R = w*sqrt(rhoJ/GGamma);
k = R*x;       
%Mode Shape Constants
c0 = cos(k);
c1= (x/k)*sin(k);
%% Elements of FTM Matrix for Torsion
F11 = c0; 
F12 = c1/GGamma;
F21 = -c1*rhoJ*w^2; %-(k*rhoJ/x)*w^2*sin(k)%
F22 = c0;  
FTM_T =  [F11, F12;
          F21, F22];

end

