function [FTM_T]  = FTM_torsion_pointmass(w,I)
%CALCULATES THE INDIVIDUAL TRANSFER MATRIXES OF POINT MASSES
    %Elements of FTM for Torsion with axial load 'T'
    %'T' is positive in compression
FTM_T =  [   1  , 0 ;
          -I*w^2, 1];

FTM_T = vpa(FTM_T,8);

end

