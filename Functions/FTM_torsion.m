function [FTM]  = FTM_torsion(w,x,GGamma,rhoJ,T,J,A)
%% CALCULATES: FIELD TRANSFER MATRIX OF SEGMENT IN TORSION
% This is based on the paper cited below:
% Chatterjee, P., & Bryant, M. (2019). 
% 'Analysis of Tension-Tunable Clamped-Clamped Piezoelectric Beams for 
% Harvesting Energy from Wind and Vibration',
% Journal of Intelligent Material Systems and Structures. 
% https://doi.org/10.1177/1045389X19862390

%--------------------------------INPUT-------------------------------------
% w      - symbolic variable for circular natural frequency
% x      - symbolic variable for local beam length function
% GGamma - scalar, effective torsional rigidity for segment
% rhoJ   - scalar, effective MOI/length for segment
% T      - scalar, applied axial load (+ve = tension, -ve = compression)
% A      - cross sectional area of the segment
%---------------------------------OUTPUT-----------------------------------
% FTM  - [2x2] Field transfer matrix for the segment
%--------------------------------------------------------------------------

GGamma = GGamma + T*J/A; %effective torsional rigidity (with Tension effects)
R = w*sqrt(rhoJ/GGamma);
k = R*x; 

% Mode Shape Constants
c0 = cos(k);
c1 = (x/k)*sin(k);

% Elements of FTM Matrix for Torsion
F11 = c0; 
F12 = c1/GGamma;
F21 = -c1*rhoJ*w^2; 
F22 = c0;  
FTM =  [F11, F12;
        F21, F22];

end

