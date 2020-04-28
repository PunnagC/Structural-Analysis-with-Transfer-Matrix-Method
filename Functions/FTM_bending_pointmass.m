function [FTM]  = FTM_bending_pointmass(w,I,m,theta,type)
%CALCULATES THE INDIVIDUAL FIELD TRANSFER MATRIXES OF POINT MASS
    %Elements of FTM from WickenHeiser Paper
    %Title: "Eigensolution of piezoelectric energy harvesters with
    %geometric discontinuities: Analytical modeling and validation"
%--------------------------------INPUT-------------------------------------
% w     - symbolic variable for circular natural frequency
% I     -  moment of inertia of the point mass
% m     - mass of point mass
% theta - angle between subsequent beam segmentgs
% type  - 1 = 4x4 matrix output, 2 = 6x6 matrix output 
%         if theta == 0 deg then it is safe to get 4x4 else please use
%         6x6 matrix output for accurate results
%---------------------------------OUTPUT-----------------------------------
% FTM  - [4x4] OR [6x6] Field transfer matrix for point mass
%--------------------------------------------------------------------------
narginchk(4,5);
if nargin < 5
    type = 1; %if there are four input args then default output is 4x4 matrix
end

%% Elements of FTM Matrix
cth = cosd(theta);
mw2 = m*w^2;

if type == 1
    if theta ~= 0
        warning('It is not safe to use type == 1, please change to type == 2')
        disp('where; 1 = 4x4 matrix output, 2 = 6x6 FTM matrix output')
    end
    FTM =  [   cth   ,    0  , 0,   0  ;
                0    ,    1  , 0,   0  ;
                0    , -I*w^2, 1,   0  ;
             mw2*cth ,    0  , 0,  cth];
%         -mw2*cth ,    0  , 0,  cth];
    
elseif type == 2
    sth = sind(theta);
    FTM = [   cth,     0 ,    sth  ,      0  ,   0,   0  ;
        mw2*cth,  cth,  mw2*sth,      0  ,   0,  sth ;
        -sth  ,   0 ,    cth  ,      0  ,   0,   0  ;
        0   ,   0 ,      0  ,      1  ,   0,   0  ;
        0   ,   0 ,      0  ,   -I*w^2,   1,   0  ;
        -mw2*sth, -sth, mw2*cth,      0  ,  0,  cth];
else
    error('Type can be either 1 or 2 [1 = 4x4 matrix output, 2 = 6x6 FTM matrix output]')
end
% else
%     if type == 1
%         FTM = eye(4);
%     elseif type ==2
%         FTM = eye(6);
%     end
% end
FTM = vpa(FTM,10);

end

