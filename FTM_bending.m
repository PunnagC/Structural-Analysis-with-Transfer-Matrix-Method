function [FTM]  = FTM_bending(w,x,EI,rhoA,T,type)
%% CALCULATES THE INDIVIDUAL FIELD TRANSFER MATRIXES OF SECTIONS
    %Elements of FTM from WickenHeiser Paper
    %Title: "Eigensolution of piezoelectric energy harvesters with
    %geometric discontinuities: Analytical modeling and validation"
%--------------------------------INPUT-------------------------------------
% w    - symbolic variable for circular natural frequency
% x    - symbolic variable for local beam length function
% EI   - scalar, effective bending rigidity value for segment
% rhoA - scalar, effective mass/length value for segment
% T    - scalar, applied axial load (+ve = tension, -ve = compression)
%---------------------------------OUTPUT-----------------------------------
% FTM  - [4x4] OR [6x6] Field transfer matrix for the beam segment
%--------------------------------------------------------------------------

narginchk(5,6);
if nargin < 6
    type = 1; %if there are four input args then default output is 4x4 matrix
end

Bet = (rhoA*w^2/EI)^0.25; %Beta for section
gam = (T/(2*EI))^0.5;     %Gamma for section
% assume(Bet,'real');
% assume(gam,'real');
M = (gam^2 + sqrt(gam^4 + Bet^4))^0.5;
N = (-gam^2 + sqrt(gam^4 + Bet^4))^0.5;

%% Mode Shape Constants
coshMx = cosh(M*x);
sinhMx = sinh(M*x);
cosNx  = cos(N*x);
sinNx  = sin(N*x);

c0     = (N^2*coshMx + M^2*cosNx)/(N^2 + M^2);
c1     = ((N^2/(N^2*M + M^3))*sinhMx + (M^2/(N^3 + N*M^2))*sinNx)/x;
c2     = (1/(x^2*(N^2 + M^2)))*(coshMx - cosNx);
c3     = ((1/(N^2*M + M^3))*sinhMx - (1/(N^3 + N*M^2))*sinNx) /x^3;

% Betx = Bet*x;
% c0 =      0.5      *(cosh(Betx) + cos(Betx));
% c1 = (0.5/(Betx))  *(sinh(Betx) + sin(Betx));
% c2 = (0.5/(Betx)^2)*(cosh(Betx) - cos(Betx));
% c3 = (0.5/(Betx)^3)*(sinh(Betx) - sin(Betx));
% 
c0 = vpa(c0,5);
c1 = vpa(c1,5);
c2 = vpa(c2,5);
c3 = vpa(c3,5);

%% checking for high values
% M_val = double(subs(M,[w,x],[100,dl]));
% N_val = double(subs(N,[w,x],[100,dl]));
% c0_val =(N_val^2*cosh(M_val*dl) + M_val^2*cos(N_val*dl))/(N_val^2 + M_val^2);
% c0_val = double(c0_val);
% c1_val = ((N_val^2/(N_val^2*M_val + M_val^3))*sinh(M_val*dl) + (M_val^2/(N_val^3 + N_val*M_val^2))*sin(N_val*dl))/dl;
% c1_val = double(c1_val);
% c2_val =(1/(dl^2*(N_val^2 + M_val^2)))*(cosh(M_val*dl) - cos(N_val*dl));
% c2_val = double(c2_val);
% c3_val =((1/(N_val^2*M_val + M_val^3))*sinh(M_val*dl) - (1/(N_val^3 + N_val*M_val^2))*sin(N_val*dl)) /dl^3;
% c3_val = double(c3_val);
% c_val = [c0_val, c1_val, c2_val, c3_val];
% MN_val = [M_val, N_val];
%% Normalizing values 
% c0 = 1;
% c1 = c1/c0;
% c2 = c2/c0;
% c3 = c3/c0;
%% Elements of FTM Matrix
F33 = c0;
F34 = (T*c3*x^3)/EI + c1*x;
F35 = (c2*x^2)/EI;
F36 = -(c3*x^3)/EI;

F43 = (c3*x^3*rhoA*w^2)/EI;
F44 = (T*c2*x^2)/EI + c0 ;
F45 = (c1*x)/EI + (T*c3*x^3)/EI^2;
F46 = -(c2*x^2)/EI;

F53 = c2*x^2*rhoA*w^2;
F54 = (c3*T^2*x^3)/EI + c1*T*x + c3*rhoA*x^3*w^2;
F55 = (T*c2*x^2)/EI + c0;
F56 = - c1*x - (T*c3*x^3)/EI;

F63 = -c1*x*rhoA*w^2;
F64 = -c2*x^2*rhoA*w^2;
F65 = -(c3*x^3*rhoA*w^2)/EI;
F66 = c0;

FTM_temp =  [F33, F34, F35, F36;
             F43, F44, F45, F46;
             F53, F54, F55, F56;
             F63, F64, F65, F66];
         
if type == 1
   FTM = FTM_temp;
elseif type == 2
   FTM1      = [     1    , 0;
                x*rhoA*w^2, 1];
   FTM       = blkdiag(FTM1, FTM_temp);
else
    error('Type can be either 1 or 2 [1 = 4x4 matrix output, 2 = 6x6 matrix output]');
end
% keyboard    
FTM = vpa(FTM,10);
end

