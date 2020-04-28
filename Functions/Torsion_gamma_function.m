function gamma = Torsion_gamma_function(width,t,n,optn)
%m^4, torsion constant for rectangle
a = width;
b = t;
% keyboard
if optn == 1 % Using Roark's formula
    % Young, W. C., Budynas, R. G.(2002). Roark's Formulas for Stress and Strain . 7nd Edition, McGraw-Hill, Chapter 10 , pp 401
    
    gamma =  (0.5*a)*(0.5*b)^3*(16/3  - 3.36*(b/a)*(1 - ((b/a)^4)/12));
else
    summation = 0;
    
    for i = 1:n
        summation = summation + (1/(2*i - 1)^5)*tanh(pi*b*0.5*(2*i - 1)/a);
    end
    
    gamma =  (a^3*b/3)*(1 - (192/pi^5)*(a/b)*summation);
end
% keyboard
end