function [A, B, D] = laminate_matrices(Q_bar,z)
% Q
% Source: http://nptel.ac.in/courses/112104168/L26.pdf
nQ = numel(Q_bar);
% keyboard
for i = 1:3
    for j = 1:3
        A_temp = 0;
        B_temp = 0;
        D_temp = 0;
         for k = 2:numel(z)
%              keyboard
             A_temp = A_temp + ( z(k)   - z(k-1)  )*Q_bar{k-1}(i,j);
             B_temp = B_temp + ((z(k)^2 - z(k-1)^2)*Q_bar{k-1}(i,j))/2;
             D_temp = D_temp + ((z(k)^3 - z(k-1)^3)*Q_bar{k-1}(i,j))/3;
         end
         A(i,j) = A_temp; %extensional stiffness matrix
         B(i,j) = B_temp; %bending-extension coupling matrix
         D(i,j) = D_temp; %bending stiffness matrix
    end
end
end