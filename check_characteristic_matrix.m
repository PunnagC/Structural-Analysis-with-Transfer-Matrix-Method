function modified_char_matrix = check_characteristic_matrix(Char_Matrix_w,w)
%% This function checks if the characteristic matrix has nonzero elements on the top most rows
% If not, it reduces the 3x3 matrix to 2x2 matrix with non-zero elements
% such that successful evaluation of determinant can take place.
% this is useful for bending TMM as dending TMM can have 3x3 and 2x2
% characteristic matrices based on the type of FTM size (6x6 OR 4x4)
%----------------------------------INPUT-----------------------------------
% Char_Matrix_w - characteristic matrix, used for evaluating characteristic
%                 equation (6x6 OR 4x4)
% w             - symbolic variable, circular frequency (rad/s)
%---------------------------------OUTPUT-----------------------------------
% modified_char_matrix - modified characteristic matrix
%--------------------------------------------------------------------------
%% Actual code below
assume(w > 0)
[nrow,~] = size(Char_Matrix_w);
if nrow > 2
    FTM_size = 2;
else
    FTM_size = 1;
end

if FTM_size == 2
    Char_Matrix_w = vpa(Char_Matrix_w,5);
    Det_F         = det(Char_Matrix_w);
    if Det_F == 0
        modified_char_matrix = Char_Matrix_w([2,3],[2,3]);
        disp('Modified characteristic matrix evaluated')
    else
        modified_char_matrix = Char_Matrix_w;
    end
else
    modified_char_matrix = Char_Matrix_w;
    
end

end