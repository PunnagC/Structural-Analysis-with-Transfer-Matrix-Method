function PHI_shape = get_sym_bending_modeshape(w_eval,w,FTM_x_w,FTM_x_w_B,x,l_s,BC,FTM_size)
%% EVALUATES THE SYMBOLIC MODE SHAPE
%This function plots the mode shape corresponding to a natural frequency
%----------------------------------INPUT-----------------------------------
% w_eval   - scalar, current natural freq. (rad/s) for which mode shape is reqd.
% w        - scalar, symbolic natural freq.
% FTM_x_w  - Field transfer matrix (symbolic in 'x' and 'w')
% x        - symbolic beam local coordinate along length
% l_s      - [1xns] length of each segment
% BC       - string, geometric BC on either side ('clamped','pinned','free')
%            example: 'clamped-clamped' OR 'clamped-free'
% FTM_size - whether to perform TMM in 4x4 or more general 6x6 format
%---------------------------------OUTPUT-----------------------------------
% PHI_shape - [1xns] segmental mode shape corresponding to 'w_eval' natural freq.
%--------------------------------------------------------------------------
% keyboard
[nrow, ncol, ns] = size(FTM_x_w);
FTM_w            = sym(NaN(nrow,ncol,ns));
% [N,~]            = size(FTM_x_w);
mul              = 1;
for j = 1:ns
    FTM_w(:,:,j) = subs(FTM_x_w(:,:,j),x,l_s(j));
    termj = subs(FTM_w(:,:,j),w,w_eval);%( FTM_w,freq_w,w,ns,BC)
    mul   = termj*mul;
end
FTM = subs(FTM_w,w,w_eval); %3D matrix 'FTM' is just numbers
GTM = double(mul); 
% GTM = vpa(GTM,5);
   
[U, ~,root_cols] = TMM_characteristic_equation(BC,GTM,'bend',FTM_size);
PHI_shape        = sym(NaN(1,ns));

if nrow == 4 %for 4x4 matrix
    if U(1,2) ~= 0
        k = - U(1,1)/U(1,2);
    elseif U(2,2) ~= 0
        k = - U(2,1)/U(2,2);
    else
        k = 0;
    end
    k = double(k);
elseif nrow == 6 %for 6x6 matrix
    if U(2,1) ~= 0
        k1    = U(2,2)/U(2,1);
        k2    = U(2,3)/U(2,1);
        sigma = (U(3,1)*U(2,2) - U(3,2)*U(2,1))/(U(3,3)*U(2,1)  - U(3,1)*U(2,3));
    elseif U(3,1) ~= 0
        k1    = U(3,2)/U(3,1);
        k2    = U(3,3)/U(3,1);
        sigma = (U(3,1)*U(2,2) - U(3,2)*U(2,1))/(U(3,3)*U(2,1)  - U(3,1)*U(2,3));
    else
        k1    = 0;
        k2    = 0;
        sigma = - U(3,2)/U(3,3);
    end
    %       if U(3,3) ~= 0
    %           k3 = -(U(3,1) + U(3,2))/U(3,3);
else
    error('Please change/modify the code')
end
% end
    
%% Current segment is just the beam and all previous segments are beam + point mass FTM's
% keyboard
for seg = 1:ns
    GTM_seg_x = eye(nrow); %4x4 matrix
    if seg > 1
        for i = 1: seg - 1 %segmental sweep before current segment
            GTM_seg_x = FTM(:,:,i)*GTM_seg_x; %at this point this matrix 'GTM_seg_x' is just numbers
        end
        GTM_seg_x = subs(FTM_x_w_B(:,:,seg),w,w_eval)*GTM_seg_x; %at this point this matrix 'GTM_seg_x' is just a function of 'x'
    else
        GTM_seg_x = subs(FTM_x_w_B(:,:,seg),w,w_eval); %for a structure with just 1 segment
        GTM_seg_x = vpa(GTM_seg_x); %4x4 matrix OR 6x6 matrix 
    end
    if nrow == 4 %for 4x4 matrix
        PHI_shape(seg) = GTM_seg_x(1,root_cols)*[1; k];
    elseif nrow == 6 %for 6x6 matrix
        PHI_shape(seg) = -GTM_seg_x(3,root_cols(1))*(k1 + sigma*k2) + GTM_seg_x(3,root_cols(2)) + GTM_seg_x(3,root_cols(3))*sigma;
    end
end
 
end


