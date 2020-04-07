function [PHI_shape_x] =  get_sym_torsional_modeshape(w_eval,w,FTM_x_w,FTM_x_w_B,x,l_s)
%Mode_shape_calc_heave_V1(plot_freq,f,FTM_x_f,x,FTM_f,Li,divisions)
%This function finds the mode shape corresponding to a natural frequency
[nrow, ncol, ns] = size(FTM_x_w);
% ns = length(l_seg);     %number of segments in the structure
PHI_shape_x(1:ns) = sym(NaN);
[N,~] = size(FTM_x_w);

FTM_w            = sym(NaN(nrow,ncol,ns));

mul              = 1;
for j = 1:ns
    FTM_w(:,:,j) = subs(FTM_x_w(:,:,j),x,l_s(j));
    termj = subs(FTM_w(:,:,j),w,w_eval);%( FTM_w,freq_w,w,ns,BC)
    mul   = termj*mul;
end
FTM = subs(FTM_w,w,w_eval); %3D matrix 'FTM' is just numbers


for seg = 1:ns %segmental sweep
    GTM_seg_x = eye(N);
    if seg > 1
        for i = 1: seg - 1 %segmental sweep before current segment
            GTM_seg_x = FTM(:,:,i)*GTM_seg_x; %at this point this matrix 'GTM_seg_x' is just numbers
        end
        GTM_seg_x = subs(FTM_x_w_B(:,:,seg),w,w_eval)*GTM_seg_x; %at this point this matrix 'GTM_seg_x' is just a function of 'x'
    else
        GTM_seg_x = subs(FTM_x_w_B(:,:,seg),w,w_eval); %for a structure with just 1 segment
        GTM_seg_x = vpa(GTM_seg_x); %4x4 matrix OR 6x6 matrix 
    end
    PHI_shape_x(seg) = GTM_seg_x(1,2);
end
end
    
