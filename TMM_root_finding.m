function [w_h, nroots] = TMM_root_finding(w_guess,w,Char_Matrix_w,n_roots,w_inc,method, w_guess111)
assume(w > 0)
% keyboard
fprintf('Root finding method selected = %d \n',method);
% Char_Matrix_w = Char_Matrix_check(Char_Matrix_w,w);
Det_F         = Det_F_calc(Char_Matrix_w,w);
Det_F         = vpa(Det_F,5);

w_h    = [];

% keyboard
w_guess111 = [w_guess111, w_guess111 + 5];
tim       = 0;
long_flag = 0;
while numel(w_h) < n_roots
    tic
    switch method
        case 1 %Root finding within given range
           clear  mag_mF
%             mF         = matlabFunction(vpa(Char_Matrix_w,8));%simplify(expand(

            x_vec      = w_guess(1):0.1:w_guess(2);
            nx         = numel(x_vec);
            mag_mF     = NaN(1,nx);
            root_count = 0;
%             dummy_start_val = 
            for i = 1:nx
                  tempcalc  = subs(Det_F,w,x_vec(i));
                  mag_mF(i) = double(tempcalc);%det(mF(x_vec(i))) 
                  
                  % This part of the code is very useful if just 1 root is
                  % required
                  if i == 1
                      signi(1) = sign(mag_mF(i));
                  elseif i > 1
                      signi(i) = sign(mag_mF(i));
                  end
                  nUniq = numel(unique(signi(1:i)));
                  if nUniq > 1
                      x_vec  = x_vec(1:i);
                      mag_mF = mag_mF(1:i);
                      break;
                  end
                  
            end
%             keyboard
            [~,t_temp] = crossing(mag_mF,x_vec);%, zeros(numel(x_vec),1)
            w_h        = uniquetol(sort([w_h',t_temp]),1e-3);
            w_h        = nonzeros(w_h);
%                         plot(x_vec,mag_mF)
%                                         keyboard
            w_guess(1) = w_guess(2);
            w_guess(2) = w_guess(2) + w_inc;%w_inc
            
        case 2 %Root finding with vpasolve
%             keyboard
            w_h_temp   = vpasolve(Det_F == 0,w,[w_guess(1),w_guess(2)],'random',true);
            w_h_temp   = double(abs(w_h_temp));
            w_h        = uniquetol(sort([w_h,w_h_temp]),1e-3);
            w_h        = nonzeros(w_h);
            w_h        = transpose(nonzeros(w_h));
            w_guess(1) = w_guess(2);
            w_guess(2) = w_guess(2) + w_inc;
            
        case 3 %Root finding with fzero
%             keyboard
%             [w_h_temp,fval] = fzero(@(w0) fzero_roots(Det_F,w,w0),w_guess111(1));
            [w_h_temp,fval] = fzero(@(w0) double(subs(Det_F,w,w0)),w_guess111(1));
%             fval
            
            w_h_temp      = double(abs(w_h_temp));
            w_h           = uniquetol(sort([w_h,w_h_temp]),1e-3);
            w_h           = nonzeros(w_h);
            w_h           = transpose(nonzeros(w_h));
            w_guess111(1) = w_guess111(2);
            w_guess111(2) = w_guess111(1) + w_inc*1.5;
            
        otherwise
            error('Method can be either 1 or 2')
    end
    tstep = toc;
    tim   = tim + tstep;
    if tim >= 60 && long_flag == 0
       cprintf('SystemCommands','TMM_root_finding.m running for >= 60 (secs)\n');
       cprintf('SystemCommands','Try changing the root finding method(options are 1, 2, 3)\n');
       long_flag = long_flag + 1;
    end
       
end

w_h    = (w_h(1:n_roots))';
nroots = numel(w_h);



if isempty(w_h) == 1
    error('No frequencies found within search range !')
end
disp('Finished solving Eigenvalue problem.')
% w_h

end


function Det_F = Det_F_calc(Char_Matrix_w,w)
assume(w > 0)
[nrow,~] = size(Char_Matrix_w);
if nrow > 2
    FTM_size = 2;
elseif nrow == 2
    FTM_size = 1;
else
    FTM_size = 0;
end
% keyboard
Char_Matrix_w = vpa(Char_Matrix_w,5);
if FTM_size == 1
    Det_F = Char_Matrix_w(1,1)*Char_Matrix_w(2,2) - Char_Matrix_w(1,2)*Char_Matrix_w(2,1);
else
    Det_F = det(vpa(Char_Matrix_w,5));
end
Det_F = vpa(Det_F,5);
% Det_F = simplify(expand(Det_F));
disp('Characteristic determinant evaluated')
end



