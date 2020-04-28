function plot_modeS = mode_shape_plot_panel(figno,panel,mode_type,x_vec,BC,freq,T)
%% Error checking at INPUT
argmin = 6;
argmax = 7;
narginchk(argmin,argmax);

%%
if strcmp('Pitch', mode_type) == 1
   [~, n] = size(panel.PSI);
   mode_shape = panel.PSI;
   legend_type = 'na = ';
   mode_symbol = '\psi(x)';
elseif strcmp('Heave', mode_type) == 1
   [~, n] = size(panel.PHI);
   mode_shape = panel.PHI;
   legend_type = 'nh = ';
   mode_symbol = '\phi_w(x)';
   elseif strcmp('Lag', mode_type) == 1
   [~, n] = size(panel.PHI);
   mode_shape = panel.PHI;
   legend_type = 'nlag = ';
   mode_symbol = '\phi_u(x)';
else
    disp('Mode type can only be "Pitch" or "Heave" or "Lag"!');
end
% keyboard
plot_modeS = figure(figno);
set(plot_modeS,'defaulttextinterpreter','latex')
for i = 1:n %number of modes
    plot(x_vec, mode_shape(:,i),'LineWidth',2);
%     legend_stringh(i,:) = [legend_type, num2str(i), ';f=',num2str(freq(i)),'(Hz)'] ;
    legend_stringh{i} = [legend_type, num2str(i), ';f=',num2str(freq(i)),'(Hz)'] ;
    hold on
end
grid on
legend(legend_stringh{:});
if nargin < argmax
    T = 0;
end
title(['[',BC,']',char(mode_type),'\,\,Mode\,\,Shape\,\,Plot',' at T = ', num2str(T),' N'])
xlabel('$Structure\,\,Length(m)$')
ylabel(['$',char(mode_type),'\,\,mode\,\,shape\,\,',mode_symbol,'$']) %'$Normalized\,\,',
set(legend,'Location','Best','LineWidth',1,'FontSize',9)
hold off
    
end