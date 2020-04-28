
function [matl] = read_material_data(mat_pzt)
fid   = fopen(mat_pzt,'r');
fgetl(fid);
mat_type = fgetl(fid);
delimiter  = {'\t','='};
% formatSpec = '%s%f%s%*s%[^\n\r]';
if strcmp(mat_type,'Orthotropic') == 1
    M = NaN(1,9);
    for i = 1:9
        stri         = fgetl(fid);
        split_string = strsplit(stri, delimiter);
        M(i)         = str2double(split_string{2});
    end
    
    
elseif strcmp(mat_type,'Isotropic') == 1
    M = NaN(1,2);
    for i = 1:2
        stri         = fgetl(fid);
        split_string = strsplit(stri, delimiter);
        M(i)         = str2double(split_string{2});
    end
    
end

rho_temp = strsplit(fgetl(fid), delimiter);
rho      = str2double(rho_temp{2});

fgetl(fid);

cE_temp  = strsplit(fgetl(fid), delimiter);
cE       = str2double(cE_temp{2});

d31_temp = strsplit(fgetl(fid), delimiter);
d31      = str2double(d31_temp{2});

eS_temp  = strsplit(fgetl(fid), delimiter);
er       = str2double(eS_temp{2});



matl.M    = M;
matl.rho  = rho;
matl.type = mat_type;

matl.cE   = cE;
matl.d31  = d31;
matl.er   = er;


fclose('all');

end
