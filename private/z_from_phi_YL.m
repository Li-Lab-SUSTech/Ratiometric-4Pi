function [ z_phi,dphi ] = z_from_phi_YL( z_ast, phi, k, z0, dz )
    lambdanm = 2 * pi / k;
    z_ast_nm = (z_ast - z0) * dz;
    z_ast_phi = z_ast_nm*2*k;
    dphi = -(phi-z_ast_phi)/pi/2;
    period = round(dphi);
    z_phi = period*lambdanm/2 + phi/2/k;
    
    
    
    
%     
%     num_of_periods = fix(z_ast * 2 / lambdanm);
%     num_of_periods(z_ast < 0) = num_of_periods(z_ast < 0) - 1;
%     z_phi =  num_of_periods  * lambdanm / 2 + phi / ( 2 * k);
end

