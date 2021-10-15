function [dudt,model] =  kernel_DerivativeSpline_v3(xc,yc,zc,xsize,ysize,zsize,delta_f,delta_dxf,delta_dyf,delta_dzf,PSF,theta, NV, phi0)
dudt = zeros(NV, 1);
dudtI = zeros(3, 1); %derivetives of I with respect to x, y and z
dudtA = zeros(3, 1);
dudtB = zeros(3, 1);
model = zeros(1, 3);
xc = max(xc, 0);
xc = min(xc, xsize - 1);

yc = max(yc, 0);
yc = min(yc, ysize - 1);

zc = max(zc, 0);
zc = min(zc, zsize - 1);

temp = PSF.Ispline(xc + 1, yc + 1, zc + 1, :);
dudtI(1) = -1 * sum(delta_dxf .* (temp(:)));
dudtI(2) = -1 * sum(delta_dyf .* (temp(:)));
dudtI(3) = sum(delta_dzf.*(temp(:)));
I = sum(delta_f .* (temp(:)));


temp = PSF.Aspline(xc + 1, yc + 1, zc + 1,:);
dudtA(1) = -1 * sum(delta_dxf .* (temp(:)));
dudtA(2) = -1 * sum(delta_dyf .* (temp(:)));
dudtA(3) = sum(delta_dzf.*(temp(:)));
A = sum(delta_f.*(temp(:)));

temp = PSF.Bspline(xc + 1, yc + 1, zc + 1, :);
dudtB(1) = -1 * sum(delta_dxf .* (temp(:)));
dudtB(2) = -1 * sum(delta_dyf .* (temp(:)));
dudtB(3) = sum(delta_dzf.*(temp(:)));
B = sum(delta_f .* (temp(:)));


for ii = 1:2
    dudt(1, ii) = theta(3) * (dudtI(1) + dudtA(1) * cos(theta(6) + phi0(ii)) + dudtB(1) * sin(theta(6) + phi0(ii)));
    dudt(2, ii) = theta(3) * (dudtI(2) + dudtA(2) * cos(theta(6) + phi0(ii)) + dudtB(2) * sin(theta(6) + phi0(ii)));
    dudt(5, ii)= theta(3) * (dudtI(3) + dudtA(3) * cos(theta(6) + phi0(ii)) + dudtB(3) * sin(theta(6) + phi0(ii)));
    dudt(3, ii) = I +  A * cos(theta(6) + phi0(ii)) + B * sin(theta(6) + phi0(ii));
    dudt(4, ii) = 1;
    dudt(6, ii) = theta(3) * (- A * sin(theta(6) + phi0(ii)) + B * cos(theta(6) + phi0(ii)));
    model (ii)= theta(4) + theta(3) .* dudt(3, ii);
end
for ii = 3:3
    dudt(1, ii) = theta(7) * (dudtI(1) + dudtA(1) * cos(theta(6) + phi0(ii)) + dudtB(1) * sin(theta(6) + phi0(ii)));
    dudt(2, ii) = theta(7) * (dudtI(2) + dudtA(2) * cos(theta(6) + phi0(ii)) + dudtB(2) * sin(theta(6) + phi0(ii)));
    dudt(5, ii)= theta(7) * (dudtI(3) + dudtA(3) * cos(theta(6) + phi0(ii)) + dudtB(3) * sin(theta(6) + phi0(ii)));
    dudt(7, ii) = I +  A * cos(theta(6) + phi0(ii)) + B * sin(theta(6) + phi0(ii));
    dudt(4, ii) = 1;
    dudt(6, ii) = theta(7) * (- A * sin(theta(6) + phi0(ii)) + B * cos(theta(6) + phi0(ii)));
    model (ii)= theta(4) + theta(7) .* dudt(7, ii);
end



