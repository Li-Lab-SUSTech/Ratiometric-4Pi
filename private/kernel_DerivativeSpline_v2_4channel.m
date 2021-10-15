function [dudt,model] =  kernel_DerivativeSpline_v2_4channel(xc,yc,zc,xsize,ysize,zsize,delta_f,delta_dxf,delta_dyf,delta_dzf,PSF,theta, NV, phi0)
dudt = zeros(NV, 1);
dudtI = zeros(3, 1); %derivetives of I with respect to x, y and z
dudtA = zeros(3, 1);
dudtB = zeros(3, 1);

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



    dudt(1) = theta(3) * (dudtI(1) + dudtA(1) * cos(theta(6) + phi0) + dudtB(1) * sin(theta(6) + phi0));
    dudt(2) = theta(3) * (dudtI(2) + dudtA(2) * cos(theta(6) + phi0) + dudtB(2) * sin(theta(6) + phi0));
    dudt(5)= theta(3) * (dudtI(3) + dudtA(3) * cos(theta(6) + phi0) + dudtB(3) * sin(theta(6) + phi0));
    dudt(3) = I +  A * cos(theta(6) + phi0) + B * sin(theta(6) + phi0);
    dudt(4) = 1;
    dudt(6) = theta(3) * (- A * sin(theta(6) + phi0) + B * cos(theta(6) + phi0));

model = theta(4) + theta(3) .* dudt(3);