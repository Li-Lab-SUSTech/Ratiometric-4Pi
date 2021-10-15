% Copyright (c) 2021 Li Lab, Southern University of Science and Technology, Shenzhen & Ries Lab, European Molecular Biology Laboratory, Heidelberg.
% author: Yiming Li & Jianwei Chen 
% email: liym2019@sustech.edu.cn & 12149038@mail.sustech.edu.cn
% date: 2021.10.07
% Tested with CUDA 11.3 (Express installation) and Matlab 2019a
clear
close all
clc
%% hyper parameters for PSF model used for fit
paraSim.NA = 1.35;                                                % numerical aperture of obj             
paraSim.refmed = 1.450;                                            % refractive index of sample medium
paraSim.refcov = 1.518;                                           % refractive index of converslip
paraSim.refimm = 1.40;                                           % refractive index of immersion oil
paraSim.lambda = 668;                                             % wavelength of emission
paraSim.objStage0_upper = -0;                                        % nm, initial objStage0 position,relative to focus at coverslip
paraSim.objStage0_lower = -0;                                        % nm, initial objStage0 position,relative to focus at coverslip
paraSim.zemit0_upper = -1*paraSim.refmed/paraSim.refimm*(paraSim.objStage0_upper);  % reference emitter z position, nm, distance of molecule to coverslip
paraSim.zemit0_lower = -1*paraSim.refmed/paraSim.refimm*(paraSim.objStage0_lower);  % reference emitter z position, nm, distance of molecule to coverslip


paraSim. pixelSizeX = 120;                                        % nm, pixel size of the image
paraSim. pixelSizeY = 120;                                        % nm, pixel size of the image
paraSim.Npupil = 64;                                             % sampling at the pupil plane

paraSim.aberrations(:,:,1) = [2,-2,0.0; 2,2,-0.06; 3,-1,0.0; 3,1,0.0; 4,0,0; 3,-3,0.0; 3,3,0.0; 4,-2,0.0; 4,2,0.00; 5,-1,0.0; 5,1,0.0; 6,0,0.0; 4,-4,0.0; 4,4,0.0;  5,-3,0.0; 5,3,0.0;  6,-2,0.0; 6,2,0.0; 7,1,0.0; 7,-1,0.00; 8,0,0.0];
paraSim.aberrations(:,:,2) = [2,-2,0.0; 2,2,0.06; 3,-1,0.0; 3,1,0.0; 4,0,0.0; 3,-3,0.0; 3,3,0.0; 4,-2,0.0; 4,2,0.00; 5,-1,0.0; 5,1,0.0; 6,0,0.0; 4,-4,0.0; 4,4,0.0;  5,-3,0.0; 5,3,0.0;  6,-2,0.0; 6,2,0.0; 7,1,0.0; 7,-1,0.00; 8,0,0.0];

paraSim.aberrations(:,3,:) =  paraSim.aberrations(:,3,:)*paraSim.lambda;

paraSim.offset = [0 0];
paraSim.phaseshift = [0 ,pi/2, pi, 3*pi/2];

%% parameters for molecules for simulation
Nmol = 101;
Npixels = 31;
% Nphotons = 5000 +10000*rand(1,Nmol);
% bg = 10 + 10*rand(1,Nmol);
paraSim.Nmol = Nmol;
paraSim.sizeX = Npixels;
paraSim.sizeY = Npixels;

paraSim.xemit = (-200+400*rand(1,Nmol))*0;                             %nm
paraSim.yemit = (-200+400*rand(1,Nmol))*0;                             %nm
paraSim.zemit = linspace(-1000,1000,Nmol)*1;                                      %nm
paraSim.objStage = linspace(-1000,1000,Nmol)*0;                                  %nm

[PSFs,PSFsUpper,PSFsLower,WaberrationUpper, WaberrationLower] = vectorPSF_4Pi(paraSim);


% data = PSFs;
% for i = 1:Nmol
%     data(:,:,i) = Nphotons(i)*PSFs(:,:,i)+bg(i);
% end
% data = poissrnd(data,size(data));

%% generate IAB model
% PSF = struct('ipalm_im', nan, 'I', nan, 'A', nan, 'B', nan, 'Ispline', nan, ...
%     'Aspline', nan, 'Bspline', nan, 'k', nan, 'zstep', nan, 'zcand', nan, 'imsz', nan);
% 
% ipalm_im = allPSFs * 4 / sum(reshape(allPSFs(:, :, round(length(zcand) / 2), :), 1, [])); %normalize PSF
% imageslicer(ipalm_im);
% PSF.ipalm_im=ipalm_im;


%phaseshift, lambdanm and zcand (z range in nm) are 
%saved in the file with the zstack
ipalm_im  = PSFs;

phaseshift = paraSim.phaseshift;
% k = 2 * pi*paraSim.refimm / paraSim.lambda; %lambda in nm
k = 2 * pi / paraSim.lambda; %lambda in nm
zcand = paraSim.zemit;% if move obj stage use paraSim.objStage
zstep = zcand(2) - zcand(1);

imsz = paraSim.sizeX;

%PSF.ipalm_im = PSF.ipalm_im * 4 / sum(reshape(PSF.ipalm_im(:, :, round(length(PSF.zcand) / 2), :), 1, [])); %normalize PSF
%imageslicer(PSF.ipalm_im);

%make I, A, B and their splines
I = squeeze((ipalm_im(:, :, :, 1) + ipalm_im(:, :, :, 3)) / 2);
% I2 = squeeze((ipalm_im(:, :, :, 2) + ipalm_im(:, :, :, 4)) / 2);
% imageslicer(PSF.I)
% imageslicer(I2)
% imageslicer(I2 - PSF.I)

kz2 = 2 * k * zcand';
kz2 = permute(repmat(kz2, 1, imsz, imsz), [2, 3, 1]);

F_phi1 = squeeze(ipalm_im(:, :, :, 1)) - I;
F_phi2 = squeeze(ipalm_im(:, :, :, 2)) - I;
phi1 = phaseshift(1);
phi2 = phaseshift(2);

A = (F_phi1 .* sin(kz2 + phi2) - F_phi2 .* sin(kz2 + phi1)) / sin(phi2 - phi1);
B = (-F_phi1 .* cos(kz2 + phi2) + F_phi2 .* cos(kz2 + phi1)) / sin(phi2 - phi1);
% imageslicer(PSF.A);
% imageslicer(PSF.B);

%check: restore PSF in the 4th quadrant
phi4 = phaseshift(4);
check_PSF = I + A .* cos(kz2 + phi4) + B .* sin(kz2 + phi4);
% imageslicer(squeeze(PSF.ipalm_im(:, :, :, 4)));
% imageslicer(check_PSF);
% dipshow([check_PSF PSFs(:,:,:,4) (check_PSF-PSFs(:,:,:,4))]);
Ispline = Spline3D_interp(I);
Aspline = Spline3D_interp(A);
Bspline = Spline3D_interp(B);

%% test fitter
lambdanm = paraSim.lambda;
dz = zcand(2) - zcand(1);
z0 = round(length(zcand) / 2); 
% k = 2 * pi * paraSim.refimm / lambdanm; %lambda in nm
NV = 6;
PSF.Ispline = Ispline;
PSF.Aspline = Aspline;
PSF.Bspline = Bspline;
Nfits = 2;
Npixels = 13;
bg = [20 20 20 20];
zrange = -600:20:600;

ratio=0.4377;
NP=1000;
Nphotons = [NP*ratio NP*ratio NP*ratio NP*ratio];


n = 1;
for z = zrange
    phi0 = [0, pi/2, pi, 1.5 * pi];
    
    ground_truth.x(:,1) = Npixels/2 -1 +0*rand([Nfits 1]);
    ground_truth.y(:,1) = Npixels/2 -1 +0*rand([Nfits 1]);
    ground_truth.N = repmat(Nphotons, [Nfits 1]);
    ground_truth.bg =repmat(bg, [Nfits 1]);
    ground_truth.znm = z*ones(Nfits,1);
    ground_truth.zspline = ground_truth.znm / dz + z0;
    ground_truth.phi =  wrapTo2Pi(2 * k * ground_truth.znm);
    
    scale =0;
    ground_truth.x(:,2)=ground_truth.x(:,1)+(rand(Nfits, 1))*scale;
    ground_truth.y(:,2)=ground_truth.y(:,1)+(rand(Nfits, 1))*scale;
    
    ground_truth.x(:,3)=ground_truth.x(:,1)+(rand(Nfits, 1))*scale;
    ground_truth.y(:,3)=ground_truth.y(:,1)+(rand(Nfits, 1))*scale;
    
    ground_truth.x(:,4)=ground_truth.x(:,1)+(rand(Nfits, 1))*scale;
    ground_truth.y(:,4)=ground_truth.y(:,1)+(rand(Nfits, 1))*scale;
    
    % coordinates for simulation
    coordinates = zeros(Nfits,4,length(phi0));
    for kk=1:1:length(phi0)
        coordinates(:,:,kk) = [ground_truth.x(:,kk) ground_truth.y(:,kk) ground_truth.zspline ground_truth.phi];
    end
    
    true_theta = zeros(Nfits,length(phi0),NV);
    true_theta(:,:,1) = ground_truth.x;
    true_theta(:,:,2) = ground_truth.y;
    true_theta(:,:,3) = ground_truth.N;
    true_theta(:,:,4) = ground_truth.bg;
    true_theta(:,:,5) = repmat(ground_truth.zspline,1,4);
    true_theta(:,:,6) = repmat(ground_truth.phi,1,4);
    true_theta1=squeeze(true_theta(:,1,:));
    true_theta1(:,7)=squeeze(true_theta(:,3,3));
    true_theta =permute(true_theta,[3 2 1]);
 
    shared=[1,1,1,1,1,1];

    CRLB_YL_shared = calculate_CRLB_YL_shared(1, PSF, Npixels, phi0, true_theta,shared);
    offset = 0;
    CRLBx1=CRLB_YL_shared(:,1);
    CRLBy1=CRLB_YL_shared(:,2);
    CRLBz1=CRLB_YL_shared(:,5+offset);
    CRLBphi1=CRLB_YL_shared(:,6+offset);
    
    meanCRLBx1(n,1) = mean(sqrt(CRLBx1));
    meanCRLBy1(n,1) = mean(sqrt(CRLBy1));
    meanCRLBz1(n,1) = mean(sqrt(CRLBz1))*dz;
    meanCRLBphiz1(n,1) = mean(sqrt(CRLBphi1))/(2*k);

    n = n+1;
end

ratio1=0.4264;
ratio2=0.7448;  
NP=1000;
Nphotons = [NP*ratio1 NP*ratio1 NP*ratio2 NP*ratio2];


n = 1;
for z = zrange
    phi0 = [0, pi, pi/2, 1.5 * pi];
    
    ground_truth.x(:,1) = Npixels/2 -1 +0*rand([Nfits 1]);
    ground_truth.y(:,1) = Npixels/2 -1 +0*rand([Nfits 1]);
    ground_truth.N = repmat(Nphotons, [Nfits 1]);
    ground_truth.bg =repmat(bg, [Nfits 1]);
    ground_truth.znm = z*ones(Nfits,1);
    ground_truth.zspline = ground_truth.znm / dz + z0;
    ground_truth.phi =  wrapTo2Pi(2 * k * ground_truth.znm);
    
    scale =0;
    ground_truth.x(:,2)=ground_truth.x(:,1)+(rand(Nfits, 1))*scale;
    ground_truth.y(:,2)=ground_truth.y(:,1)+(rand(Nfits, 1))*scale;
    
    ground_truth.x(:,3)=ground_truth.x(:,1)+(rand(Nfits, 1))*scale;
    ground_truth.y(:,3)=ground_truth.y(:,1)+(rand(Nfits, 1))*scale;
    
    ground_truth.x(:,4)=ground_truth.x(:,1)+(rand(Nfits, 1))*scale;
    ground_truth.y(:,4)=ground_truth.y(:,1)+(rand(Nfits, 1))*scale;
    
    % coordinates for simulation
    coordinates = zeros(Nfits,4,length(phi0));
    for kk=1:1:length(phi0)
        coordinates(:,:,kk) = [ground_truth.x(:,kk) ground_truth.y(:,kk) ground_truth.zspline ground_truth.phi];
    end
    
    true_theta = zeros(Nfits,length(phi0),NV);
    true_theta(:,:,1) = ground_truth.x;
    true_theta(:,:,2) = ground_truth.y;
    true_theta(:,:,3) = ground_truth.N;
    true_theta(:,:,4) = ground_truth.bg;
    true_theta(:,:,5) = repmat(ground_truth.zspline,1,4);
    true_theta(:,:,6) = repmat(ground_truth.phi,1,4);
    true_theta1=squeeze(true_theta(:,1,:));
    true_theta1(:,7)=squeeze(true_theta(:,3,3));
    true_theta =permute(true_theta,[3 2 1]);
   
%    
    CRLB_multicolor = calculate_CRLB(2, PSF, Npixels, phi0, true_theta1);
    CRLBx=CRLB_multicolor(:,1);
    CRLBy=CRLB_multicolor(:,2);
    CRLBz=CRLB_multicolor(:,5);
    CRLBphi=CRLB_multicolor(:,6);
    
    meanCRLBx(n,1) = mean(sqrt(CRLBx));
    meanCRLBy(n,1) = mean(sqrt(CRLBy));
    meanCRLBz(n,1) = mean(sqrt(CRLBz))*dz;
    meanCRLBphiz(n,1) = mean(sqrt(CRLBphi))/(2*k); 

    n = n+1;
end


figure,plot(zrange,meanCRLBphiz,'-','Color',[0.667 0.667 0.1]);
hold on
plot(zrange,meanCRLBphiz1,'s','Color',[0.667 0.667 0.1]);
% grid on

plot(zrange,meanCRLBx*120,'r-');
plot(zrange,meanCRLBx1*120,'rs');

plot(zrange,meanCRLBy*120,'-','Color',[0.5 0.5 0.5]);
plot(zrange,meanCRLBy1*120,'s','Color',[0.5 0.5 0.5]);

x_label = xlabel('z position (nm)');
y_label = ylabel('CRLB^{1/2} (nm)');
lgd = legend('{\itz}_\phi','{\itz}_\phi-salvaged',' {\itx}',' {\itx}-salvaged',' {\ity}',' {\ity}-salvaged','Location','north');
title('DY633','FontSize',22,'FontWeight','bold'); 
set(lgd,'FontName','time','FontSize',15,'LineWidth',1,'FontWeight','bold');   
set(x_label,'FontName','time','FontSize',22,'LineWidth',3,'FontWeight','bold');
set(y_label,'FontName','time','FontSize',22,'LineWidth',3,'FontWeight','bold');
set(gca,'FontName','time','FontSize',17,'FontWeight','bold');
axis([-600 600 0 14])
box off
legend('boxoff')
