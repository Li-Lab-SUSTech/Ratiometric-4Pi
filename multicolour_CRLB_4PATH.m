% Copyright (c) 2021 Li Lab, Southern University of Science and Technology, Shenzhen
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

ratio=1;
NP=1000;
Nphotons = [NP NP NP*ratio NP*ratio];


n = 1;
for z = zrange
    phi0 = [0, pi, pi/2, 1.5 * pi]; 
    
    ground_truth.x(:,1) = Npixels/2 -1 +2*rand([Nfits 1]);
    ground_truth.y(:,1) = Npixels/2 -1 +2*rand([Nfits 1]);
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
   
%     imstack = simSplinePSF(Npixels, PSF, ground_truth.N, ground_truth.bg, coordinates, phi0);%simulate images
%   
%     zstart = [z0-300/dz  z0+300/dz];
%   
%     initZA=repmat(zstart(1),[1 Nfits]);
%     initPhase = [ pi/2 pi/2*3];
%     initPhaseA=repmat(initPhase(1),[1 Nfits]);
%     
%     IAB(:,:,:,:,3)=PSF.Ispline;
%     IAB(:,:,:,:,1)=PSF.Aspline;
%     IAB(:,:,:,:,2)=PSF.Bspline;
    
%      [P, update, error, model] =  kernel_MLEfit_Spline_LM_newFit(imstack(:, :, :, :), PSF, Npixels, 50, phi0, zstart(1));
%     P(:, 6) = wrapTo2Pi(P(:, 6));
%      if length(zstart)>1
%         for i=2:length(zstart) 
%              initZA=repmat(zstart(i),[1 Nfits]);
%              initPhaseA=repmat(initPhase(i),[1 Nfits]);
%               [Ph, updateh, errorh, modelh] =  kernel_MLEfit_Spline_LM_newFit(imstack(:, :, :, :), PSF, Npixels, 50, phi0, zstart(i));
%               Ph(:, 6) = wrapTo2Pi(Ph(:, 6));
% %             [Ph,CRLBh, LogLh] = GPUmleFit_LM_4Pi(single(imstack(:, :, :, :)),uint32(sharedA),100,single(IAB),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
% %             %         indbettero=LogLh<LogL;
%            indbetter=error-errorh>1e-4; %copy only everything if LogLh increases by more than rounding error.
%              P(indbetter,:)=Ph(indbetter,:);
% %              indbetter=LogLh-LogL>1e-4; %copy only everything if LogLh increases by more than rounding error.
% %              P(indbetter,:)=Ph(indbetter,:);
% %              CRLB(indbetter,:)=CRLBh(indbetter,:);
% %              LogL(indbetter)=LogLh(indbetter);
%          end
%      end 
%     X=P(:,1);
%     Y=P(:,2);
%     Z=P(:,5);
%     Znm = (Z-z0)*dz;
%     N1(:,n)=P(:,3);
%     N2(:,n)=P(:,7);
%     BG=P(:,4);
%     phi =P(:,6);
%     phi = wrapTo2Pi(phi);
%     z_phi = z_from_phi_YL(Z, phi, k, z0, dz);
%     
%     fitResults(:,:,n) =[X Y Znm z_phi]; 
%     
%     indgood = (z_phi-ground_truth.znm)<100;
%     
%     s_x_found(n,1)= std(X(indgood)-coordinates(indgood,1,1));
%     s_y_found(n,1)= std(Y(indgood)-coordinates(indgood,2,1));
%     s_z_found(n,1)=std(Znm(indgood)-ground_truth.znm(indgood));
%     s_phiZ_found(n,1) = std(z_phi(indgood)-ground_truth.znm(indgood))
%     
%     RMSEX(n,1) = sqrt(mean((X(indgood)-coordinates(indgood,1,1)).^2));
%     RMSEY(n,1) = sqrt(mean((Y(indgood)-coordinates(indgood,2,1)).^2));
%     RMSEZ(n,1) = sqrt(mean((Znm(indgood)-ground_truth.znm(indgood)).^2));
%     RMSEphiZ(n,1) = sqrt(mean((z_phi(indgood)-ground_truth.znm(indgood)).^2));
   
    shared=[1,1,1,1,1,1];
%     sharedA = repmat(shared',[1 Nfits]);
%     phi0A=repmat(phi0',[1 Nfits]);
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
   
    shared1=[1,1,0,0,1,1];
     CRLB_YL_shared0 = calculate_CRLB_YL_shared(1, PSF, Npixels, phi0, true_theta,shared1);
    offset1 = 6;
    CRLBx0=CRLB_YL_shared0(:,1);
    CRLBy0=CRLB_YL_shared0(:,2);
    CRLBz0=CRLB_YL_shared0(:,5+offset1);
    CRLBphi0=CRLB_YL_shared0(:,6+offset1);
    
    meanCRLBx0(n,1) = mean(sqrt(CRLBx0));
    meanCRLBy0(n,1) = mean(sqrt(CRLBy0));
    meanCRLBz0(n,1) = mean(sqrt(CRLBz0))*dz;
    meanCRLBphiz0(n,1) = mean(sqrt(CRLBphi0))/(2*k);
    
%     CRLB_multicolor = calculate_CRLB(2, PSF, Npixels, phi0, true_theta1);
%     CRLBx=CRLB_multicolor(:,1);
%     CRLBy=CRLB_multicolor(:,2);
%     CRLBz=CRLB_multicolor(:,5);
%     CRLBphi=CRLB_multicolor(:,6);
%     
%     meanCRLBx(n,1) = mean(sqrt(CRLBx));
%     meanCRLBy(n,1) = mean(sqrt(CRLBy));
%     meanCRLBz(n,1) = mean(sqrt(CRLBz))*dz;
%     meanCRLBphiz(n,1) = mean(sqrt(CRLBphi))/(2*k); 
%     
%      phi01 = [0, pi/2, pi, 1.5 * pi];
%     CRLB_multicolor1 = calculate_CRLB(2, PSF, Npixels, phi01, true_theta1);
%     CRLBx1=CRLB_multicolor1(:,1);
%     CRLBy1=CRLB_multicolor1(:,2);
%     CRLBz1=CRLB_multicolor1(:,5);
%     CRLBphi1=CRLB_multicolor1(:,6);
%     
%     meanCRLBx2(n,1) = mean(sqrt(CRLBx1));
%     meanCRLBy2(n,1) = mean(sqrt(CRLBy1));
%     meanCRLBz2(n,1) = mean(sqrt(CRLBz1))*dz;
%     meanCRLBphiz2(n,1) = mean(sqrt(CRLBphi1))/(2*k); 
    
    n = n+1;
end

ratio=1;
NP=1000;
Nphotons = [NP NP NP*ratio NP*ratio];


n = 1;
for z = zrange
    phi0 = [0, pi, pi/2, 1.5 * pi];
    
    ground_truth.x(:,1) = Npixels/2 -1 +2*rand([Nfits 1]);
    ground_truth.y(:,1) = Npixels/2 -1 +2*rand([Nfits 1]);
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
   
    
    CRLB_multicolor = calculate_CRLB(2, PSF, Npixels, phi0, true_theta1);
    CRLBx=CRLB_multicolor(:,1);
    CRLBy=CRLB_multicolor(:,2);
    CRLBz=CRLB_multicolor(:,5);
    CRLBphi=CRLB_multicolor(:,6);
    
    meanCRLBx(n,1) = mean(sqrt(CRLBx));
    meanCRLBy(n,1) = mean(sqrt(CRLBy));
    meanCRLBz(n,1) = mean(sqrt(CRLBz))*dz;
    meanCRLBphiz(n,1) = mean(sqrt(CRLBphi))/(2*k); 
    
     phi01 = [0, pi/2, pi, 1.5 * pi];
    CRLB_multicolor1 = calculate_CRLB(2, PSF, Npixels, phi01, true_theta1);
    CRLBx1=CRLB_multicolor1(:,1);
    CRLBy1=CRLB_multicolor1(:,2);
    CRLBz1=CRLB_multicolor1(:,5);
    CRLBphi1=CRLB_multicolor1(:,6);
    
    meanCRLBx2(n,1) = mean(sqrt(CRLBx1));
    meanCRLBy2(n,1) = mean(sqrt(CRLBy1));
    meanCRLBz2(n,1) = mean(sqrt(CRLBz1))*dz;
    meanCRLBphiz2(n,1) = mean(sqrt(CRLBphi1))/(2*k); 
    
    n = n+1;
end


values1 = spcrv([[zrange(1) zrange zrange(end)];[meanCRLBphiz(1) meanCRLBphiz' meanCRLBphiz(end)]],10);
values2 = spcrv([[zrange(1) zrange zrange(end)];[meanCRLBphiz2(1) meanCRLBphiz2' meanCRLBphiz2(end)]],10);
values3 = spcrv([[zrange(1) zrange zrange(end)];[meanCRLBphiz0(1) meanCRLBphiz0' meanCRLBphiz0(end)]],10);
values4 = spcrv([[zrange(1) zrange zrange(end)];[meanCRLBphiz1(1) meanCRLBphiz1' meanCRLBphiz1(end)]],10);

figure,plot(values1(1,1:2:end),values1(2,1:2:end),'s','Color',[0.5 0.5 0.5]);
hold on
plot(values2(1,:),values2(2,:),'-','Color',[0.5 0.5 0.5]);
plot(values3(1,:),values3(2,:),'-','Color',[0.667 0.667 0.1]);
plot(values4(1,:),values4(2,:),'r-','LineWidth',1);
x_label = xlabel('z position (nm)');
y_label = ylabel('CRLB (nm)');
lgd = legend('CRLBphiZ-0+\pi','CRLBphiZ-0+\pi/2','CRLBphiZ-noshared','CRLBphiZ-allshared','Location','north');
set(lgd,'FontName','time','FontSize',15,'LineWidth',1,'FontWeight','bold');   
set(x_label,'FontName','time','FontSize',20,'LineWidth',3,'FontWeight','bold');
set(y_label,'FontName','time','FontSize',20,'LineWidth',3,'FontWeight','bold');
set(gca,'FontName','time','FontSize',15,'FontWeight','bold');
axis([-600 600 1.6 5.2])
box off

CRLBphiZ_4channel{1} = values1;
CRLBphiZ_4channel{2} = values2;
CRLBphiZ_4channel{3} = values3;
CRLBphiZ_4channel{4} = values4;

% save('CRLBphiZ_4channel.mat','CRLBphiZ_4channel')

values1 = spcrv([[zrange(1) zrange zrange(end)];[meanCRLBx(1) meanCRLBx' meanCRLBx(end)]],10);
values2 = spcrv([[zrange(1) zrange zrange(end)];[meanCRLBx0(1) meanCRLBx0' meanCRLBx0(end)]],10);
values3 = spcrv([[zrange(1) zrange zrange(end)];[meanCRLBx1(1) meanCRLBx1' meanCRLBx1(end)]],10);
values4 = spcrv([[zrange(1) zrange zrange(end)];[meanCRLBx2(1) meanCRLBx2' meanCRLBx2(end)]],10);


values11 = spcrv([[zrange(1) zrange zrange(end)];[meanCRLBy(1) meanCRLBy' meanCRLBy(end)]],10);
values22 = spcrv([[zrange(1) zrange zrange(end)];[meanCRLBy0(1) meanCRLBy0' meanCRLBy0(end)]],10);
values33 = spcrv([[zrange(1) zrange zrange(end)];[meanCRLBy1(1) meanCRLBy1' meanCRLBy1(end)]],10);
values44 = spcrv([[zrange(1) zrange zrange(end)];[meanCRLBy2(1) meanCRLBy2' meanCRLBy2(end)]],10);



CRLBphiX_4channel{1} = values1;
CRLBphiX_4channel{2} = values2;
CRLBphiX_4channel{3} = values3;
CRLBphiX_4channel{4} = values4;

CRLBphiY_4channel{1} = values11;
CRLBphiY_4channel{2} = values22;
CRLBphiY_4channel{3} = values33;
CRLBphiY_4channel{4} = values44;

% save('CRLBphiX_4channel.mat','CRLBphiX_4channel')
% save('CRLBphiY_4channel.mat','CRLBphiY_4channel')

figure,plot(values1(1,:),values1(2,:)*120,'r-');
hold on
plot(values2(1,1:5:end),values2(2,1:5:end)*120,'ro');
plot(values3(1,2:5:end),values3(2,2:5:end)*120,'k.');
plot(values4(1,4:5:end),values4(2,4:5:end)*120,'rs');
plot(values11(1,:),values11(2,:)*120,'-','Color',[0.5 0.5 0.5]);
plot(values22(1,1:5:end),values22(2,1:5:end)*120,'o','Color',[0.667 0.667 0.1]);
plot(values33(1,2:5:end),values33(2,2:5:end)*120,'k.');
plot(values44(1,4:5:end),values44(2,4:5:end)*120,'s','Color',[0.5 0.5 0.5]);
x_label = xlabel('z position (nm)');
y_label = ylabel('CRLB (nm)');
lgd = legend('CRLBX-0+\pi','CRLBX-0+\pi/2','CRLBX-noshared','CRLBX-allshared','CRLBY-0+\pi','CRLBY-0+\pi/2','CRLBY-noshared','CRLBY-allshared','Location','north');
set(lgd,'FontName','time','FontSize',15,'LineWidth',1,'FontWeight','bold');   
set(x_label,'FontName','time','FontSize',20,'LineWidth',3,'FontWeight','bold');
set(y_label,'FontName','time','FontSize',20,'LineWidth',3,'FontWeight','bold');
set(gca,'FontName','time','FontSize',15,'FontWeight','bold');
axis([-600 600 2.5 8])
box off
