% Copyright (c) 2021 Li Lab, Southern University of Science and Technology, Shenzhen
% author: Yiming Li & Jianwei Chen 
% email: liym2019@sustech.edu.cn & 12149038@mail.sustech.edu.cn
% date: 2021.10.07
% Tested with CUDA 11.3 (Express installation) and Matlab 2019a
%%
clear
close all
clc
%% hyper parameters for PSF model used for fit
paraSim.NA = 1.35;                                                % numerical aperture of obj             
paraSim.refmed = 1.40;                                            % refractive index of sample medium
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

paraSim.aberrations(:,:,1) = [2,-2,0.0; 2,2,-0.1; 3,-1,0.0; 3,1,0.0; 4,0,0; 3,-3,0.0; 3,3,0.0; 4,-2,0.0; 4,2,0.00; 5,-1,0.0; 5,1,0.0; 6,0,0.0; 4,-4,0.0; 4,4,0.0;  5,-3,0.0; 5,3,0.0;  6,-2,0.0; 6,2,0.0; 7,1,0.0; 7,-1,0.00; 8,0,0.0];
paraSim.aberrations(:,:,2) = [2,-2,0.0; 2,2,0.1; 3,-1,0.0; 3,1,0.0; 4,0,0.0; 3,-3,0.0; 3,3,0.0; 4,-2,0.0; 4,2,0.00; 5,-1,0.0; 5,1,0.0; 6,0,0.0; 4,-4,0.0; 4,4,0.0;  5,-3,0.0; 5,3,0.0;  6,-2,0.0; 6,2,0.0; 7,1,0.0; 7,-1,0.00; 8,0,0.0];

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


%make I, A, B and their splines
I = squeeze((ipalm_im(:, :, :, 1) + ipalm_im(:, :, :, 3)) / 2);


kz2 = 2 * k * zcand';
kz2 = permute(repmat(kz2, 1, imsz, imsz), [2, 3, 1]);

F_phi1 = squeeze(ipalm_im(:, :, :, 1)) - I;
F_phi2 = squeeze(ipalm_im(:, :, :, 2)) - I;
phi1 = phaseshift(1);
phi2 = phaseshift(2);

A = (F_phi1 .* sin(kz2 + phi2) - F_phi2 .* sin(kz2 + phi1)) / sin(phi2 - phi1);
B = (-F_phi1 .* cos(kz2 + phi2) + F_phi2 .* cos(kz2 + phi1)) / sin(phi2 - phi1);


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
Nfits = 1000;
Npixels = 13;
ratio=0.9357;
bg = [20 20 20 20];
phi0 = [0, pi/2, pi, 1.5 * pi];
distribution = 'Photons_AF647 52%, CF660C 75%, DY634 33%, DL650 67%, CF680 79%';
load([distribution '.mat'])
[CF680_no, CF680_xo]=hist(CF680,1000);
distribution1=zeros(1000,2);
distribution1(:,1) = CF680_xo';
distribution1(:,2) = CF680_no';
N = rand_arb_cjw(Nfits, distribution1, 0);

    ground_truth.x(:,1) = Npixels/2 -1 +2*rand([Nfits 1]);
    ground_truth.y(:,1) = Npixels/2 -1 +2*rand([Nfits 1]);
    ground_truth.N  = [N N N*ratio N*ratio];
    ground_truth.bg =repmat(bg, [Nfits 1]);
    ground_truth.znm = -500+1000*rand([Nfits 1]);
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
   
    imstack = simSplinePSF(Npixels, PSF, ground_truth.N, ground_truth.bg, coordinates, phi0);%simulate images
  
    zstart = [z0-300/dz  z0+300/dz];
   numberchannel = 4;
    dTAll=zeros(Nfits,NV,numberchannel);
    for m =1:Nfits
        for i=1:2
            for j=2:numberchannel
                
                dTAll(m,i,j)= coordinates(m,i,j)-coordinates(m,i,1);
            end
        end
    end
    

    
    dTAll=permute(dTAll,[2 3 1]);
    initZA=repmat(zstart(1),[1 Nfits]);
    initPhase = [ pi/2 pi/2*3];
    initPhaseA=repmat(initPhase(1),[1 Nfits]);
    
    IAB(:,:,:,:,3)=PSF.Ispline;
    IAB(:,:,:,:,1)=PSF.Aspline;
    IAB(:,:,:,:,2)=PSF.Bspline;
    shared=[1 1 2 0 1 1]; % shared = [x, y, NP, bg, z, phi]. when shared(i)=1, link all channels. when shared(i)=0, individual fit. when shared(i)=2, link channel 1 and 3.
    if shared(3)==0
        offset= 4*(2-(shared(3))-shared(4));
    else
    offset= 4*(2-(1/shared(3))-shared(4));
    end 
    sharedA = repmat(shared',[1 Nfits]);
    phi0A=repmat(phi0',[1 Nfits]);
    initZA=repmat(zstart(1),[1 Nfits]);
    initPhaseA=repmat(initPhase(1),[1 Nfits]);
    
    tic % CPU fitter
    [P,CRLB, LogL] = CPUmleFit_LM_4Pi(single(imstack(:, :, :, :)),uint32(sharedA),100,single(IAB(:,:,:,:,3)),single(IAB(:,:,:,:,1)),single(IAB(:,:,:,:,2)),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
    tCPU=toc;
    disp(['GPU fitter speed' num2str(Nfits/tCPU)])
%     if length(zstart)>1
%         for i=2:length(zstart)
%             initZA=repmat(zstart(i),[1 Nfits]);
%             initPhaseA=repmat(initPhase(i),[1 Nfits]);
%             [Ph,CRLBh, LogLh] = CPUmleFit_LM_4Pi(single(imstack(:, :, :, :)),uint32(sharedA),100,single(IAB(:,:,:,:,3)),single(IAB(:,:,:,:,1)),single(IAB(:,:,:,:,2)),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
%             %             [Ph,CRLBh, LogLh] = GPUmleFit_LM_4Pi_cuda11p3_sm50(single(imstack(:, :, :, :)),uint32(sharedA),100,single(IAB),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
%             %         indbettero=LogLh<LogL;
%             indbetter=LogLh-LogL>1e-4; %copy only everything if LogLh increases by more than rounding error.
%             P(indbetter,:)=Ph(indbetter,:);
%             CRLB(indbetter,:)=CRLBh(indbetter,:);
%             LogL(indbetter)=LogLh(indbetter);
%         end
%     end
    tic % GPU fitter
    [P,CRLB, LogL] = GPUmleFit_LM_4Pi(single(imstack(:, :, :, :)),uint32(sharedA),100,single(IAB),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
    tGPU=toc;
    disp(['GPU fitter speed' num2str(Nfits/tGPU)])
    if length(zstart)>1
        for i=2:length(zstart)
            initZA=repmat(zstart(i),[1 Nfits]);
            initPhaseA=repmat(initPhase(i),[1 Nfits]);
            [Ph,CRLBh, LogLh] = GPUmleFit_LM_4Pi(single(imstack(:, :, :, :)),uint32(sharedA),100,single(IAB),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
            indbetter=LogLh-LogL>1e-4; %copy only everything if LogLh increases by more than rounding error.
            P(indbetter,:)=Ph(indbetter,:);
            CRLB(indbetter,:)=CRLBh(indbetter,:);
            LogL(indbetter)=LogLh(indbetter);
        end
    end
    
     X=P(:,1);
     Y=P(:,2);
     Z=P(:,3+offset);
     Znm = (Z-z0)*dz;
     if shared(3)==0
         N11=(P(:,3)+P(:,4));
         N12=(P(:,3+(offset-(1-shared(4))*4)/2)+P(:,3+(offset-(1-shared(4))*4)/2+1));
     elseif shared(3)==2
         N11=P(:,3);
         N12=P(:,3+(offset-(1-shared(4))*4)/2);
     else
         N11=P(:,3);
         N12=P(:,3);
     end
    BG=P(:,2+offset);
    phi =P(:,4+offset);
    phi = wrapTo2Pi(phi);
    z_phi = z_from_phi_YL(Z, phi, k, z0, dz);
    
    fitResults(:,:) =[X Y Znm z_phi]; 
    
    indgood = (z_phi-ground_truth.znm)<100;
    
    s_x_found1= (X(indgood)-coordinates(indgood,1,1));
    s_y_found1= (Y(indgood)-coordinates(indgood,2,1));
    s_z_found1=(Znm(indgood)-ground_truth.znm(indgood));
    s_phiZ_found1 = (z_phi(indgood)-ground_truth.znm(indgood));
 
    
 
 ratio=0.76;   
  [DL650_no, DL650_xo]=hist(DL650,1000);
 distribution1=zeros(1000,2);
distribution1(:,1) = DL650_xo';
distribution1(:,2) = DL650_no';
N = rand_arb_cjw(Nfits, distribution1, 0);
 ground_truth.N  = [N N N*ratio N*ratio];
 
 imstack = simSplinePSF(Npixels, PSF, ground_truth.N, ground_truth.bg, coordinates, phi0);%simulate images
 initZA=repmat(zstart(1),[1 Nfits]);
initPhaseA=repmat(initPhase(1),[1 Nfits]);
%     tic % CPU fitter
%     [P,CRLB, LogL] = CPUmleFit_LM_4Pi(single(imstack(:, :, :, :)),uint32(sharedA),100,single(IAB(:,:,:,:,3)),single(IAB(:,:,:,:,1)),single(IAB(:,:,:,:,2)),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
%     tCPU=toc;
%     disp(['GPU fitter speed' num2str(Nfits/tCPU)])
%     if length(zstart)>1
%         for i=2:length(zstart)
%             initZA=repmat(zstart(i),[1 Nfits]);
%             initPhaseA=repmat(initPhase(i),[1 Nfits]);
%             [Ph,CRLBh, LogLh] = CPUmleFit_LM_4Pi(single(imstack(:, :, :, :)),uint32(sharedA),100,single(IAB(:,:,:,:,3)),single(IAB(:,:,:,:,1)),single(IAB(:,:,:,:,2)),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
%             %             [Ph,CRLBh, LogLh] = GPUmleFit_LM_4Pi_cuda11p3_sm50(single(imstack(:, :, :, :)),uint32(sharedA),100,single(IAB),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
%             %         indbettero=LogLh<LogL;
%             indbetter=LogLh-LogL>1e-4; %copy only everything if LogLh increases by more than rounding error.
%             P(indbetter,:)=Ph(indbetter,:);
%             CRLB(indbetter,:)=CRLBh(indbetter,:);
%             LogL(indbetter)=LogLh(indbetter);
%         end
%     end
%     
    tic % GPU fitter
    [P,CRLB, LogL] = GPUmleFit_LM_4Pi(single(imstack(:, :, :, :)),uint32(sharedA),100,single(IAB),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
    tGPU=toc;
    disp(['GPU fitter speed' num2str(Nfits/tGPU)])
        if length(zstart)>1
        for i=2:length(zstart)
            initZA=repmat(zstart(i),[1 Nfits]);
            initPhaseA=repmat(initPhase(i),[1 Nfits]);
            [Ph,CRLBh, LogLh] = GPUmleFit_LM_4Pi(single(imstack(:, :, :, :)),uint32(sharedA),100,single(IAB),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
            indbetter=LogLh-LogL>1e-4; %copy only everything if LogLh increases by more than rounding error.
            P(indbetter,:)=Ph(indbetter,:);
            CRLB(indbetter,:)=CRLBh(indbetter,:);
            LogL(indbetter)=LogLh(indbetter);
        end
    end

    X=P(:,1);
    Y=P(:,2);
    Z=P(:,3+offset);
    Znm = (Z-z0)*dz;
     if shared(3)==0
         N21=(P(:,3)+P(:,4));
         N22=(P(:,3+(offset-(1-shared(4))*4)/2)+P(:,3+(offset-(1-shared(4))*4)/2+1));
     elseif shared(3)==2
         N21=P(:,3);
         N22=P(:,3+(offset-(1-shared(4))*4)/2);
     else
         N21=P(:,3);
         N22=P(:,3);
     end
    BG=P(:,2+offset);
    phi =P(:,4+offset);
    phi = wrapTo2Pi(phi);
    z_phi = z_from_phi_YL(Z, phi, k, z0, dz);
    
    fitResults(:,:) =[X Y Znm z_phi]; 
    
    indgood = (z_phi-ground_truth.znm)<100;
    
    s_x_found2= (X(indgood)-coordinates(indgood,1,1));
    s_y_found2= (Y(indgood)-coordinates(indgood,2,1));
    s_z_found2=(Znm(indgood)-ground_truth.znm(indgood));
    s_phiZ_found2 = (z_phi(indgood)-ground_truth.znm(indgood));
    
    
    
 ratio=0.4511;
 [DY634_no,DY634_xo]=hist(DY634,1000);
 distribution1=zeros(1000,2);
distribution1(:,1) = DY634_xo';
distribution1(:,2) = DY634_no';
N = rand_arb_cjw(Nfits, distribution1, 0);
 ground_truth.N  = [N N N*ratio N*ratio];

 imstack = simSplinePSF(Npixels, PSF, ground_truth.N, ground_truth.bg, coordinates, phi0);%simulate images
   initZA=repmat(zstart(1),[1 Nfits]);
initPhaseA=repmat(initPhase(1),[1 Nfits]);
%     tic % CPU fitter
%     [P,CRLB, LogL] = CPUmleFit_LM_4Pi(single(imstack(:, :, :, :)),uint32(sharedA),100,single(IAB(:,:,:,:,3)),single(IAB(:,:,:,:,1)),single(IAB(:,:,:,:,2)),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
%     tCPU=toc;
%     disp(['GPU fitter speed' num2str(Nfits/tCPU)])
%     if length(zstart)>1
%         for i=2:length(zstart)
%             initZA=repmat(zstart(i),[1 Nfits]);
%             initPhaseA=repmat(initPhase(i),[1 Nfits]);
%             [Ph,CRLBh, LogLh] = CPUmleFit_LM_4Pi(single(imstack(:, :, :, :)),uint32(sharedA),100,single(IAB(:,:,:,:,3)),single(IAB(:,:,:,:,1)),single(IAB(:,:,:,:,2)),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
%             %             [Ph,CRLBh, LogLh] = GPUmleFit_LM_4Pi_cuda11p3_sm50(single(imstack(:, :, :, :)),uint32(sharedA),100,single(IAB),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
%             %         indbettero=LogLh<LogL;
%             indbetter=LogLh-LogL>1e-4; %copy only everything if LogLh increases by more than rounding error.
%             P(indbetter,:)=Ph(indbetter,:);
%             CRLB(indbetter,:)=CRLBh(indbetter,:);
%             LogL(indbetter)=LogLh(indbetter);
%         end
%     end
%     
    tic % GPU fitter
    [P,CRLB, LogL] = GPUmleFit_LM_4Pi(single(imstack(:, :, :, :)),uint32(sharedA),100,single(IAB),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
    tGPU=toc;
    disp(['GPU fitter speed' num2str(Nfits/tGPU)])
        if length(zstart)>1
        for i=2:length(zstart)
            initZA=repmat(zstart(i),[1 Nfits]);
            initPhaseA=repmat(initPhase(i),[1 Nfits]);
            [Ph,CRLBh, LogLh] = GPUmleFit_LM_4Pi(single(imstack(:, :, :, :)),uint32(sharedA),100,single(IAB),single(dTAll),single(phi0A),single(initZA),single(initPhaseA));
            indbetter=LogLh-LogL>1e-4; %copy only everything if LogLh increases by more than rounding error.
            P(indbetter,:)=Ph(indbetter,:);
            CRLB(indbetter,:)=CRLBh(indbetter,:);
            LogL(indbetter)=LogLh(indbetter);
        end
    end
    X=P(:,1);
    Y=P(:,2);
    Z=P(:,3+offset);
    Znm = (Z-z0)*dz;
     if shared(3)==0
         N31=(P(:,3)+P(:,4));
         N32=(P(:,3+(offset-(1-shared(4))*4)/2)+P(:,3+(offset-(1-shared(4))*4)/2+1));
     elseif shared(3)==2
         N31=P(:,3);
         N32=P(:,3+(offset-(1-shared(4))*4)/2);
     else
         N31=P(:,3);
         N32=P(:,3);
     end
    BG=P(:,2+offset);
    phi =P(:,4+offset);
    phi = wrapTo2Pi(phi);
    z_phi = z_from_phi_YL(Z, phi, k, z0, dz);
    
    fitResults(:,:) =[X Y Znm z_phi]; 
    
    indgood = (z_phi-ground_truth.znm)<100;
    
    s_x_found3= (X(indgood)-coordinates(indgood,1,1));
    s_y_found3= (Y(indgood)-coordinates(indgood,2,1));
    s_z_found3=(Znm(indgood)-ground_truth.znm(indgood));
    s_phiZ_found3 = (z_phi(indgood)-ground_truth.znm(indgood));
    

NF=1:Nfits;

figure,h1=plot(N11,N12,'r.');
hold on
h2=plot(N21,N22,'b.');
h3=plot(N31,N32,'g.');
xlabel('Number of Photos without filter')
ylabel('Number of Photos with filter')
legend([h1(1) h2(1) h3(1)],'CF680','DL649','DY633')

maxl1=4.3;
minl1=2.6;
maxl2=4.3;
minl2=2.6;
figure,h1=scatter(N11,N12,1,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor','r');
hold on
h2=scatter(N21,N22,1,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor',[0.5 0.5 0.5]);
h3=scatter(N31,N32,1,'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor',[0.667 0.667 0.1]);
% axis([0 2500 0 2500])
ax=gca;
% axis(ax,'equal')
ax.XScale='log';
ax.YScale='log';
xlim(10.^[minl1 maxl1])
ylim(10.^[minl2 maxl2])
x_label = xlabel('Photon number of {\itf}^1 channel');
y_label = ylabel('Photon number of {\itf}^2 channel');
set(gca,'FontName','time','FontSize',15)
% legend([h1(1) h2(1) h3(1)],'CF680','DL649','DY633',2)
lgd = legend('CF680','DL649','DY633','Location','northwest');
set(lgd,'FontName','time','FontSize',15,'LineWidth',1,'FontWeight','bold');   
set(x_label,'FontName','time','FontSize',20,'LineWidth',3,'FontWeight','bold');
set(y_label,'FontName','time','FontSize',20,'LineWidth',3,'FontWeight','bold');
set(gca,'FontName','time','FontSize',15,'FontWeight','bold');
legend('boxoff')


figure,plot(NF(1:length(s_phiZ_found1)),s_phiZ_found1,'ro');
hold on
plot(NF(1:length(s_phiZ_found2)),s_phiZ_found2,'bo');
plot(NF(1:length(s_phiZ_found3)),s_phiZ_found3,'go');
xlabel('molecus')
ylabel('localization error (nm)')
legend('locaization precision Z CF680','locaization precision Z DL649','locaization precision Z DY633')


figure,plot(NF(1:length(s_x_found1)),s_x_found1*120,'ro');
hold on
plot(NF(1:length(s_x_found2)),s_x_found2*120,'bo');
plot(NF(1:length(s_x_found3)),s_x_found3*120,'go');
xlabel('molecus')
ylabel('localization error (nm)')
legend('locaization precision Z CF680','locaization precision Z DL649','locaization precision Z DY633')

figure,plot(NF(1:length(s_y_found1)),s_y_found1*120,'rs');
hold on
plot(NF(1:length(s_x_found2)),s_y_found2*120,'bo');
plot(NF(1:length(s_x_found3)),s_y_found3*120,'go');
xlabel('molecus')
ylabel('localization error (nm)')
legend('locaization precision Z CF680','locaization precision Z DL649','locaization precision Z DY633')


