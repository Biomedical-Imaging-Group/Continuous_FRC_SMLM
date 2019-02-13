%-----------------------------------------------------------
% Continuous Fourier Ring Correlation
%
% Exemple of use of the function FRCcont with and without the use of the
% lookup table. 
%
% Copyright (2018) Thanh-An Pham (thanh-an.pham@epfl.ch)
%                  Emmanuel Soubies (esoubies@gmail.com)
%-----------------------------------------------------------
clear, close all;

%% Parameters
fov=6400;         % Field of view to generate points (nm)     
Nmol=1000;        % Number of molecules
stdMol=8;         % Standard deviation between the position of the molecules of the two sets
Nsamples = 200;   % Number of points in the FRC curve
fmax = 0.1;       % Maximal frequency. The FRC will be computed on [0,fmax]

%% Load data
gt=rand(Nmol,2)*fov;
loc=gt+randn(Nmol,2)*stdMol;

%% FRC computation
tic;[FRCcLook,freq] = FRCcont(loc,gt,fov,Nsamples,fmax,1);
disp(['Computation with Lookup table : ',num2str(toc),' s']);
tic;[FRCcNoLook,freq] = FRCcont(loc,gt,fov,Nsamples,fmax,0);
disp(['Computation without Lookup table : ',num2str(toc),' s']);

%% Plots
figure;hold all;
plot(freq,FRCcLook,'linewidth',1.5);
plot(freq,FRCcNoLook,'linewidth',1.5,'linestyle','--');grid
xlabel('f'); ylabel('FRC'); set(gca,'FontSize',14);
legend('Lookup Table','No Lookup Table');axis([0 0.1 0 1]);

