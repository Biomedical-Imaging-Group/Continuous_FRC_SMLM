function [FRC,freq] = FRCcont(S1,S2,distMax,Nsamples,fmax,useLookUp)
%-----------------------------------------------------------
% [FRC,freq] = FRCcont(S1,S2,distMax,Nsamples,fmax,useLookUp)
%
% Compute the FRC in continuous domain [1]
%
% Input:  S1        -> First set of points of size (N1 x dim)
%         S2        -> Second set of points of size (N2 x dim)
%         distMax   -> Maximal distance between two points (used if useLookUp=true)
%         Nsamples  -> Number of samples for the FRC curve (default 1000)
%         fmax      -> Maxiamal frequency for the FRC curve (default 0.1)
%         useLookUp -> Boolean true to use a lookup table for Besselj (default false)
%
% Output: FRC  -> FRC curve
%         freq -> frequencies where the FRC is evaluated 
%
% Reference:
% [1] Closed-Form Expression of the Fourier Ring-Correlation for Single-Molecule Localization Microscopy. 
%     Proceeding of ISBI, 2019.
%     Thanh-an Pham, Emmanuel Soubies, Daniel Sage, and Michael Unser.
%
% Copyright (2018) Thanh-An Pham (thanh-an.pham@epfl.ch)
%                  Emmanuel Soubies (esoubies@gmail.com)
%-----------------------------------------------------------

%% Set default parameters
if nargin < 4 || isempty(Nsamples), Nsamples = 1000; end
if nargin < 5 || isempty(fmax), fmax = 0.1; end
if nargin < 6 || isempty(useLookUp), useLookUp = false; end

%% Compute Distances between each pair of points
D = size(S1,2);
dist1 = 0;
dist2 = 0;
dist = 0;
for d=1:D
    dist1 = dist1 + (S1(:,d) - S1(:,d)').^2;
    dist2 = dist2 + (S2(:,d) - S2(:,d)').^2;
    dist = dist + (S1(:,d) - S2(:,d)').^2;
end
dist1 = single(sqrt(dist1));
dist2 = single(sqrt(dist2));
dist = single(sqrt(dist));

%% Compute FRC
% - rho discretization
k = 2*pi*fmax;
rho_set_orig = linspace(0,k, Nsamples);
FRC = zeros(Nsamples,1);
% - Lookup Table
step = 0.001;
xquery = (0:step:distMax)';
lookupB = besselj(0,xquery);
max_len = floor(2e9/(4*max(numel(dist),max(numel(dist1),numel(dist2)))));
% - Main Loop
for kk = 1:ceil(Nsamples/max_len)
    cind = 1 + (kk - 1)*max_len:min(kk*max_len,Nsamples);
    rho_set = rho_set_orig(cind);
    if useLookUp
        J0r = sum(sum(mybesselj(reshape(rho_set,[1,1,length(rho_set)]).*dist),1),2);
        S1N = sum(sum(mybesselj(reshape(rho_set,[1,1,length(rho_set)]).*dist1),1),2);
        S2N = sum(sum(mybesselj(reshape(rho_set,[1,1,length(rho_set)]).*dist2),1),2);
    else
        J0r = sum(sum(besselj(0,reshape(rho_set,[1,1,length(rho_set)]).*dist),1),2);
        S1N = sum(sum(besselj(0,reshape(rho_set,[1,1,length(rho_set)]).*dist1),1),2);
        S2N = sum(sum(besselj(0,reshape(rho_set,[1,1,length(rho_set)]).*dist2),1),2);
    end
    FRCdenom = sqrt(S1N.*S2N);
    FRC(cind) = squeeze(J0r./FRCdenom);
end
freq=linspace(0,fmax,Nsamples);

%% Evaluate Bessel with a lookup table
    function out = mybesselj(x)
        sz_in = size(x);
        try
            flx=floor(x(:)/step);
            out = reshape(lookupB(1 + flx),sz_in); % Faster but noisy
        catch
            fprintf('Increase the maximal distance for the lookup table as the lookup table limit is reached\n');
        end
    end
end