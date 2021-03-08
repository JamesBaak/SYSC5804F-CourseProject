%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              YOLOImageGen.m                             %
%                             ----------------                            %
%   Script attempting to generate the image used by YOLO following steps  %
% from their paper.                                                       %
%                                                                         %
% By: Ben Earle (BenEarle@cmail.carleton.ca)                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

%DEFN FROM PAPER ON NOMP
c = physconst('lightspeed'); % speed of light in m/s
fftPoints = 2048;
df = 75; %khz
BW = 90 * 10^6; %hz
fc = 3.5 * 10^9; %hz
pt = -20;% dbm (transmit power)
lambda = c / fc; % carrier wavelength in m
d = lambda / 2; % distance between antennas in m

N = 128; % Sub-carrier count
M = 128; % Antenna count

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple ULA transmission

v = 15.0;                    % UE velocity in km/h
%fc = 4e9;                    % carrier frequency in Hz

fd = (v*1000/3600)/c*fc;     % UE max Doppler frequency in Hz
 
CUSTOM = 1;

cdl = nrCDLChannel;

x = nrCDLChannel;
x.DelayProfile = 'CDL-E';
starterCDL = info(x);

if CUSTOM
%     cdl.DelayProfile = 'custom';
%     cdl.PathDelays = starterCDL.PathDelays; % [0 0.022 0.214 0.95] .* 1.0e-03;
% %     cdl.AveragePathGains = starterCDL.AveragePathGains + ones(size(starterCDL.AveragePathGains)); % [1 1 1 1];
%     cdl.AveragePathGains = ones(size(starterCDL.AveragePathGains)); % [1 1 1 1];
%     cdl.AnglesAoA = starterCDL.AnglesAoA; % [12 36 -27 -123];
%     cdl.AnglesZoA = starterCDL.AnglesZoA; % [12 -36 27 123];
%     cdl.AnglesAoD = starterCDL.AnglesAoD; % [12 36 27 123];
%     cdl.AnglesZoD = starterCDL.AnglesZoD; % [-12 36 27 -123];
%     cdl.HasLOSCluster = false; % default
%     % cdl.KFactorFirstCluster = 13.3; % default
%     cdl.AngleSpreads = [90.0 90.0 40.0 13.0]; % [ASD ASA ZSD ZSA] 
%     cdl.XPR = 10; % default, crosspolarization power in db
%     cdl.NumStrongestClusters = 0; % default
    
    
    cdl.DelayProfile = 'custom';
    cdl.PathDelays = [0 5 10 12] .* 10^-5;
%     cdl.AveragePathGains = starterCDL.AveragePathGains + ones(size(starterCDL.AveragePathGains)); % [1 1 1 1];
    cdl.AveragePathGains = [15 -1 -12 -5]; % [1 1 1 1];
    cdl.AnglesAoA = [15 45 45 -45]; % [12 36 -27 -123];
    cdl.AnglesZoA = [15 45 90 -45]; % [12 -36 27 123];
    cdl.AnglesAoD = [0 0 0 0]; % [12 36 27 123];
    cdl.AnglesZoD = [0 0 0 0]; % [-12 36 27 -123];
    cdl.HasLOSCluster = false; % default
    % cdl.KFactorFirstCluster = 13.3; % default
    cdl.AngleSpreads = [0 0 0.0 0]; % [ASD ASA ZSD ZSA] 
    cdl.XPR = 0; % default, crosspolarization power in db
    cdl.NumStrongestClusters = 0; % default
else
    cdl.DelayProfile = 'CDL-A';
    cdl.DelaySpread = 10e-9;
    % cdl.DelaySpread = 0;
    
end

% cdl.MaximumDopplerShift = 0;
cdl.CarrierFrequency = fc;
cdl.MaximumDopplerShift = fd;


cdl.TransmitAntennaArray.Size = [1 1 1 1 1];
cdl.ReceiveAntennaArray.Size = [1 M 1 1 1];

SR = 15.36e6;
T = SR * 1e-3;
cdl.SampleRate = SR;
cdlinfo = info(cdl);
Nt = cdlinfo.NumTransmitAntennas;
Nr = cdlinfo.NumReceiveAntennas;


c = nrCarrierConfig;
c.NSizeGrid = N;
c.SubcarrierSpacing = 15;
c.NSlot = 1;
c.NFrame = 0;

txOne= ones(N*12,Nt);
txZero= zeros(N,Nt);

txWaveform = nrOFDMModulate(c,txOne);
% maxChDelay = ceil(max(chInfo.PathDelays*cdl.SampleRate)) + chInfo.ChannelFilterDelay;

chInfo = info(cdl);


% The first time we pass the signal through the CDL channel the first 7 
% subcarriers receive zero value on each antenna, I don't really understand
% why that is... Running it the second time does not have this issue.
[rxWaveform, pathGains] = cdl(txWaveform);

rx = nrOFDMDemodulate(c,rxWaveform);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.
% YOLO PAPER IMAGE GENERATION

% Oversampling factors
alpha = 16;
beta = 1;

% Max after normalization
delta = 10000;

% Received waveform Y 
Y =  squeeze(rx).';

% DFT Matricies
U_theta = dftmtx(alpha * M);
U_theta = U_theta(1:M,:);
U_T = dftmtx(beta * N*12);
U_T = U_T(1:N*12,:);

Y_bar = U_theta.' * Y * U_T;

Y_tilda = (delta/max(max(abs(Y_bar)))) .* abs(Y_bar);
% Y_tilda = abs(Y_bar) .* 50;
figure()

image(Y_tilda);
title(cdl.DelayProfile)

