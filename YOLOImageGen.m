%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              YOLOImageGen.m                             %
%                             ----------------                            %
%   Script attempting to generate the image used by YOLO following steps  %
% from their paper.                                                       %
%                                                                         %
% By: Ben Earle (BenEarle@cmail.carleton.ca)                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFN FROM PAPER ON NOMP
c = physconst('lightspeed'); % speed of light in m/s
fftPoints = 2048;
df = 75 * 10^3; %hz
BW = 90 * 10^6; %hz
fc = 3.5 * 10^9; %hz
pt = -20;% dbm (transmit power)
lambda = c / fc; % carrier wavelength in m
d = lambda / 2; % distance between antennas in m

N = 32; % Sub-carrier count
M = 32; % Antenna count

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple ULA transmission

v = 15.0;                    % UE velocity in km/h
%fc = 4e9;                    % carrier frequency in Hz

fd = (v*1000/3600)/c*fc;     % UE max Doppler frequency in Hz
 
cdl = nrCDLChannel;
cdl.DelayProfile = 'CDL-A';
cdl.DelaySpread = 10e-9;
% cdl.DelaySpread = 0;
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

txWaveform = ones(N,Nt);

chInfo = info(cdl);

% The first time we pass the signal through the CDL channel the first 7 
% subcarriers receive zero value on each antenna, I don't really understand
% why that is... Running it the second time does not have this issue.
[rxWaveform, pathGains] = cdl(txWaveform);
[rxWaveform, pathGains] = cdl(txWaveform);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.
% YOLO PAPER IMAGE GENERATION

% Oversampling factors
alpha = 16;
beta = 16;

% Max after normalization
delta = 255;

% Received waveform Y 
Y = rxWaveform.';

% DFT Matricies
U_theta = dftmtx(alpha * M);
U_theta = U_theta(1:M,:);
U_T = dftmtx(beta * N);
U_T = U_T(1:N,:);

Y_bar = U_theta.' * Y * U_T;

Y_tilda = (delta/max(max(abs(Y_bar)))) .* abs(Y_bar);

image(Y_tilda);


