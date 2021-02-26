%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               SimpleULA.m                               %
%                             ---------------                             %
%   Script attempting to recreate NOMP channel reconstruction. I want to  %
% compare the channel generated using the path components to the actual   %
% channel used by 5G toolbox.                                             %
%                                                                         %
% By: Ben Earle (BenEarle@cmail.carleton.ca)                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFN FROM PAPER ON NOMP
fftPoints = 2048;
df = 75 * 10^3; %hz
BW = 90 * 10^6; %hz
fc =3.5 * 10^9; %hz
pt = -20;% dbm (transmit power)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple ULA transmission

v = 15.0;                    % UE velocity in km/h
%fc = 4e9;                    % carrier frequency in Hz
c = physconst('lightspeed'); % speed of light in m/s
fd = (v*1000/3600)/c*fc;     % UE max Doppler frequency in Hz
 
cdl = nrCDLChannel;
cdl.DelayProfile = 'CDL-D';
cdl.DelaySpread = 10e-9;
cdl.CarrierFrequency = fc;
cdl.MaximumDopplerShift = fd;

cdl.TransmitAntennaArray.Size = [1 1 1 1 1];
cdl.ReceiveAntennaArray.Size = [64 1 1 1 1];

SR = 15.36e6;
T = SR * 1e-3;
cdl.SampleRate = SR;
cdlinfo = info(cdl);
Nt = cdlinfo.NumTransmitAntennas;
 
txWaveform = ones(T,Nt);

chInfo = info(cdl);

rxWaveform = cdl(txWaveform);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General defs
N = 52; % Sub-carrier count
M = 64; % Antenna count
df = 15000 % Distance between subcarriers
% Get the set of path components from the channel model.
L = length(chInfo.PathDelays);
tau = chInfo.PathDelays;
% I do not know what anlge I should be using... hmmm
% Likely CDL channel definition with antennas along x axis means I should
% use a specific one, but which one? IDK yet
theta = chInfo.AnglesAoA;
gul = chInfo.AveragePathGains;

% Calculate the stacked vectors


% a(theta)
p(tau(10), N, df)



%hul = 





function out = a(theta, M, dF)
    m = -M/2:(M/2-1);
    
    
    out = 1;
end
function out = p(tau, N, df)
    n = -N/2:(N/2-1);
    out = exp(j*2*pi*n*df*tau)';
end