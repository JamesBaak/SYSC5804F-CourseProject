%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             ChRecFromPaths.m                            %
%                            ------------------                           %
%   Script attempting to recreate NOMP channel reconstruction. I want to  %
% compare the channel generated using the path components to the actual   %
% channel used by 5G toolbox.                                             %
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple ULA transmission

v = 15.0;                    % UE velocity in km/h
%fc = 4e9;                    % carrier frequency in Hz

fd = (v*1000/3600)/c*fc;     % UE max Doppler frequency in Hz
 
cdl = nrCDLChannel;
cdl.DelayProfile = 'CDL-D';
% cdl.DelaySpread = 10e-9;
cdl.DelaySpread = 0;
cdl.CarrierFrequency = fc;
% cdl.MaximumDopplerShift = fd;
cdl.MaximumDopplerShift = 0;

cdl.TransmitAntennaArray.Size = [1 1 1 1 1];
cdl.ReceiveAntennaArray.Size = [1 8 1 1 1];

SR = 15.36e6;
T = SR * 1e-3;
cdl.SampleRate = SR;
cdlinfo = info(cdl);
Nt = cdlinfo.NumTransmitAntennas;
Nr = cdlinfo.NumReceiveAntennas;
 
txWaveform = ones(52,Nt);

chInfo = info(cdl);

% The first time we pass the signal through the CDL channel the first 7 
% subcarriers receive zero value on each antenna, I don't really understand
% why that is... Running it the second time does not have this issue.
[rxWaveform, pathGains] = cdl(txWaveform);
[rxWaveform, pathGains] = cdl(txWaveform);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel Reconstruction
N = 52; % Sub-carrier count
M = Nr; % Antenna count

% Get the set of path components from the channel model.
L = length(chInfo.PathDelays);
tau = chInfo.PathDelays;
% I do not know what anlge I should be using... hmmm
% Likely CDL channel definition with antennas along x axis means I should
% use a specific one, but which one? IDK yet
theta = wrapTo180(chInfo.AnglesAoA - 180 .* ones(size(chInfo.AnglesAoA)));
gul = squeeze(pathGains(1,:,1,:));

% Reconstruct the channel using method from NOMP paper


% Temporarily initialize the h to the exact size I am getting from kron
h = zeros(N, M);
for i = 1 : length(theta)
     h = h + gul(i,:) .* kron(p(tau(i), N, df), a(theta(i), M, d, lambda)');
%    h = h + gul(i) .* a(theta(i), M, d, lambda) * p(tau(i), N, df)';
end

disp("RMSE: ");
mean(mean(abs((h-rxWaveform).^2)))^1/2
% k = kron(p(tau(i), N, df), a(theta(i), M, d, lambda));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function defs
% a() returns the steering vector for a given angle
% Input:
%    theta - angle of the received path
%    M     - antenna count
%    d     - distance between antennas
%    labda - wavelength of the carrier
% Output:
%    out   - 

function out = a(theta, M, d, lambda)
    m = -M/2:(M/2-1);
    out = exp(1j*2*pi*m*d*lambda*sin(theta))';
end

% p() returns the steering vector for a given angle
% Input:
%    tau   - received path
%    M     - antenna count
%    df    - subcarrier spacing
% Output:
%    out   - steering vector of the ULA
function out = p(tau, N, df)
    n = -N/2:(N/2-1);
    out = exp(1j*2*pi*n*df*tau)';
end
