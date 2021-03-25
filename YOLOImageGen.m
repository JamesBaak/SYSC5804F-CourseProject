%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              YOLOImageGen.m                             %
%                             ----------------                            %
%   Script attempting to generate the image used by YOLO following steps  %
% from their paper.                                                       %
%                                                                         %
% By: Ben Earle (BenEarle@cmail.carleton.ca)                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
outputFile = "pathData.csv";
% outputDir = "C:\Users\BenEa\Desktop\Datasets\ChannelImages\";

%DEFN FROM PAPER ON NOMP
c = physconst('lightspeed'); % speed of light in m/s
fftPoints = 2048;
df = 75000; %hz
BW = 90 * 7^6; %hz
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
 
MAX_PATH_COUNT = 10;
CHANNEL_COUNT = 500;

SR = 15.36e6;
T = SR * 1e-3;

Nt = 1;
Nr = M;

cdlList = randCdlProfile(CHANNEL_COUNT, MAX_PATH_COUNT, M, fc, fd, SR);

AoA_str = "";
AoZ_str = "";
T_str = "";
g_str = "";

for i = 1:MAX_PATH_COUNT 
    AoA_str = strcat(AoA_str, "AoA_"+i, ",");
    AoZ_str = strcat(AoZ_str, "AoZ_"+i, ",");
    T_str   = strcat(T_str, "T_"+i, ",");
    g_str   = strcat(g_str, "g_"+i, ",");
end


x = [strcat("Index,", AoA_str, AoZ_str, T_str, g_str, "NumPaths,MaxPixel")];

for i = 1:CHANNEL_COUNT
    cdl = cdlList{i};
    c = nrCarrierConfig;
    c.NSizeGrid = N;
    c.SubcarrierSpacing = 15;
    c.NSlot = 1;
    c.NFrame = 0;

    txOne = ones(N*12,Nt);
    txZero= zeros(N,Nt);
 
    txWaveform = nrOFDMModulate(c,txOne);
    % maxChDelay = ceil(max(chInfo.PathDelays*cdl.SampleRate)) + chInfo.ChannelFilterDelay;

    chInfo = info(cdl);


    % The first time we pass the signal through the CDL channel the first 7 
    % subcarriers receive zero value on each antenna, I don't really understand
    % why that is... Running it the second time does not have this issue.
    [rxWaveform, pathGains] = cdl(txWaveform);
    [rxWaveform, pathGains] = cdl(txWaveform);
    noise = db2mag(-80) .* complex(randn(size(rxWaveform)),randn(size(rxWaveform)));
%     rxWaveform = rxWaveform + noise;

    rx = nrOFDMDemodulate(c,rxWaveform);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.
    % YOLO PAPER IMAGE GENERATION

    % Oversampling factors
    alpha = 12;
    beta = 1;

    % Max after normalization
    delta = 2^10;

    % Received waveform Y 
    Y =  squeeze(rx).';

    % DFT Matricies
    U_theta = dftmtx(alpha * M);
    U_theta = U_theta(1:M,:);
    U_T = dftmtx(beta * N*12);
    U_T = U_T(1:N*12,:);

    Y_bar = U_theta.' * Y * U_T;
    
    [max_Y_bar,y_coord] = max(abs(Y_bar));
    [max_Y_bar,x_coord] = max(max_Y_bar);

    Y_tilda = (delta/max_Y_bar) .* abs(Y_bar);
    % Y_tilda = abs(Y_bar) .* 50;
%     figure()
% 
%     image(Y_tilda);
%     title(i)

%     imwrite(repmat(uint8(Y_tilda), [1 1 3]), 'example.png');
    imwrite(uint8(Y_tilda), sprintf('ChannelImages\\%d.png',i));
    
    delay_pad = padarray(chInfo.PathDelays, [0,MAX_PATH_COUNT-size(chInfo.PathDelays, 2)],0,"post");
    AoA_pad   = padarray(chInfo.AnglesAoA, [0,MAX_PATH_COUNT-size(chInfo.AnglesAoA, 2)],0,"post");
    ZoA_pad   = padarray(chInfo.AnglesZoA, [0,MAX_PATH_COUNT-size(chInfo.AnglesZoA, 2)],0,"post");
    gain_pad  = padarray(chInfo.AveragePathGains, [0,MAX_PATH_COUNT-size(chInfo.AveragePathGains, 2)],-1000,"post");
    
    x = [x strcat(sprintf('%d,', i), sprintf('%f,' , AoA_pad), sprintf('%f,' , ZoA_pad), sprintf('%f,' , delay_pad.*1000), sprintf('%f,' , gain_pad), sprintf('%d, %d %d',  size(chInfo.AnglesZoA, 2), x_coord, y_coord(x_coord)))];
    
    % Test dominant path recovery
%     [height, width] = size(Y_tilda);
%     THETA_hat = 1-(y_coord(x_coord)/height);
%     T_hat = (x_coord/width);
    
%     Y(X)
%     X
%      theta_hat = lambda / d * asind(THETA_hat)
%      t_hat = T_hat / df
% 
%     theta = lambda / d * asind(THETA_hat);
%     tau = T_hat / df;
%     gul = squeeze(pathGains(1,:,1,:));

%     h = zeros(N, M);
%     for i = 1 : length(theta)
%          h = h + gul(i,:) .* kron(p(tau(i), N, df), a(theta(i), M, d, lambda)');
%     end
    
%     h = gul .* kron(p(tau, N*12, df), a(theta, M, d, lambda)');
%     mean(mean((h'-Y)))
%     mean(mean((Y)))
%     mean(mean((h')))
%     size(chInfo.AnglesZoA, 2)
end


fid = fopen( outputFile, 'w' );
for i = 1 : length(x)
    fprintf(fid, x(i) + "\n");
end

fclose(fid);


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

function cdl = randCdlProfile(n, LMax, M, fc, fd, SR)
    cdl = cell(1,n);
    for i = 1:n
        % Add 1/2 before rounding to int causing it to always round up.
        % [0,1) -> 1, [1,2) -> 2, ..., [L-1, L) -> L
%         L = int64((rand(1,1)*LMax)+0.5);
        L = 5;
        cdl{i} = nrCDLChannel;

        cdl{i}.CarrierFrequency = fc;
        cdl{i}.MaximumDopplerShift = fd;
        cdl{i}.SampleRate = SR;

        cdl{i}.TransmitAntennaArray.Size = [1 1 1 1 1];
        cdl{i}.ReceiveAntennaArray.Size = [1 M 1 1 1];
        
        cdl{i}.DelayProfile = 'custom';
        cdl{i}.PathDelays = rand(1,L) .* 10^-4;
    %     cdl.AveragePathGains = starterCDL.AveragePathGains + ones(size(starterCDL.AveragePathGains)); % [1 1 1 1];
        cdl{i}.AveragePathGains = (rand(1,L)-0.5) * 30; % [1 1 1 1];
        cdl{i}.AnglesAoA = (rand(1,L)-0.5)*360;
        cdl{i}.AnglesZoA = (rand(1,L)-0.5)*360;
%         cdl{i}.AnglesZoA = (ones(1,L))*90;
        cdl{i}.AnglesAoD = (rand(1,L)-0.5)*360; 
        cdl{i}.AnglesZoD = (rand(1,L)-0.5)*360; 
        cdl{i}.HasLOSCluster = false; % default
        % cdl.KFactorFirstCluster = 13.3; % default
        cdl{i}.AngleSpreads = [0 0 0 0]; %rand(1,4)*0; % [ASD ASA ZSD ZSA] 
        cdl{i}.XPR = 0; % default, crosspolarization power in db
        cdl{i}.NumStrongestClusters = 0; % default
        
    end
    
    
end
