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
 
MAX_PATH_COUNT = 10;
CHANNEL_COUNT = 50;

SR = 15.36e6;
T = SR * 1e-3;

Nt = 1;
Nr = M;

cdlList = randCdlProfile(CHANNEL_COUNT, MAX_PATH_COUNT, M, fc, fd, SR);

x = ["Image, AOA, ZOA, Delays, Gains, NumPaths"];

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
    noise = db2mag(-65) .* complex(randn(size(rxWaveform)),randn(size(rxWaveform)));
    rxWaveform = rxWaveform + noise;

    rx = nrOFDMDemodulate(c,rxWaveform);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.
    % YOLO PAPER IMAGE GENERATION

    % Oversampling factors
    alpha = 12;
    beta = 1;

    % Max after normalization
    delta = 2^12;

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
%     figure()
% 
%     image(Y_tilda);
%     title(i)

%     imwrite(repmat(uint8(Y_tilda), [1 1 3]), 'example.png');
    imwrite(uint8(Y_tilda), sprintf('ChannelImages\\%d.png',i));
    
    x = [x sprintf('%d, %s, %s, %s, %s, %d',i, num2str(chInfo.AnglesAoA), num2str(chInfo.AnglesZoA), num2str(chInfo.PathDelays), num2str(chInfo.AveragePathGains), size(chInfo.AnglesZoA, 2))];
    
end


fid = fopen( outputFile, 'w' );
for i = 1 : length(x)
    fprintf(fid, x(i) + "\n");
end

fclose(fid);

function cdl = randCdlProfile(n, LMax, M, fc, fd, SR)
    cdl = cell(1,n);
    for i = 1:n
        % Add 1/2 before rounding to int causing it to always round up.
        % [0,1) -> 1, [1,2) -> 2, ..., [L-1, L) -> L
        L = int64((rand(1,1)*LMax)+0.5);
        cdl{i} = nrCDLChannel;

        cdl{i}.CarrierFrequency = fc;
        cdl{i}.MaximumDopplerShift = fd;
        cdl{i}.SampleRate = SR;

        cdl{i}.TransmitAntennaArray.Size = [1 1 1 1 1];
        cdl{i}.ReceiveAntennaArray.Size = [1 M 1 1 1];
        
        cdl{i}.DelayProfile = 'custom';
        cdl{i}.PathDelays = rand(1,L) .* 10^-4;
    %     cdl.AveragePathGains = starterCDL.AveragePathGains + ones(size(starterCDL.AveragePathGains)); % [1 1 1 1];
        cdl{i}.AveragePathGains = rand(1,L) *10; % [1 1 1 1];
        cdl{i}.AnglesAoA = (rand(1,L)-0.5)*360;
        cdl{i}.AnglesZoA = (rand(1,L)-0.5)*360;
        cdl{i}.AnglesAoD = (rand(1,L)-0.5)*360; 
        cdl{i}.AnglesZoD = (rand(1,L)-0.5)*360; 
        cdl{i}.HasLOSCluster = false; % default
        % cdl.KFactorFirstCluster = 13.3; % default
        cdl{i}.AngleSpreads = rand(1,4)*2; % [ASD ASA ZSD ZSA] 
        cdl{i}.XPR = 0; % default, crosspolarization power in db
        cdl{i}.NumStrongestClusters = 0; % default
        
    end
    
    
end
