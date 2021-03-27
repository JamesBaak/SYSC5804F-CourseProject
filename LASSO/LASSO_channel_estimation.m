% Implementing Joint Burst LASSO for Spare Channel Estimation
% [1] An Liu1, Vincent Lau1, and Wei Dai. "Joint Burst LASSO for Sparse Channel Estimation
% in Multi-user Massive MIMO". Wireless Communications Symposium. IEEE ICC 2016
% Using Clustered Delay Line (CDL) channel model with a CDL-D delay profile

% Constants ===============================================
v = 15.0;                    % UE velocity in km/h
fc = 4e9;                    % carrier frequency in Hz
c = physconst('lightspeed'); % speed of light in m/s
fd = (v*1000/3600)/c*fc;     % UE max Doppler frequency in Hz
subcarriers = 52;
antennas = 64;
nPilots = 50;
% =========================================================
range = 25:50;

for k = range
    nPilots = k;
    cdl = nrCDLChannel;
    cdl.DelayProfile = 'CDL-A';
    cdl.DelaySpread = 10e-9;
    cdl.CarrierFrequency = fc;
    cdl.MaximumDopplerShift = fd;

    cdl.TransmitAntennaArray.Size = [1 antennas 1 1 1];
    cdl.ReceiveAntennaArray.Size = [1 1 1 1 1];

    SR = 15.36e6;
    T = SR * 1e-3;
    cdl.SampleRate = SR;
    cdlinfo = info(cdl);
    Nt = cdlinfo.NumTransmitAntennas;
    Nr = cdlinfo.NumReceiveAntennas;

    % Generate waveform with pilots
    %refSig = zeros(subcarriers,Nt);
    %refSig(1:nPilots) = 1;
    %txWaveform = complex(refSig);
    txWaveform = complex(ones(nPilots,Nt));
    %txWaveform = complex(repmat([1:nPilots].',1,Nt));

    % Calculate H
    % Get the Inverse Discrete Fourier Transform of the tx waveform
    %ift = dsp.IFFT;
    %txWaveformA = ift(txWaveform);
    %txWaveformA = txWaveform;
    A = (1/sqrt(Nt)) * exp(-1j*2*pi/Nt .* ((0:Nt-1)' * (0:Nt-1)));
    txWaveformA = txWaveform * A';

    % Send Signal through the channel
    chInfo = info(cdl);
    [rxWaveform, pathGains, sampleTimes] = cdl(txWaveformA);

    % Lasso in attempt to find H
    % Complex LASSO code taken from Mark Schmidt ==============================
    % Link: https://www.cs.ubc.ca/~schmidtm/Software/code.html
    % Efficient projection onto this complex L1-Ball is described in:
    % van den Berg and Friedlander.  <http://www.optimization-online.org/DB_FILE/2008/01/1889.pdf Probing the Pareto Frontier for Basis
    % Pursuit Solutions>.  SIAM Journal of Scientific Computing (2008).

    % Functions for switching between the two complex representations
    makeReal = @(z)[real(z);imag(z)];
    makeComplex = @(zRealImag)zRealImag(1:Nt) + 1i*zRealImag(Nt+1:end);

    % Initial guess of parameters
    zRealImag = makeReal(zeros(Nt,1));

    % Set up Objective Function
    XRealImag = [real(txWaveformA) -imag(txWaveformA);imag(txWaveformA) real(txWaveformA)];
    yRealImag = [real(rxWaveform);imag(rxWaveform)];
    funObj = @(zRealImag)SquaredError(zRealImag,XRealImag,yRealImag);

    % Set up Complex L1-Ball Projection
    tau = 1;
    funProj = @(zRealImag)complexProject(zRealImag,tau);

    % Solve with PQN
    fprintf('\nComputing optimal Lasso parameters...\n');
    zRealImag = minConf_PQN(funObj,zRealImag,funProj);
    h_1 = makeComplex(zRealImag);
    % =========================================================================

    % Another method to complex LASSO regression ==============================
    % https://stats.stackexchange.com/questions/469653/implementing-complex-lasso-in-matlab
    % lasso(txWaveformA,rxWaveform); % Can't due to the values being complex
    [h, FitInfo] = lasso(XRealImag,yRealImag);
    hl = h(:,50);
    h_2 = makeComplex(hl);
    % =========================================================================

    estiY = txWaveformA * h_2;
    MSEs(k - min(range) + 1) = abs(mean((rxWaveform - estiY).^2));
%     disp(MSE)
end

plot(range,MSEs)