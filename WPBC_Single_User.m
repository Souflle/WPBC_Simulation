%% 
% *Summary*
%% 
% * *Fading Coefficient h = complex constatnt*
% * *MMSE Estimation using Orthogonal (for time) Pilot Sequence*
%% *System Parameter*

Fc = 915e6;
Ae = (3e8)^2 / (4*pi*Fc^2);  % Effective aperture of isotropic tag antenna

T = 200;                    % 200 symbol for coherence time
M = 8;   % # of Transmit Antenna
R = 4;   % # of Receive Antenna
L = 4;    % pilot length per antenna
K = 1 ;                      % Single User Case

del = 0.25;
eff = 0.65;    % rectifier efficiency
fb_bit = 8;     % # of feedback bits    
distance = 4;  % distance from reader

noise_dBm = -20;  % dB Noise
sigma = 1e-3 * 10 ^ (noise_dBm/10);            % Noise Power
%% Orthogonal Pilot Sequence

D = L * M;                    % Pilot slot length
T_WET = T - K * D;           % WET slot length
  
X_pilot= eye(M);   
 
X_T = repmat(sqrt(2) * X_pilot, 1, L);     % Transmit Power : 2 W  
X = X_T.';
%% Tag Model
% $$\delta \;:\mathrm{reflection}\;\mathrm{coefficient}:\left\lbrack 0,1\right\rbrack 
% \;1\to \mathrm{reflect},0\to \mathrm{absorb}$$
% 
% $$\eta :\mathrm{rectifier}\;\mathrm{efficiency}$$

total_state = 2 ^ fb_bit;              % # of state
%% *Channel Model*

PL = Ae / (4 * pi * distance^2) ;     % Path loss from Distance
sigma_h = sqrt(PL) * 1;               % Multi-path fading variance

% Fading Coefficient of F & B channel
h_f = sigma_h / sqrt(2) * (randn(M, 1) + 1j * randn(M, 1));
h_b = sigma_h / sqrt(2) * (randn(R, 1) + 1j * randn(R, 1));

% Noise of F & B channel for each slot
n = sigma / sqrt(2) * (randn(1, D) + 1j * randn(1, D));
w = sigma / sqrt(2) * (randn(R, D) + 1j * randn(R, D));

n_WET = sigma / sqrt(2) * (randn(1, T_WET) + 1j * randn(1, T_WET));
%% 
%% Forward Channel
% Received signal in Tag: $b={\left(h^{\;f} \right)}^{\;T} X+n\in C^{\;1\textrm{XD}}$
% 
% Received Power in Tag: $P=\eta \left(1-\delta \;\right)E\left\lbrace b^{\;2} 
% \right\rbrace$
% 
% Incident Power (compensated): $P_{\textrm{inc}} =\left|b{{\left|\right.}^2 
% =}^{\;} \frac{P}{\eta \;\left(1-\delta \;\right)\;}\right.$

b = h_f.' * X_T + n;  
P_inc = abs(b).^ 2;     % Compensated incident power

% Quantizing Incident Power in tag 
bin = 6 * PL/ (total_state - 1);
partition = [bin: bin: 6 * PL];
codebook = [partition, 9 * PL];
[idx, rn] = quantiz(P_inc, partition, codebook);     % rn = quantized incident power
%% Backward Channel
% Received Signal in HAP: $Y=\sqrt{\;\delta \;}h^{\;b} \;\left({\left(h^f \right)}^T 
% X^T +n\right)+w=\sqrt{\;\delta \;}{\textrm{HX}}^T +N$
% 
% Total Noise: $N=\sqrt{\;\delta \;}h^b n+w\in C^{\textrm{RXD}} \;$
% 
% Total Noise Variance: $\sigma_N^2 =E\left\lbrace N^H N\right\rbrace =R\sigma^2 
% \;\left(1+\delta \;\sigma_{\;h}^2 \right)I_{\;D}$

Y = sqrt(del) * h_b * b + w;
sigma_N = sqrt(R * sigma ^ 2 * (1 + del * sigma_h ^ 2));
%% MMSE Estimation
% $$\sigma {\;}_{\mathrm{BS}}^2 =R\delta \;\sigma {\;}_h^4 \;$$
% 
% $$\hat{H} =\delta {\;}^{-\frac{1}{2}} \;Y_{\mathrm{CE}} {X\left(R_{\mathrm{hh}} 
% X^T X+\sigma_N^2 I_M \right)}^{\;-1} R_{\mathrm{hh}}$$
% 
% $$\mathrm{Error}\;\mathrm{Covariance}:R_E ={\left(R_{\mathrm{hh}}^{-1} +\delta 
% \;R_{\mathrm{NN}}^{-1} X^T X\right)}^{-1}$$ 

sigma_BS = sqrt(R * sigma_h ^ 4 * del);

H_mmse = (del)^(-1/2) * Y * X * inv(sigma_BS * X_T * X + sigma_N^2 * eye(M)) * sigma_BS % estimated BS-CSI
H = h_b * h_f.'     % real BS-CSI

MSE = 1 / (1 / (R * del * sigma_h ^ 4) + L * 1 / (sigma_N ^ 2 / del ));
error = H_mmse - H;
error_cov = error' * error;
%% 
% 
%% F & B-CSI Magnitude Estimation

mag_f_est = zeros(M, 1);
mag_b_est = zeros(R, 1);

for i = 1:M
    ith_path = repmat(X_pilot(:, i), L, 1);
    mag_f_est(i) = sqrt(rn * ith_path / L / 2);  % Because 2 W Transmit Power 
end

for i = 1: M
    mag_b_est = mag_b_est + abs(H_mmse(:, i)) / mag_f_est(i) / M;
end
%% F & B-CSI Phase Estimation

phase_H = angle(H_mmse);
phase_f_est = phase_H(1, :);

phase_b_tilda = angle(H_mmse ./ exp(1j * ones(R, 1) * phase_f_est));
phase_b_est = phase_b_tilda(:, 1);
%% Energy Harvesting
% Estimated path gain: $g^f ={\left(\hat{\left|h^f \right.} \left|\right.\right)}^2$
% 
% TX Beamforming: $\Phi =\frac{g^f }{\left|g^f \right|}e^{{j\left({\hat{\Phi 
% \;} }^f \right)}^* }$    -> Supply high power to good path & normalize to 1W
% 
% WET harvest power = $\eta \;\left(1-\delta \;\right)\cdot \left({\left|\left(h^f 
% \right)\right.}^T \Phi {\;+\;n\left|\right.}^2 \right)\cdot \;\left(2W\right)$                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

TX_beamforming = mag_f_est / norm(mag_f_est) .* exp(1j * phase_f_est.');
bs_1 = conj(H_mmse(1, :)) ./ norm(H_mmse(1, :));


MAX_harvest = 2 * eff * (1 - del) * (h_f' * h_f)  % Maximum Harvestable power corresponds to fading channel
BS_harvest = 2 * eff * (1 - del) * abs(bs_1 * h_f) .^ 2     % Just using BS channel
WET_harvest = 2 * eff * (1 - del) * abs(TX_beamforming' * h_f) .^ 2     % Proposition


WET_vs_BS = WET_harvest / BS_harvest    % larger than 1, Proposition is more efficient
WET_vs_MAX = WET_harvest / MAX_harvest

Pilot_harvest = 2 * eff * (1-del) * mean(P_inc);
harvest_rate = ((T-T_WET) * Pilot_harvest + T_WET * WET_harvest) / T;
%% Verification

mag_f_real = abs(h_f);
mag_b_real = abs(h_b);

mag_f_acc = mag_f_est ./ mag_f_real    % Accurately estimated when accuracy closed to 1
mag_b_acc = mag_b_est ./ mag_b_real ;   

phase_f_real = angle(h_f);
phase_b_real = angle(h_b);

phase_f_acc = angle(exp(1j * phase_f_est.') ./ exp(1j * phase_f_real))  % Well alligned when each components are similar
phase_b_acc = angle(exp(1j * phase_b_est) ./ exp(1j * phase_b_real));