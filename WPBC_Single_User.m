%% 
% *Summary*
%% 
% * *Fading Coefficient h = complex constatnt*
% * *MMSE Estimation using Orthogonal (for time) Pilot Sequence*
%% *System Parameter*

Fc = 915e6;
Ae = (3e8)^2 / (4*pi*Fc^2);  % Effective aperture of isotropic tag antenna
K = 1;                      % Single User Case
T = 200;                    % 200 symbol for coherence time
M = 8;   % # of Transmit Antenna
R = 4;   % # of Receive Antenna
%% Orthogonal Pilot Sequence

L = 4;    % pilot length per antenna
D = L * M;                    % Pilot slot length
T_WET = T - K * D;           % WET slot length
  

X_pilot= eye(M);   
 
X_T = repmat(sqrt(2) * X_pilot, 1, L);     % Transmit Power : 2 W  
X = X_T.';
%% Tag Model
% $$\delta \;:\textrm{reflection}\;\textrm{coefficient}:\left\lbrack 0,1\right\rbrack 
% \;1\to \textrm{reflect},0\to \textrm{absorb}$$
% 
% $$\eta :\textrm{rectifier}\;\textrm{efficiency}$$

del = 0.25;
eff = 0.65;    % rectifier efficiency
fb_bit = 8;     % # of feedback bits      
total_state = 2 ^ fb_bit;              % # of state
distance = 4;  % distance from reader
%% *Channel Model*

PL = Ae / (4 * pi * distance^2);      % Path loss from Distance
sigma_h = sqrt(PL) * 1;               % Multi-path fading variance
SNRdB = 100;  % dB
sigma = 10 ^ -(SNRdB/10);            % Noise 

% Fading Coefficient of F & B channel
h_f = sigma_h / sqrt(2) * (randn(M, 1) + 1j * randn(M, 1))
h_b = sigma_h / sqrt(2) * (randn(R, 1) + 1j * randn(R, 1))

% Noise of F & B channel for each slot
n = sigma / sqrt(2) * (randn(1, D) + 1j * randn(1, D));
w = sigma / sqrt(2) * (randn(R, D) + 1j * randn(R, D));

n_WET = sigma / sqrt(2) * (randn(1, T_WET) + 1j * randn(1, T_WET));
w_WET = sigma / sqrt(2) * (randn(R, T_WET) + 1j * randn(R, T_WET));
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
bin = 1e-3 / (total_state - 1);
partition = [bin: bin: 1e-3];
codebook = [partition, 1.5e-3];
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
% $$\sigma {\;}_{\textrm{BS}}^2 =R\delta \;\sigma {\;}_h^4 \;$$
% 
% $$\hat{H} =\delta {\;}^{-\frac{1}{2}} \;Y_{\textrm{CE}} {X\left(R_{\textrm{hh}} 
% X^T X+\sigma_N^2 I_M \right)}^{\;-1} R_{\textrm{hh}}$$
% 
% $$\textrm{Error}\;\textrm{Covariance}:R_E ={\left(R_{\textrm{hh}}^{-1} +\delta 
% \;R_{\textrm{NN}}^{-1} X^T X\right)}^{-1}$$ 

sigma_BS = sqrt(R * sigma_h ^ 4 * del);
H_mmse = (del)^(-1/2) * Y * X * inv(sigma_BS * X_T * X + sigma_N^2 * eye(M)) * sigma_BS  % estimated BS-CSI
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
% TX Beamforming: $\Phi =\frac{g^f }{\left|g^f \right|}$

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
Maximum_harvesting_power = 2 * PL * M * eff * (1-del)  % fading channel power가 1일 때 Maximum power
%% Verification

mag_f_real = abs(h_f);
mag_b_real = abs(h_b);

mag_f_acc = mag_f_est ./ mag_f_real;    % 1에 가까울 수록 정확하게 추정한 것
mag_b_acc = mag_b_est ./ mag_b_real;    % 맞는지 잘 모르겠음 ... 맞는것같은데..... 왜 각 요소가 같징

phase_f_real = angle(h_f);
phase_b_real = angle(h_b);

phase_f_acc = angle(exp(1j * phase_f_est.') ./ exp(1j * phase_f_real));  % 각 요소가 같으면 잘 allign된 것
phase_b_acc = angle(exp(1j * phase_b_est) ./ exp(1j * phase_b_real));