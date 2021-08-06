%% signal processing
% baseline correction, noise reduction;
% algorithms: 
% 1) detrend
% 2) fits of polynominal
% 3) beads; 
% 4) wavelet transform; 
% 5) experience mode decomposition;
    
    clc
    clear all
    load sin-response;
    set(0,'defaultfigurecolor','w');
    x=x(600:900,:);

    % requirements for beads: array values should be positive integers;
    p1=ceil(x(:,1)); p2=ceil(x(:,2)); p3=ceil(x(:,3)); p4=ceil(x(:,4)); p5=ceil(x(:,5));
    % parameter configuration
    d = 1;          % d : filter order parameter (d = 1 or 2)
    r = 6;          % r : asymmetry parameter
    fc = 0.006;     % fc : cut-off frequency (cycles/sample)
    amp = 0.8; lam0 = 0.5*amp; lam1 = 5*amp; lam2 = 4*amp;
    [P1, f1, cost1] = beads(p1, d, fc, r, lam0, lam1, lam2);
    [P2, f2, cost2] = beads(p2, d, fc, r, lam0, lam1, lam2);
    [P3, f3, cost3] = beads(p3, d, fc, r, lam0, lam1, lam2);
    [P4, f4, cost4] = beads(p4, d, fc, r, lam0, lam1, lam2);
    [P5, f5, cost5] = beads(p5, d, fc, r, lam0, lam1, lam2);
    % display
    figure(1);  plot(x);  
    figure(2);  % plot for the estimated sparse-derivative signals
    subplot(5,1,1); plot(P1); subplot(5,1,2); plot(P2); subplot(5,1,3); plot(P3); subplot(5,1,4); plot(P4); subplot(5,1,5); plot(P5);    
    figure(3);  % plot for the estimated baseline signals
    subplot(5,1,1); plot(f1); subplot(5,1,2); plot(f2); subplot(5,1,3); plot(f3); subplot(5,1,4); plot(f4); subplot(5,1,5); plot(f5);     
    figure(4);  % plot signals that the estimated sparse-derivative signals subtract baselines
    subplot(5,1,1); plot(P1 - f1); subplot(5,1,2); plot(P2 - f2); subplot(5,1,3); plot(P3 - f3); subplot(5,1,4); plot(P4 - f4); subplot(5,1,5); plot(P5 - f5);    
    grid;
    drawnow;
    hold on;    

%% function of beads toolbox
function [x, f, cost] = beads(y, d, fc, r, lam0, lam1, lam2)
% [x, f, cost] = beads(y, d, fc, r, lam0, lam1, lam2)
% Baseline estimation and denoising using sparsity (BEADS)
% INPUT
%   y: Noisy observation
%   d: Filter order (d = 1 or 2)
%   fc: Filter cut-off frequency (cycles/sample) (0 < fc < 0.5)
%   r: Asymmetry ratio
%   lam0, lam1, lam2: Regularization parameters
% OUTPUT
%   x: Estimated sparse-derivative signal
%   f: Estimated baseline
%   cost: Cost function history
% Reference:
% Chromatogram baseline estimation and denoising using sparsity (BEADS)
% Xiaoran Ning, Ivan W. Selesnick, Laurent Duval
% Chemometrics and Intelligent Laboratory Systems (2014)
% doi: 10.1016/j.chemolab.2014.09.014
% The following parameter may be altered.
Nit = 100;       % Nit: Number of iterations
pen = 'L1_v2';  % pen : penalty function for sparse derivative ('L1_v1' or 'L1_v2')
EPS0 = 1e-6;    % cost smoothing parameter for x (small positive value)
EPS1 = 1e-6;    % cost smoothing parameter for derivatives (small positive value)
switch pen
    case 'L1_v1'
        phi = @(x) sqrt(abs(x).^2 + EPS1);
        wfun = @(x) 1./(sqrt(abs(x).^2 + EPS1));
    case 'L1_v2'
        phi = @(x) abs(x) - EPS1 * log(abs(x) + EPS1);
        wfun = @(x) 1./( abs(x) + EPS1);
    otherwise
        disp('penalty must be L1_v1, L1_v2')
        x = []; cost = []; f = [];
        return
end
theta = @(x) sum(x(x>EPS0)) - r * sum(x(x<-EPS0)) ...
    + sum( (1+r)/(4*EPS0)*x(abs(x)<=EPS0).^2 ...
    + (1-r)/2 * x(abs(x)<=EPS0) + EPS0*(1+r)/4 );
y = y(:);
x = y;
cost = zeros(1, Nit);
N = length(y);
[A, B] = BAfilt(d, fc, N);
H = @(x) B*(A\x);
e = ones(N-1, 1);
D1 = spdiags([-e e], [0 1], N-1, N);
D2 = spdiags([e -2*e e], 0:2, N-2, N);
D = [D1;  D2];
BTB = B'*B;
w = [lam1 * ones(N-1, 1); lam2 * ones(N-2, 1)];
b = (1-r)/2 * ones(N, 1);
d = BTB * (A\y) - lam0 * A' * b;
gamma = ones(N, 1);

for i = 1:Nit    
    Lambda = spdiags(w.*wfun(D*x), 0, 2*N-3, 2*N-3);    
    k = abs(x) > EPS0;
    gamma(~k) = ((1 + r)/4) / abs(EPS0);
    gamma(k) = ((1 + r)/4) ./  abs(x(k));
    Gamma = spdiags(gamma, 0, N, N);   
    M = 2 * lam0 * Gamma + D' * Lambda * D;
    x = A * ((BTB + A'*M*A)\d);    
    cost(i) = 0.5 * sum(abs(H(y - x)).^2) + lam0 * theta(x) ...
        + lam1 * sum(phi(diff(x))) + lam2 * sum(phi(diff(x, 2)));
end
f = y - x - H(y-x);
end

% --- local function ----
function [A, B] = BAfilt(d, fc, N)
% [A, B] = BAfilt(d, fc, N)
%
% Banded matrices for zero-phase high-pass filter.
% The matrices are 'sparse' data type in MATLAB.
% INPUT
%   d  : degree of filter is 2d (use d = 1 or 2)
%   fc : cut-off frequency (normalized frequency, 0 < fc < 0.5)
%   N  : length of signal
b1 = [1 -1];
for i = 1:d-1
    b1 = conv(b1, [-1 2 -1]);
end
b = conv(b1, [-1 1]);
omc = 2*pi*fc;
t = ((1-cos(omc))/(1+cos(omc)))^d;
a = 1;
for i = 1:d
    a = conv(a,[1 2 1]);
end
a = b + t*a;
A = spdiags( a(ones(N, 1), :), -d:d, N, N);   % A: Symmetric banded matrix
B = spdiags(b(ones(N, 1), :), -d:d, N, N);    % B: banded matrix
end

%% function of wavelet transform


%% function of experience mode decomposition
