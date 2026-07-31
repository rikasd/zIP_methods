function [pxy,f_Hz,ixy,fft_x_mean,fft_y_mean] = ...
    cpsd_custom(x,y,window,n_overlap,nfft,mimo_option,fs,detrend_option)
% CPSD_CUSTOM Custom function to compute the cross power spectral density.
% (CPSD)
%
% Inputs: 
% x,y = signals: Nt x m, Nt x n
%       Nt: number of time samples
%       m:  number of x signals
%       n:  number of y signals
% window = data window [default: hamming window]
%          If window is an integer, assumed to indicate data window length
% n_overlap = number of overlapping time samples
% nfft  = number of FFT points used to calculate the CPSD estimate
% mimo_option = 'mimo'
%       - If m == n, default is not mimo, and output gives CPSD between 
%         each column vector of x and y: pxy = Nf x m
%       - If 'mimo' or m ~= n, output gives CPSD matrix, for all combos
%         of x and y: pxy = Nf x m x n                       
% fs    = sampling frequency of x,y
% detrend_option = 1: no detrend (default), 2: zeromean, 3: linear detrend
% 
% Outputs:
% pxy = CPSD between input signals x and y; nfft x 1
% f_Hz = frequency of CPSD estimates; nfft x 1
% ixy = CPSD estimates from each data window;
%       nfft x ( Nt/(length(window)-n_overlap) + 1 )
%
%
% Rika Sugimoto-Dimitrova (rikasd@mit.edu)
% 2023-10-26
% Last updated: 2024-02-16

% CPSD detrend settings
iRAW = 1;
iZEROMEAN = 2;
iDETREND = 3;

if nargin < 8
    detrend_option = iRAW; % default detrend option
end

Nt = length(x);
if length(y) ~= Nt
    error('Error using cpsd_custom\n x and y must have the same length.');
end

L = length(window); % window size
if  L == 1
    L = window;
    window = hamming(L);
end
window_norm_factor = window'*window;

if isempty(n_overlap)
    n_overlap = round(L/2);
end

if isempty(nfft)
    nfft = L*2;
end

if size(x,2) > size(x,1)
   x = x';
end
if size(y,2) > size(y,1)
   y = y';
end

bool_cpsd = false;
if strcmp(mimo_option,'mimo') && size(x,2) == size(y,2) 
    bool_cpsd = true;
elseif strcmp(mimo_option,'mimo') && size(x,2) ~= size(y,2) 
    error('Cannot compute mimo cpsd if inputs have different number of columns');
elseif size(x,2) ~= size(y,2) 
    bool_cpsd = true;
end

if size(window,2) ~= 1
   window = window';
   if size(window,2) ~= 1
       error('Make sure the window vector past is a column vector');
   end
end

K = max(1, floor(Nt / (L - n_overlap) - 1)); % number of segments
Nx = size(x,2);
Ny = size(y,2);

ixy = cell(Nx,Ny);

for i = 1:K

    x_seg = x((1:L)+(i-1)*(L - n_overlap),:);
    y_seg = y((1:L)+(i-1)*(L - n_overlap),:);
    switch detrend_option
        case iZEROMEAN
            x_seg = detrend(x_seg,0);
            y_seg = detrend(y_seg,0);
        case iDETREND
            x_seg = detrend(x_seg,1);
            y_seg = detrend(y_seg,1);
    end
    x_window = (x_seg.*window);
    y_window = (y_seg.*window);

    fft_x = fft(x_window,nfft);
    fft_y = fft(y_window,nfft);
    
    L_fft = round(length(fft_x)/2);
    f_Hz = fs*(0:(nfft/2))'/nfft;

    for m = 1:Nx
    for n = 1:Ny
        ixy{m,n}(:,i) = ...
            fft_x(1:L_fft+1,m).*conj(fft_y(1:L_fft+1,n)) ./window_norm_factor*2/fs;
    end
    end

    fft_x_all(:,:,i) = fft_x;
    fft_y_all(:,:,i) = fft_y;

end

if bool_cpsd
    for m = 1:Nx
    for n = 1:Ny
        pxy(:,m,n) = mean(ixy{m,n},2); % pxy(:,m,n) = cpsd(x(:,m),y(:,n))
    end
    end
else
    for m = 1:Nx
        pxy(:,m) = mean(ixy{m,m},2);
    end
end


for m = 1:Nx
    fft_x_mean(:,m) = mean(fft_x_all(:,m,:),3);
end

for n = 1:Ny
    fft_y_mean(:,n) = mean(fft_y_all(:,n,:),3);
end

end