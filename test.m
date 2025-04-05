
Fs = 48*10^3; % Sampling frequency
t = 0:1/Fs:1; % Time vector


function output= analyse_t(t, x)
    plot(t, abs(x));
    xlim([0, 0.003]);
    grid on;
end
    
    function output= analyse_f(Fs,X)
        N = length(X);
        f = (0:N/2-1)*(Fs/N);
        plot(f, abs(X(1:floor(N/2))));
        % xlim([0, 3*10^4]);
        grid on;
    end

% Defining 10 sinusoidal signals of different frequencies
x = sin(2*pi*t) + sin(2*pi*10*t) +sin(2*pi*10^2*t) +sin(2*pi*500*t) +sin(2*pi*10^3*t) +sin(2*pi*2500*t) +sin(2*pi*7400*t) +sin(2*pi*5900*t) +sin(2*pi*10^4*t) + sin(2*pi*1.5*10^4*t);

%Computing FFT
X = fft(x);

%Loading the filter
load('filter.mat');

%Filtering the generated signals using filter coefficients
y = filter(Num, 1, x);

% analyse_t(t, x);
% hold on;
% analyse_t(t, y);
% hold off;
 
%Computing the FFT of the filtered signal.
Y = fft(y);

analyse_f(Fs, X);
hold on;
analyse_f(Fs, Y);
hold off;