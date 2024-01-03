clc;
clear;
%% Generating a random signal of length 4096
L = 4096;
random_signal = rand(1, L);
%% DFT execution time
tic;
X_dft = myDFT(random_signal);
time_dft = toc;
%% FFT execution time
tic;
X_fft = fft(random_signal);
time_fft = toc;
%% Implementation check
tolerance = 1e-4;
check_for_equality = isequal(size(X_dft), size(X_fft)) && all(abs(X_dft - X_fft) < tolerance);
disp('--------------- Implementation check ---------------');
fprintf('Are DFT and FFT results the same? \n %d \n', check_for_equality);
%% Displaying the execution times
disp('----------------- Execution times ------------------');
fprintf('DFT execution time: %.7f seconds \n', time_dft);
fprintf('FFT execution time: %.7f seconds \n', time_fft);
%% Comments
disp('-------------------- Comments ----------------------');
disp('There is a significant speed advantage and superior performance for the');
disp('built-in MATLAB function fft over the custom implementation of the DFT.');
disp('The reason for the great difference in execution times is that the FFT ');
disp('algorithm reduces the time complexity from O(N^2) to O(N log N), making');
disp('it much more efficient for larger input sizes.');
%% DFT function
function X = myDFT(x)
    N = length(x);
    X = zeros(1, N);
    for k = 0:N-1
        summation = 0;
        for n = 0:N-1
            summation = summation + x(n+1) * exp(-1j * 2 * pi * k * n/N);
        end
        X(k+1) = summation;
    end
end