%% random data
nb = 256;
data_in = randi([0 1],1,nb);
%% Interleaver
m = 16; L = 16; % Number of rows(m) and columns(L)
nb = m*L;  % Size of interleaver
k = 0:nb-1;  % interleaver output indexing vector
i = L*mod(k,m)+floor(k/m);  % intermediate variable
j = (nb/2)*floor(2*i/nb)+mod((i+nb),nb/2); % interleaver input indexing vector
data_out = data_in(j+1);    % interleaver equation
data_rec(j+1) = data_out;   % de-interleaver equation
sum(abs(data_in-data_rec))  % number of errors