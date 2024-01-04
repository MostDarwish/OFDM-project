%% TRANSMITTER
% Different modulation schemes will all have the same average energy per bit
% 
% For the fair comparison between their BER later

N_Bits = 4e6;
BitStream_Data = randi([0 1],[1, N_Bits]);



Eb = 4;
Modulation_Order = [2 4 16, 2 4 16];
RepCode_En = 0;

SymbolStream = {};
Constellation = {};

for i = 1:6
    if (i > 3)
        RepCode_En(i) = 1;
        BitStream = repelem(BitStream_Data,3);
        Eb = Eb/3;
    else 
        RepCode_En(i) = 0;
        BitStream = BitStream_Data;   
    end


    [SymbolStream{i}, Constellation{i}] = Modulate(BitStream,Eb,Modulation_Order(i));

    if(RepCode_En(i) ~= 1)
        scatterplot(SymbolStream{i},1,0,'r*')
        grid on
        axis([-6 6 -6 6])
    end

end
%% CHANNEL
% We will define the channel where it has the following proberties
% 
% The Real Component of the resultant vector has a gaussian ditribution:
% 
% $$h_R ~\aleph \left(\mu_R ,\sigma_R \right)$$
% 
% The Imaginary Component of the resultant vector has a gaussian ditribution:
% 
% $$h_I ~\aleph \left(\mu_I ,\sigma_I \right)$$
% 
% Notes:
% 
% 1- Zero mean, Equal Variance channel is a Rayleigh channel where no LoS exists.
% 
% 2-If any mean has a Non-Zero value, but still equal variances are used Ricean 
% channel will be created.
% 
% 3- The equal variances and zero mean proberty is what gives the phase of the 
% channel the uniform distribution.

SNR_dB =linspace(0,15,50);
SNR = 10.^(SNR_dB/10);
SNR = diag(SNR);

N_Symbols = {};   h = {};  N = {}; Channel_OUT = {};
for i=1:6
N_Symbols{i} = length(SymbolStream{i});

U_R = 0;       Sigma_R = sqrt(0.5);
U_I = 0;       Sigma_I = sqrt(0.5);
h_R = randn([1 N_Symbols{i}])*Sigma_R + U_R;
h_I = randn([1 N_Symbols{i}])*Sigma_I + U_I;

h{i} = h_R + h_I*1j;
if (i == 1)
figure
H_mag = histogram(abs(h{1}),'Normalization','probability');
title('PDF of Channel Magnitude')
figure
H_arg = histogram(angle(h{1}),'Normalization','probability');
title('PDF of Channel Phase')
end

%Noise power will be determined by the SNR value
%SNR = Eb/No
No = Eb*SNR^-1;

N_c = randn([1 N_Symbols{i}]);
N_s = randn([1 N_Symbols{i}]);
N{i} = N_c + N_s*1j;
N{i} = (No/2)*repmat(N{i},length(SNR_dB),1);
%% Sending the signal through the channel
% Since the channel is represented by delta (One Resorvable Path)
% 
% Then the channel output will be simple multiplication with the channel value, 
% where each symbol will only be scaled by a gain and added phase.
% 
% Then the White Gaussian Noise is added. 

Channel_OUT{i} = SymbolStream{i}.*h{i};
Channel_OUT{i} = Channel_OUT{i} + N{i};
end
%% RECEIVER
%% 1) Equalizing the Channel 
% *Since the channel gain is less than one it's expected that the noise is going 
% to be amplfied during equalization increasing symbol errors.*

BER{6} = {};
for i = 1:6
BER{i} = zeros(length(SNR_dB),N_Bits);
for l = 1:length(SNR_dB)
ReceivedSymbolStream = Channel_OUT{i}(l,:)./h{i};

% scatterplot(ReceivedSymbolStream,1,0,'r*');
% grid on
% axis([-6 6 -6 6])
%% 2) Applying Maximum Apsteriori decoding Plus Hard Decision Decoding

ReceivedBitStream = Decode(ReceivedSymbolStream, Constellation{i}, RepCode_En(i));
%% Evaluating Bit Errors

BER{i}(l,:) = xor(ReceivedBitStream,BitStream_Data);
end

BER{i} = sum(BER{i},2)/N_Bits;
end

%%
%%PLotting The 6 Curves over laid
figure
hold on
for i = 1:6
semilogy(SNR_dB, BER{i})
end
hold off
title('BER vs SNR')
xlabel('SNR(db)')
ylabel('BER')
ylim([0,0.002])
legend({'BPSK NoCode','QPSK NoCode','16-QAM NoCode', ...
    'BPSK WithCode','QPSK WithCode','16-QAM WithCode'});
%% 
% We note here that the case when we use encoding is worse than that with coding.
% 
% This is because we're using the *same energy per information*  for fair comparison 
% reasons, therefore each bit in the 1/3 repition encoding case will carry the 
% third power in the encoding case.
% 
% Evaluating that scenario on BPSK case:
% 
% $$\begin{array}{l}{\textrm{BER}}_{\textrm{Code}} =3P_e^2 -P_e^3 \;\;\;\textrm{where}\;P_e 
% \;\textrm{is}\;\textrm{the}\;\textrm{SER}\;\textrm{which}\;\textrm{is}\;\textrm{equal}\;\textrm{to}\;\textrm{error}\;\textrm{per}\;\textrm{individual},\textrm{not}\;\textrm{information},\textrm{bit}\;\textrm{in}\;\textrm{BPSK}\;\textrm{case}\\P_e 
% =\frac{1}{2}\textrm{erfc}\left(\sqrt{\left(\frac{\left(\frac{\textrm{Eb}}{3}\right)}{N_o 
% }\right)}\right)=\frac{1}{2}\textrm{erfc}\left(\sqrt{\left(\frac{E_b }{{3N}_o 
% }\right)}\right)\\\textrm{Then}:{\textrm{BER}}_{\textrm{CODE}} \approx \frac{3}{2}{\textrm{erfc}}^2 
% \left(\sqrt{\left(\frac{E_b }{{3N}_o }\right)}\right)\to \textrm{The}\;\textrm{decrease}\;\textrm{in}\;\textrm{Energy}\;\textrm{and}\;\textrm{multiplication}\;\textrm{by}\;3\;\textrm{will}\;\textrm{try}\;\textrm{to}\;\textrm{increases}\;\textrm{the}\;\textrm{BER}\;\textrm{while}\\\textrm{the}\;\textrm{squaring}\;\textrm{will}\;\textrm{try}\;\textrm{to}\;\textrm{decrease}\;\textrm{it}\;\textrm{since}\;\textrm{the}\;\textrm{outcome}\;\textrm{of}\;\textrm{erfc}\;\textrm{function}\;\textrm{is}\;\textrm{less}\;\textrm{than}\;\textrm{one}\\\textrm{The}\;\textrm{result}\;\textrm{is}\;\textrm{overall}\;\textrm{increase}\;\textrm{in}\;\textrm{the}\;\textrm{BER}\;\textrm{than}\;\textrm{the}\;\textrm{nocode}\;\textrm{case}\end{array}$$

function [SymbolStream, S] = Modulate(BitStream, Eb, Order)
%This function takes the bit stream and modulate it according
%to modulation order
%Modulation order varies from 2:BPSK 4:QPSK 16:16-QAM
%Eb: Energy per Bit

N_Bits = length(BitStream);

switch (Order)
    case 2
        A = sqrt(Eb);
        SymbolStream = (2*BitStream-1)*A;
        S = 2;
    
    case 4
        A = sqrt(Eb);
        %Creating the constellation:
        %Index represents the mapped bits - 1
        %ex: 00-> index:1 -> (-1-j)*A
        S(1) = - 1 - 1j;
        S(2) = - 1 + 1j;
        S(3) =   1 - 1j;
        S(4) =   1 + 1j;
        S = S*A;
        
        SymbolStream = reshape(BitStream,[2, N_Bits/2])';
        SymbolStream = num2str(SymbolStream);
        SymbolStream = bin2dec(SymbolStream)';
        SymbolStream = S(SymbolStream+1);
        
    case 16
        % Eb = 2.5E -> A = sqrt(E)
        A = sqrt(0.4*Eb);
        %Creating the constellation:
        %Index represents the mapped bits - 1
        %ex: 0000-> index:1 -> (-3-3j)*A
        S(1) = - 3 - 3j;
        S(2) = - 3 - 1j;
        S(3) = - 3 + 3j;
        S(4) = - 3 + 1j;
        S(5) = - 1 - 3j;
        S(6) = - 1 - 1j;
        S(7) = - 1 + 3j;
        S(8) = - 1 + 1j;
        S(9) =   3 - 3j;
       S(10) =   3 - 1j;
       S(11) =   3 + 3j;
       S(12) =   3 + 1j;
       S(13) =   1 - 3j;
       S(14) =   1 - 1j;
       S(15) =   1 + 3j;
       S(16) =   1 + 1j;
        S = S*A;
        
        SymbolStream = reshape(BitStream,[4, N_Bits/4])';
        SymbolStream = num2str(SymbolStream);
        SymbolStream = bin2dec(SymbolStream)';
        SymbolStream = S(SymbolStream+1);
        
        
        
    
end

end

function [Out_BitStream] = Decode(SymbolStream, S, RepCode_En)
%This function takes the Symbol Stream after channel inversion
%And Decodes it according maximum Liklihood decoding
%aka: minimum distance to contellation point,
%since we are using eqiprobable signals
%Decoded symbols are directly converted to bits
%S is an array that carries the constellation points
%The Function then applies  Hard Decision Decoding
%Where its index is the equivalent binary value + 1

N_Symbols = length(SymbolStream);
Order = length(S);

if(S == 2)
    %BPSK case
    SymbolStream = real(SymbolStream);
    Out_BitStream(SymbolStream>0) = 1;
    Out_BitStream(SymbolStream<0) = 0;
else
    SS_mod = repmat(SymbolStream,Order,1);
    Rep_con = transpose(repmat(S,N_Symbols,1));
    Distance = SS_mod - Rep_con;
    Distance = abs(Distance);
    
    [~, Index] = min(Distance);
    %Output Bit stream in decimal
    Out_BitStream = Index-1;
    Out_BitStream = decimalToBinaryVector(Out_BitStream,log2(Order),'MSBFirst')';
    Out_BitStream = reshape(Out_BitStream,1,N_Symbols*log2(Order));
end
%decoding the repetition-3 code if used 
if(RepCode_En == 1)
        Out_BitStream = reshape(Out_BitStream,3,[]);
        Out_BitStream = mode(Out_BitStream);
end 
end