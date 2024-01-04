OFDM_Symbol_Stream = [1 1 1 2 2 2 5]'
Ch = [1 1 1 2 4;1 2 4 1 1;1 1 1 1 2;2 2 5 5 5;1 1 1 1 2;1 1 1 1 2;1 2 4 1 2]
CONV(OFDM_Symbol_Stream,Ch)

function [Channel_Out] = CONV(OFDM_Symbol_Stream,Ch)
    Channel_Out = OFDM_Symbol_Stream .* Ch
    Channel_Out = [Channel_Out;zeros(size(Ch,2)-1,size(Ch,2))]
    
    for i = 2:size(Ch,2)
       Channel_Out(:, i) = circshift(Channel_Out(:, i), i-1)
    end
    Channel_Out = sum(Channel_Out, 2)
end