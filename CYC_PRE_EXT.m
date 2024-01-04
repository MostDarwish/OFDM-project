
function [Pre_Ext_SYM] = CYC_PRE_EXT(Symbol_Array, N_SubCarrier, Ext_Size)
%SymbolSize is represented by the Size of the symbol from IFFT
%aka: it's the #subcarriers per symbol
%Ext_Size: is the #samples to be cyclically extended

N_Symbols = length(Symbol_Array)/N_SubCarrier;

Pre_Ext_SYM = reshape(Symbol_Array,[N_SubCarrier, N_Symbols])';

Ext_Sym = Pre_Ext_SYM(:,[N_SubCarrier-Ext_Size+1,N_SubCarrier]);
if Ext_Size == 1
    Ext_Sym = Ext_Sym(:,1);
end

Pre_Ext_SYM = [Ext_Sym, Pre_Ext_SYM];

[X, Y] = size(Pre_Ext_SYM);
Pre_Ext_SYM = reshape(Pre_Ext_SYM', [1, X*Y]);
end

%%
function [EXT_REM] = REM_PRE_EXT(Symbol_Array, N_SubCarrier, Ext_Size)
%SymbolSize is represented by the Size of the symbol from IFFT
%aka: it's the #subcarriers per symbol
%Ext_Size: is the #samples to be cyclically extended

Ext_Sym_Size = N_SubCarrier + Ext_Size;
N_Symbols = length(Symbol_Array)/Ext_Sym_Size;

EXT_REM = reshape(Symbol_Array, [Ext_Sym_Size, N_Symbols])';
EXT_REM = EXT_REM(:,Ext_Size+1:end);

EXT_REM = reshape(EXT_REM', [1, N_Symbols*N_SubCarrier]);
end