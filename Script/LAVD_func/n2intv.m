function [intv] = n2intv(n,li)
%n2intv dato un numero e la lunghezza degli intervalli restituisce la
%matrice degli intervalli
%   esempio:
%   n = 25; li = 6;
%   intv = [1     7    13    19    25;
%           6    12    18    24    25]

ini = 1:li:n;
fin = [ini + li - 1]; fin(end) = n;

intv = [ini; fin];


end

