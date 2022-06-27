function [sorted_line] = lineSorter(line)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dmin = 999999999999999;

% inizializzazione con i due punti piu vicini
p0 = line(1,:);
line(1,:) = [];
k1 = dsearchn(line,p0)
p1 = line(k1,:);
line(k1,:) = [];
sorted_line = [p0; p1];

n = size(line,1);
for i = 1:n
    i
    p0 = sorted_line(1,:);
    p1 = sorted_line(end,:);
    
    [k0,d0] = knnsearch(line,p0);
    [k1,d1] = knnsearch(line,p1);
    
    [~, ind] = min([d0 d1 dmin]);

    if ind == 1
        sorted_line = [line(k0,:); sorted_line];
        line(k0,:) = [];
    elseif ind == 2
        sorted_line = [sorted_line; line(k1,:)];
        line(k1,:) = [];
    else
        % line(k0,:) = [];
    end
    
    
    
    
end





end

