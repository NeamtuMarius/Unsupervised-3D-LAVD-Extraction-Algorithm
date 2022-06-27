function [line] = cluster2line(line)
%ridgesOutliersDetector find point not belonging to the minimum spanning
%tree longest branch
%   

stopControl = 3;
while stopControl > 2
    
    D = pdist2(line,line);
    g = graph(D);
    T = minspantree(g);
    
    T1 = T;
    
    idx = find(T1.degree==1);
    line1 = line; line1(idx,:) = [];
    
    
    
    stopControl = size(line,1) - size(line1,1);
    line = line1;
    % disp(num2str(stopControl))
    % pause(.1)
end

line = lineSorting(line);


end

