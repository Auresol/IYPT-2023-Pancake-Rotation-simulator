function output = magnitude(vec)
%MAGNITUDE Summary of this function goes here
%   Detailed explanation goes here
    arguments
        vec (1,3)
    end
    
    output = sqrt(vec(1) * vec(1) + vec(2) * vec(2) + vec(3) * vec(3));
end

