function out = qtan(vector)

    rawangle = atan(vector(2)/vector(1));
    out = rawangle;

    if vector(1) < 0
        out = out + pi;
    elseif vector(2) < 0
        out = out + 2*pi;
    end

end