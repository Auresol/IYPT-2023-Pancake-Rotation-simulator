function out = calculate_angle(A,B)
    if A(1) > 0 && A(2) > 0 && B(1) > 0 && B(2) < 0 % A = Q1, B = Q4
        out = -2*pi + qtan(A) - qtan(B);
    elseif A(1) > 0 && A(2) < 0 && B(1) > 0 && B(2) > 0 % A = Q4, B = Q1
        out = 2*pi + qtan(A) - qtan(B);
    else
        out = qtan(A) - qtan(B);
    end
end