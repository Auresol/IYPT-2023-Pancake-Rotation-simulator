clear all;
n = 20;

autorunner(20);

omegaarray = zeros(n);

[result,omegaarray(2)] = autorunner(2);
sz = size(result);
resultarr = zeros(sz(2),n-1);
resultarr(:,1) = result;

for nin = 2:n
    [result,omegaarray(nin)] = autorunner(nin);
    resultarr(:,nin-1) = result(:);
    disp("PASS " + nin);
end

plot(resultarr);