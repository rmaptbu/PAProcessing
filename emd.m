function imf = emd(x)
% Empiricial Mode Decomposition (Hilbert-Huang Transform)
% imf = emd(x)
% Func : findpeaks

x0=x(:)';
N = length(x);
n = floor(N/2.0);
x = extend(x0);
imf = [];
while ~ismonotonic(x)
    x1 = x;
    sd = Inf;
    while (sd > 0.25) ||  ~isimf(x1)
        s1 = getspline(x1);
        s2 = -getspline(-x1);
        x2 = x1-(s1+s2)/2;
        
        sd = sum((x1-x2).^2)/sum(x1.^2);
        x1 = x2;
    end
    imf{end+1} = x1(n+1:end-n);
    x          = x-x1;
end
imf{end+1} = x(n+1:end-n);

%% FUNCTIONS
function u = ismonotonic(x)
u1 = length(findpeaks(x))*length(findpeaks(-x));
if u1 > 0
    u = 0;
else %there are no peaks or troughs
    u = 1;
end

function u = isimf(x) %condition 1 of IMF
N  = length(x);
u1 = sum(x(1:N-1).*x(2:N) < 0); %number of zero crossings
u2 = length(findpeaks(x))+length(findpeaks(-x)); %number of extrema
if abs(u1-u2) > 1 %zero crossings and extrema are NOT the same
    u = 0;
else
    u = 1;
end

function s = getspline(x)
N = length(x);
[~, p] = findpeaks(x);
%p=[1, p, N]; %clamp spline fitting
s = spline(p,x(p),1:N);

function xe = extend(x)
N = length(x);
n = floor(N/2.0);
e1 = x(n+1:-1:2).*[0:n-1]/(n-1);
e2 = x(end-1:-1:end-n).*[n-1:-1:0]/(n-1);
xe = [e1, x, e2];