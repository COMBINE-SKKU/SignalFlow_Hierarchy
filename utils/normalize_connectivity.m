function A = normalize_connectivity(Connmat,flag)
if nargin < 2 || isempty(flag)
    flag = 0;
end

eigvals = eig(Connmat);
radius = max(abs(eigvals));
A = Connmat/(radius+1);

if flag
    figure; imagesc(A);
end
