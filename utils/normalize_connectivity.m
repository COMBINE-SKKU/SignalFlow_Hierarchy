function A = normalize_connectivity(Connmat,flag)

[~,D] = eig(Connmat);
Ddiag= diag(D);
maxval = max(real(Ddiag));
A = Connmat/(maxval+1);

if flag
    figure; imagesc(A);
end
