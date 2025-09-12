%% Estimate tail indices
%
% First version: Richard Tol, 8 November 2011
% This version: Richard Tol, 9 March 2020

display('Estimate tail indices');

N = length(SCCs);
K = round(N/5);
weight = 'quality';

switch weight
    case 'author'
        for k=10:K+1,
            [Hill(k) Hillsd(k) DJV1(k) DJV2(k) AM(k) AMsd(k) T1(k) T1sd(k) T2(k) T3(k) D(k) Dsd(k)] = TailWHill(SCCs,AuthorWeight,k);
            [KR(k) KRsd(k) SS1(k) SS2(k) SS2sd(k) BF(k) BFsd(k) AM2(k) AM2sd(k) GI(k) GIsd(k)] = TailWZipf(SCCs,AuthorWeight,k);
        end
    case 'quality'
        for k=10:K+1,
            [Hill(k) Hillsd(k) DJV1(k) DJV2(k) AM(k) AMsd(k) T1(k) T1sd(k) T2(k) T3(k) D(k) Dsd(k)] = TailWHill(SCCs,QualWeight,k);
            [KR(k) KRsd(k) SS1(k) SS2(k) SS2sd(k) BF(k) BFsd(k) AM2(k) AM2sd(k) GI(k) GIsd(k)] = TailWZipf(SCCs,QualWeight,k);
        end
    case 'censor'
        SCCsc = SCCs(Censored>0);
        Censc = Censored(Censored>0);
        for k=10:K+1,
            [Hill(k) Hillsd(k) DJV1(k) DJV2(k) AM(k) AMsd(k) T1(k) T1sd(k) T2(k) T3(k) D(k) Dsd(k)] = TailWHill(SCCsc,Censc,k);
            [KR(k) KRsd(k) SS1(k) SS2(k) SS2sd(k) BF(k) BFsd(k) AM2(k) AM2sd(k) GI(k) GIsd(k)] = TailWZipf(SCCsc,Censc,k);
        end
    otherwise
        for k=10:K+1,
            [Hill(k) Hillsd(k) DJV1(k) DJV2(k) AM(k) AMsd(k) T1(k) T1sd(k) T2(k) T3(k) D(k) Dsd(k)] = TailHill(SCCs,k);
%ignore DJV because no standard deviation
%ignore T2, T3 because variants of T1
            [KR(k) KRsd(k) SS1(k) SS2(k) SS2sd(k) BF(k) BFsd(k) AM2(k) AM2sd(k) GI(k) GIsd(k)] = TailZipf(SCCs,k);
%ignore KR, SS because variants of BF
        end
end

%% Huisman
[HH HHsd] = TailHuisman(Hill(10:K)',10,K);
[HA HAsd] = TailHuisman(AM(10:K)',10,K);
[HT HTsd] = TailHuisman(T1(10:K)',10,K);
[HD HDsd] = TailHuisman(D(10:K)',10,K);
[HB HBsd] = TailHuisman(BF(10:K)',10,K);

%% plot
m = 100;

ind=m;
for i=2:K-m+1,
    ind(i)=1+ind(i-1);
end
marty = ones(K-m+1,1);
plot(ind,marty,'-k','HandleVisibility','off');
hold on
plot(ind,Hill(m:K),'-c')
plot(ind,Hill(m:K)+2*Hillsd(m:K),'--c','HandleVisibility','off')
plot(ind,Hill(m:K)-2*Hillsd(m:K),'--c','HandleVisibility','off')
plot(ind,AM(m:K),'-r')
plot(ind,AM(m:K)+2*AMsd(m:K),'--r','HandleVisibility','off')
plot(ind,AM(m:K)-2*AMsd(m:K),'--r','HandleVisibility','off')
plot(ind,T1(m:K),'-g')
plot(ind,T1(m:K)+2*T1sd(m:K),'--g','HandleVisibility','off')
plot(ind,T1(m:K)-2*T1sd(m:K),'--g','HandleVisibility','off')
plot(ind,D(m:K),'-m')
plot(ind,D(m:K)+2*Dsd(m:K),'--m','HandleVisibility','off')
plot(ind,D(m:K)-2*Dsd(m:K),'--m','HandleVisibility','off')
plot(ind,BF(m:K),'-b')
plot(ind,BF(m:K)+2*BFsd(m:K),'--b','HandleVisibility','off')
plot(ind,BF(m:K)-2*BFsd(m:K),'--b','HandleVisibility','off')
hold off
xlabel('number of highest observations used')
ylabel('tail index')
legend('Maximum likelihood','Best linear unbiased','Least squares','Moment','Quantile-quantile')