function [VAR,F] = StateSpaceGC(ts,moselect,lambda)


regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'OLS';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 10;     % maximum model order for model order estimation

if moselect
    [~,~,moAIC,moBIC] = tsdata_to_infocrit(ts,momax,icregmode);
    % Select model order.

    if     strcmpi(morder,'actual')
        morder = amo;
        fprintf('\nusing actual model order = %d\n',morder);
    elseif strcmpi(morder,'AIC')
        morder = moAIC;
        fprintf('\nusing AIC best model order = %d\n',morder);
    elseif strcmpi(morder,'BIC')
        morder = moBIC;
        fprintf('\nusing BIC best model order = %d\n',morder);
    else
        fprintf('\nusing specified model order = %d\n',morder);
    end
else
    morder = 1;
end

[A,SIG] = tsdata_to_var(ts,morder,regmode,lambda); 
VAR = A(:,:,1); 
VAR = VAR - diag(diag(VAR));
% assert(~isbad(A),'VAR estimation failed - bailing out');
% 
% info = var_info(A,SIG);
% assert(~info.error,'VAR error(s) found - bailing out');

[F] = var_to_pwcgc(A,SIG);
% assert(~isbad(F,false),'GC calculation failed - bailing out');

F(isnan(F)) = 0;
% 
% nvars = size(ts,1);
% nobs = size(ts,2);
% ntrials = size(ts,3);
% tstat = '';
% alpha = 0.05;
% 
% pval_t = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat);
% sig_t  = significance(pval_t,alpha,'NONE');
% 
% VAR(sig_t==0)=0;
% F(sig_t==0)=0;
