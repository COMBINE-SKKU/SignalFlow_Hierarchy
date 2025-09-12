bold_sig = BOLDs{1,1}(1:1200,1:360);

nobs=size(bold_sig,1);  
TR = .72;
bands=[0.01 0.1]; %bandpass filter lower and upper bound
data = rsHRF_band_filter(bold_sig,TR,bands);
sigma = std(data);

para.TR = TR;
BF = {'Canonical HRF (with time derivative)'
'Canonical HRF (with time and dispersion derivatives)'              
'Gamma functions'
'Fourier set'
'Fourier set (Hanning)'};
% choose the set of basis functions THIS MUST BE AN INPUT
bf_id = 2;
para.name = BF{bf_id}; % Gamma functions
para.order = 3; % for Gamma functions or Fourier set

temporal_mask = []; % without mask, it means temporal_mask = logical(ones(nobs,1)); i.e. all time points included. nobs: number of observation = size(data,1). if want to exclude the first 1~5 time points, let temporal_mask(1:5)=0;
% temporal_mask = logical(ones(nobs,1));  temporal_mask(5:15)=0;

para.T  = 5; % magnification factor of temporal grid with respect to TR. i.e. para.T=1 for no upsampling, para.T=3 for 3x finer grid
para.T0 = 1; % position of the reference slice in bins, on the grid defined by para.T. For example, if the reference slice is the middle one, then para.T0=fix(para.T/2)
if para.T==1
    para.T0 = 1;
end

para.dt  = para.TR/para.T; % fine scale time resolution.
para.AR_lag = 1; % AR(1) noise autocorrelation.
para.thr = 1; % (mean+) para.thr*standard deviation threshold to detect event.

para.len = 24; % length of HRF, in seconds

min_onset_search = 3; % minimum delay allowed between event and HRF onset (seconds)
max_onset_search = 6; % maximum delay allowed between event and HRF onset (seconds)
para.lag  = fix(min_onset_search/para.dt):fix(max_onset_search/para.dt);


data_deconv = zeros(nobs,40);
event_bold = cell(360,1);
T = round(para.len/TR);

for roi = 1:40
    [beta_hrf, bf, event_bold{roi,1}] = rsHRF_estimation_temporal_basis(data(:,roi),para,temporal_mask);
    hrf= bf*beta_hrf(1:size(bf,2),:); %HRF
    PARA = rsHRF_get_HRF_parameters(hrf,para.dt);
    hrf = resample(hrf,1,para.T);
    zdata = zscore(bold_sig(:,roi));
    sigma_ind = sigma(roi);
    data_deconv(:,roi) = zscore(rsHRF_iterative_wiener_deconv(zdata,hrf,100));  
end

event_number=length(event_bold{2,1}{1});
event_plot=nan(1,nobs);
event_plot(event_bold{2,1}{1})=1;
figure(1);plot((1:length(hrfa(:,1)))*TR,hrfa(:,1),'b');xlabel('Time (s)')
title(['HRF (',BF{bf_id},')'])
figure(2);plot((1:nobs)*TR,zscore(data(:,2)));
hold on;plot((1:nobs)*TR,zscore(data_deconv(:,2)),'r');
stem((1:nobs)*TR,event_plot,'k');legend('BOLD','Deconvolved BOLD','BOLD events');xlabel('Time (s)')

data_deconv = zscore(data_deconv,0,2);
rDCM2 = run_rdcm(bold_sig, 1.6); 
var2 = var_ols_ridge(data_deconv',100000000);
figure; imagesc(var2);