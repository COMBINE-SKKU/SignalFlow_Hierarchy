function [Wwp,Wsp,Wssp]=geomsurr(W,D,nmean,nstd)
% [Wwp,Wsp,Wssp]=geomsurr(W,D,[nmean],[nstd])
% 
% Random graphs that preserve distance effect
%
% Inputs
% W         = weighted connectivity matrix
% D         = physical distance between nodes
% nmean     = order to preserve mean
% nstd      = order to preserve std
%
% Outputs
% Wwp       = random surrogate preserving weights but not strengths
% Wsp       = preserves strength distribution but not sequence
% Wssp      = preserves strength sequence
%
% Note: Wsp and Wssp are generated assuming that W is undirected.
%
% Reference: Roberts et al. (2016) NeuroImage 124:379-393.
%
% M Breakspear, J Roberts
% QIMR Berghofer, 2015

if nargin<3, nmean=3; end
if nargin<4, nstd =2; end

%Check if directed
drct=1;
if max(max(W-W'))==0, drct=0; end

%Preliminaries
N = size(W,1);
Wwp = zeros(N,N); 

W(1:N+1:end)=0;         %Ensure there's no self connections
if drct==0,             %Only do one triangle if undirected, then replicate
    W = triu(W);
end

nz= find(W);            %Nonzero entries
w = W(nz);              %Vector of non-zero connections
d = D(nz);              %Corresponding distances
logw = log(w);          %log-weights

%1. remove mean to nmean order
pfit1=polyfit(d,logw,nmean);
mnlogw=logw-polyval(pfit1,d,nmean);

%2. adjust variance to nstd order
pfit2=polyfit(d,abs(mnlogw),nstd);
stdlogw=mnlogw./polyval(pfit2,d,nstd);

%3. Now create surrogate data, adjusted for mean and std
%Shuffle the old ones
shuff=randperm(length(stdlogw));
surr=stdlogw(shuff);

%4. Now put the geometry back in
%4.1 Invert
stdsurr=surr.*polyval(pfit2,d,nstd);     % std
mnsurr=stdsurr+polyval(pfit1,d,nmean);   % and mean 
%4.2 Use surrogate weights as a scaffold to reorder the original weights,
%    thus preserving the original set of weights but in the new distance-
%    preserving random order (cf. the "amplitude adjustment" used in
%    Fourier surrogates for time series analysis)
surrlogw=rankreorder(logw,mnsurr);
%4.3 Undo logarithm and put into weight-preserving surrogate matrix
Wwp(nz)=exp(surrlogw);

%Make undirected if W is undirected
if drct==0;
    Wwp=(Wwp+Wwp');
    W=W+W';
end

%5. Adjust node strengths - NOTE: assumes W is undirected
strengthsW=sum(W);     % original node strengths
strengthsWwp=sum(Wwp); % new node strengths
%5.1 Re-order the old strengths to match new random order in Wwp
strengthsWsp=rankreorder(strengthsW,strengthsWwp);
%5.2 Adjust strengths to give both original and random strength sequences
Wsp=strengthcorrect(Wwp,strengthsWsp); % orig strengths in new sequence
Wssp=strengthcorrect(Wwp,strengthsW);  % orig strengths in orig sequence

end

function out=rankreorder(x,scaffold)
% reorder vector x to have same rank order as vector scaffold
y(:,1)=scaffold;
y(:,2)=1:length(x);
y=sortrows(y,1);
y(:,1)=sort(x);
y=sortrows(y,2);
out=y(:,1);
end

function sW=strengthcorrect(W,ss,nreps)
% sW=strengthcorrect(W,ss,[nreps])
%
% Iteratively rescale weights in connectivity matrix W to make the strength
% sequence converge onto desired strength sequence ss
%
% Inputs
% W         = weighted connectivity matrix
% ss        = desired strength sequence
% nreps     = number of iterations in the rescaling algorithm 
%             (optional, default=9)
%
% Outputs
% sW        = network with desired strength sequence ss
%
% Note: assumes W is undirected (symmetric)
%
% Reference: Roberts et al. (2016) NeuroImage 124:379-393.
%
% M Breakspear, J Roberts
% QIMR Berghofer, 2015

if nargin<3,
    nreps=9;
end

ss=ss(:)';
N=length(ss);
discnodes=sum(W)==0;
sW=W.*repmat(ss./sum(W),N ,1);
sW(:,discnodes)=0; %ensure disconnected nodes don't introduce NaNs
sW=(sW+sW')/2;
for j=1:nreps         
    sW=sW.*repmat(ss./sum(sW),N ,1);
    sW(:,discnodes)=0;
    sW=(sW+sW')/2;
end
end