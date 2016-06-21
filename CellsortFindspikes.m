function [spmat, spt, spc] = CellsortFindspikes(ica_sig, thresh, dt, deconvtau, normalization)
% [spmat, spt, spc, zsig] = CellsortFindspikes(ica_sig, thresh, dt, deconvtau, normalization)
%
% CELLSORT
%  Deconvolve signal and find spikes using a threshold
%
% Inputs:
%   ica_sig - nIC x T matrix of ICA temporal signals
%   thresh - threshold for spike detection
%   dt - time step
%   deconvtau - time constant for temporal deconvolution (Butterworth
%   filter); if deconvtau=0 or [], no deconvolution is performed
%   normalization - type of normalization to apply to ica_sig; 0 - no
%   normalization; 1 - divide by standard deviation and subtract mean
%
% Outputs:
%   spmat - nIC x T sparse binary matrix, containing 1 at the time frame of each
%   spike
%   spt - list of all spike times
%   spc - list of the indices of cells for each spike
%
% Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, 2009
% Email: eran@post.harvard.edu, mschnitz@stanford.edu
%

if size(ica_sig,2)==1
    ica_sig = ica_sig';
end

if (nargin>=3)&&(deconvtau>0)
    dsig = diff(ica_sig,1,2);
    sig = ica_sig/deconvtau + [dsig(:,1),dsig]/dt;
else
    sig = ica_sig;
end

if (nargin<2)
    thresh=3;
    fprintf('Using threshold = 3 s.d. \n')
end
switch normalization
    case 0 % Absolute units
        zsig = sig';
    case 1 % Standard-deviation
        zsig = zscore(sig');
end
pp1=[zsig(1,:);zsig(1:end-1,:)];
pp2=[zsig(2:end,:);zsig(end,:)];
spmat = sparse((zsig>=thresh)&(zsig-pp1>=0)&(zsig-pp2>=0));

if nargout>1
    [spt,spc] = find(spmat);
    spt = spt*dt;
end
