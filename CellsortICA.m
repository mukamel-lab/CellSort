function [ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, ...
    mixedfilters, CovEvals, PCuse, mu, nIC, ica_A_guess, termtol, maxrounds)
% [ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, mixedfilters, CovEvals, PCuse, mu, nIC, ica_A_guess, termtol, maxrounds)
%
%CELLSORT
% Perform ICA with a standard set of parameters, including skewness as the
% objective function
%
% Inputs:
%   mixedsig - N x T matrix of N temporal signal mixtures sampled at T
%   points.
%   mixedfilters - N x X x Y array of N spatial signal mixtures sampled at
%   X x Y spatial points.
%   CovEvals - eigenvalues of the covariance matrix
%   PCuse - vector of indices of the components to be included. If empty,
%   use all the components
%   mu - parameter (between 0 and 1) specifying weight of temporal
%   information in spatio-temporal ICA
%   nIC - number of ICs to derive
%   termtol - termination tolerance; fractional change in output at which
%   to end iteration of the fixed point algorithm.
%   maxrounds - maximum number of rounds of iterations
%
% Outputs:
%     ica_sig - nIC x T matrix of ICA temporal signals
%     ica_filters - nIC x X x Y array of ICA spatial filters
%     ica_A - nIC x N orthogonal unmixing matrix to convert the input to output signals
%     numiter - number of rounds of iteration before termination
%
% Routine is based on the fastICA package (Hugo Gävert, Jarmo Hurri, Jaakko Särelä, and Aapo
% Hyvärinen, http://www.cis.hut.fi/projects/ica/fastica)
%
% Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, 2009
% Email: eran@post.harvard.edu, mschnitz@stanford.edu
%

fprintf('-------------- CellsortICA %s -------------- \n', date)

if (nargin<4) || isempty(PCuse)
    PCuse = [1:size(mixedsig,1)];
end
if (nargin<6) || isempty(nIC)
    nIC = length(PCuse);
end
if (nargin<7) || isempty(ica_A_guess)
    ica_A_guess = randn(length(PCuse), nIC);
end
if (nargin<8) || isempty(termtol)
    termtol = 1e-6;
end
if (nargin<9) || isempty(maxrounds)
    maxrounds = 100;
end
if isempty(mu)||(mu>1)||(mu<0)
    error('Spatio-temporal parameter, mu, must be between 0 and 1.')
end

% Check that ica_A_guess is the right size
if size(ica_A_guess,1)~= length(PCuse) || size(ica_A_guess,2)~=nIC
    error('Initial guess for ica_A is the wrong size.')
end
if nIC>length(PCuse)
    error('Cannot estimate more ICs than the number of PCs.')
end
    
[pixw,pixh] = size(mixedfilters(:,:,1));
npix = pixw*pixh;

% Select PCs
if mu > 0 || ~isempty(mixedsig)
    mixedsig = mixedsig(PCuse,:);
end
if mu < 1 || ~isempty(mixedfilters)
    mixedfilters = reshape(mixedfilters(:,:,PCuse),npix,length(PCuse));
end
CovEvals = CovEvals(PCuse);

% Center the data by removing the mean of each PC
mixedmean = mean(mixedsig,2);
mixedsig = mixedsig - mixedmean * ones(1, size(mixedsig,2));

% Create concatenated data for spatio-temporal ICA
nx = size(mixedfilters,1);
nt = size(mixedsig,2);
if mu == 1
    % Pure temporal ICA
    sig_use = mixedsig;
elseif mu == 0
    % Pure spatial ICA
    sig_use = mixedfilters';
else
    % Spatial-temporal ICA
    sig_use = [(1-mu)*mixedfilters', mu*mixedsig];
    sig_use = sig_use / sqrt(1-2*mu+2*mu^2); % This normalization ensures that, if both mixedfilters and mixedsig have unit covariance, then so will sig_use
end

% Perform ICA
[ica_A, numiter] = fpica_standardica(sig_use, nIC, ica_A_guess, termtol, maxrounds);

% Sort ICs according to skewness of the temporal component
ica_W = ica_A';

ica_sig = ica_W * mixedsig;
ica_filters = reshape((mixedfilters*diag(CovEvals.^(-1/2))*ica_A)', nIC, nx);  % This is the matrix of the generators of the ICs
ica_filters = ica_filters / npix^2;

icskew = skewness(ica_sig');
[icskew, ICord] = sort(icskew, 'descend');
ica_A = ica_A(:,ICord);
ica_sig = ica_sig(ICord,:);
ica_filters = ica_filters(ICord,:);
ica_filters = reshape(ica_filters, nIC, pixw, pixh);

% Note that with these definitions of ica_filters and ica_sig, we can decompose
% the sphered and original movie data matrices as:
%     mov_sphere ~ mixedfilters * mixedsig = ica_filters * ica_sig = (mixedfilters*ica_A') * (ica_A*mixedsig),
%     mov ~ mixedfilters * pca_D * mixedsig.
% This gives:
%     ica_filters = mixedfilters * ica_A' = mov * mixedsig' * inv(diag(pca_D.^(1/2)) * ica_A'
%     ica_sig = ica_A * mixedsig = ica_A * inv(diag(pca_D.^(1/2))) * mixedfilters' * mov

    function [B, iternum] = fpica_standardica(X, nIC, ica_A_guess, termtol, maxrounds)
        
        numSamples = size(X,2);
        
        B = ica_A_guess;
        BOld = zeros(size(B));
        
        iternum = 0;
        minAbsCos = 0;
        
        errvec = zeros(maxrounds,1);
        while (iternum < maxrounds) && ((1 - minAbsCos)>termtol)
            iternum = iternum + 1;
            
            % Symmetric orthogonalization.
            B = (X * ((X' * B) .^ 2)) / numSamples;
            B = B * real(inv(B' * B)^(1/2));
            
            % Test for termination condition.
            minAbsCos = min(abs(diag(B' * BOld)));
            
            BOld = B;
            errvec(iternum) = (1 - minAbsCos);
        end
        
        if iternum<maxrounds
            fprintf('Convergence in %d rounds.\n', iternum)
        else
            fprintf('Failed to converge; terminating after %d rounds, current change in estimate %3.3g.\n', ...
                iternum, 1-minAbsCos)
        end
    end
end

