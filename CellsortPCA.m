function [mixedsig, mixedfilters, CovEvals, covtrace, movm, ...
    movtm] = CellsortPCA(fn, flims, nPCs, dsamp, outputdir, badframes)
% [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(fn, flims, nPCs, dsamp, outputdir, badframes)
%
% CELLSORT
% Read TIFF movie data and perform singular-value decomposition (SVD)
% dimensional reduction.
%
% Inputs:
%   fn - movie file name. Must be in TIFF format.
%   flims - 2-element vector specifying the endpoints of the range of
%   frames to be analyzed. If empty, default is to analyze all movie
%   frames.
%   nPCs - number of principal components to be returned
%   dsamp - optional downsampling factor. If scalar, specifies temporal
%   downsampling factor. If two-element vector, entries specify temporal
%   and spatial downsampling, respectively.
%   outputdir - directory in which to store output .mat files
%   badframes - optional list of indices of movie frames to be excluded
%   from analysis
%
% Outputs:
%   mixedsig - N x T matrix of N temporal signal mixtures sampled at T
%   points.
%   mixedfilters - N x X x Y array of N spatial signal mixtures sampled at
%   X x Y spatial points.
%   CovEvals - largest eigenvalues of the covariance matrix
%   covtrace - trace of covariance matrix, corresponding to the sum of all
%   eigenvalues (not just the largest few)
%   movm - average of all movie time frames at each pixel
%   movtm - average of all movie pixels at each time frame, after
%   normalizing each pixel deltaF/F
%
% Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, 2009
% Email: eran@post.harvard.edu, mschnitz@stanford.edu
%

%
% June-August 2016 (Cellsort 1.41): Fixed an inconsistency between create_tcov and
% creat_xcov. In the current version, both methods agree with the previous
% version of create_tcov.
%

tic
fprintf('-------------- CellsortPCA %s: %s -------------- \n', date, fn)

%-----------------------
% Check inputs
if isempty(dir(fn))
    error('Invalid input file name.')
end
if (nargin<2)||(isempty(flims))
    nt_full = tiff_frames(fn);
    flims = [1,nt_full];
end
if nargin<6
    badframes = [];
end

useframes = setdiff((flims(1):flims(2)), badframes);
nt = length(useframes);

if nargin<3 || isempty(nPCs)
    nPCs = min(150, nt);
end
if nargin<4 || isempty(dsamp)
    dsamp = [1,1];
end
if nargin<5 || isempty(outputdir)
    outputdir = [pwd,'/cellsort_preprocessed_data/'];
end
if isempty(dir(outputdir))
    mkdir(pwd, '/cellsort_preprocessed_data/')
end
if outputdir(end)~='/';
    outputdir = [outputdir, '/'];
end

[~, fname] = fileparts(fn);
if isempty(badframes)
    fnmat = [outputdir, fname, '_',num2str(flims(1)),',',num2str(flims(2)), '_', date,'.mat'];
else
    fnmat = [outputdir, fname, '_',num2str(flims(1)),',',num2str(flims(2)),'_selframes_', date,'.mat'];
end
if ~isempty(dir(fnmat))
    fprintf('CELLSORT: Movie %s already processed;', ...
        fn)
    forceload = input(' Re-load data? [0-no/1-yes] ');
    if isempty(forceload) || forceload==0
        load(fnmat)
        return
    end
end

[pixw,pixh] = size(imread(fn,1));
npix = pixw*pixh;

fprintf('   %d pixels x %d time frames;', npix, nt)

% Create covariance matrix
if nt < npix 
    fprintf(' using temporal covariance matrix.\n')
    [covmat, mov, movm, movtm] = create_tcov(fn, pixw, pixh, useframes, nt, dsamp);
else
    fprintf(' using spatial covariance matrix.\n')
    [covmat, mov, movm, movtm] = create_xcov(fn, pixw, pixh, useframes, nt, dsamp);
end

covtrace = trace(covmat) / npix;
movm = reshape(movm, pixw, pixh);

if nt < npix 
    % Perform SVD on temporal covariance
    [mixedsig, CovEvals] = cellsort_svd(covmat, nPCs, nt, npix);

    % Load the other set of principal components
    [mixedfilters] = reload_moviedata(pixw*pixh, mov, mixedsig, CovEvals);
else
    % Perform SVD on spatial components
    [mixedfilters, CovEvals] = cellsort_svd(covmat, nPCs, nt, npix);
    mixedfilters = mixedfilters' * npix;
    
    % Load the other set of principal components
    [mixedsig] = reload_moviedata(nt, mov', mixedfilters', CovEvals);
    mixedsig = mixedsig' / npix^2;
end
mixedfilters = reshape(mixedfilters, pixw,pixh,nPCs);

%------------
% Save the output data
save(fnmat,'mixedfilters','CovEvals','mixedsig', ...
    'movm','movtm','covtrace')
fprintf(' CellsortPCA: saving data and exiting; ')
toc

    function [covmat, mov, movm, movtm] = create_xcov(fn, pixw, pixh, useframes, nt, dsamp)
        %-----------------------
        % Load movie data to compute the spatial covariance matrix

        npix1 = pixw*pixh;

        % Downsampling
        if length(dsamp)==1
            dsamp_time = dsamp(1);
            dsamp_space = 1;
        else
            dsamp_time = dsamp(1);
            dsamp_space = dsamp(2); % Spatial downsample
        end

        if (dsamp_space==1)
            mov = zeros(pixw, pixh, nt);
            for jjind=1:length(useframes)
                jj = useframes(jjind);
                mov(:,:,jjind) = imread(fn,jj);
                if mod(jjind,500)==1
                    fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
                    toc
                end
            end
        else
            [pixw_dsamp,pixh_dsamp] = size(imresize( imread(fn,1), 1/dsamp_space, 'bilinear' ));
            mov = zeros(pixw_dsamp, pixh_dsamp, nt);
            for jjind=1:length(useframes)
                jj = useframes(jjind);
                mov(:,:,jjind) = imresize( imread(fn,jj), 1/dsamp_space, 'bilinear' );
                if mod(jjind,500)==1
                    fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
                    toc
                end
            end
        end

        fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
        toc
        mov = reshape(mov, npix1, nt);

        % DFoF normalization of each pixel
        movm = mean(mov,2); % Average over time
        movmzero = (movm==0);
        movm(movmzero) = 1;
        mov = mov ./ (movm * ones(1,nt)) - 1; % Compute Delta F/F
        mov(movmzero, :) = 0;

        if dsamp_time>1
            mov = filter(ones(dsamp,1)/dsamp, 1, mov, [], 2);
            mov = downsample(mov', dsamp)';
        end

        movtm = mean(mov,1); % Average over space
        mov = mov - ones(size(mov,1),1)*movtm; 

        covmat = (mov*mov')/size(mov,2);
        covmat = covmat * size(mov,2)/size(mov,1); % Rescale to gree with create_tcov
        toc
    end

    function [covmat, mov, movm, movtm] = create_tcov(fn, pixw, pixh, useframes, nt, dsamp)
        %-----------------------
        % Load movie data to compute the temporal covariance matrix
        npix1 = pixw*pixh;

        % Downsampling
        if length(dsamp)==1
            dsamp_time = dsamp(1);
            dsamp_space = 1;
        else
            dsamp_time = dsamp(1);
            dsamp_space = dsamp(2); % Spatial downsample
        end

        if (dsamp_space==1)
            mov = zeros(pixw, pixh, nt);
            for jjind=1:length(useframes)
                jj = useframes(jjind);
                mov(:,:,jjind) = imread(fn,jj);
                if mod(jjind,500)==1
                    fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
                    toc
                end
            end
        else
            [pixw_dsamp,pixh_dsamp] = size(imresize( imread(fn,1), 1/dsamp_space, 'bilinear' ));
            mov = zeros(pixw_dsamp, pixh_dsamp, nt);
            for jjind=1:length(useframes)
                jj = useframes(jjind);
                mov(:,:,jjind) = imresize( imread(fn,jj), 1/dsamp_space, 'bilinear' );
                if mod(jjind,500)==1
                    fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
                    toc
                end
            end
        end

        fprintf(' Read frame %4.0f out of %4.0f; ', jjind, nt)
        toc
        mov = reshape(mov, npix1, nt);

        % DFoF normalization of each pixel
        movm = mean(mov,2); % Average over time
        movmzero = (movm==0); % Avoid dividing by zero
        movm(movmzero) = 1;
        mov = mov ./ (movm * ones(1,nt)) - 1;
        mov(movmzero, :) = 0;

        if dsamp_time>1
            mov = filter(ones(dsamp,1)/dsamp, 1, mov, [], 2);
            mov = downsample(mov', dsamp)';
        end

        c1 = (mov'*mov)/npix1;
        movtm = mean(mov,1); % Average over space
        covmat = c1 - movtm'*movtm;
        clear c1
    end

    function [mixedsig, CovEvals, percentvar] = cellsort_svd(covmat, nPCs, nt, npix1)
        %-----------------------
        % Perform SVD

        covtrace1 = trace(covmat) / npix1;

        opts.disp = 0;
        opts.issym = 'true';
        if nPCs<size(covmat,1)
            [mixedsig, CovEvals] = eigs(covmat, nPCs, 'LM', opts);  % pca_mixedsig are the temporal signals, mixedsig
        else
            [mixedsig, CovEvals] = eig(covmat);
            CovEvals = diag( sort(diag(CovEvals), 'descend'));
            nPCs = size(CovEvals,1);
        end
        CovEvals = diag(CovEvals);
        if nnz(CovEvals<=0)
            nPCs = nPCs - nnz(CovEvals<=0);
            fprintf(['Throwing out ',num2str(nnz(CovEvals<0)),' negative eigenvalues; new # of PCs = ',num2str(nPCs),'. \n']);
            mixedsig = mixedsig(:,CovEvals>0);
            CovEvals = CovEvals(CovEvals>0);
        end

        mixedsig = mixedsig' * nt;
        CovEvals = CovEvals / npix1;

        percentvar = 100*sum(CovEvals)/covtrace1;
        fprintf([' First ',num2str(nPCs),' PCs contain ',num2str(percentvar,3),'%% of the variance.\n'])
    end

    function [mixedfilters] = reload_moviedata(npix1, mov, mixedsig, CovEvals)
        %-----------------------
        % Re-load movie data
        nPCs1 = size(mixedsig,1);

        Sinv = inv(diag(CovEvals.^(1/2)));

        movtm1 = mean(mov,1); % Average over space
        movuse = mov - ones(npix1,1) * movtm1;
        mixedfilters = reshape(movuse * mixedsig' * Sinv, npix1, nPCs1);
    end

    function j = tiff_frames(fn)
        %
        % n = tiff_frames(filename)
        %
        % Returns the number of slices in a TIFF stack.
        %
        % Modified April 9, 2013 for compatibility with MATLAB 2012b

        j = length(imfinfo(fn));
    end
end