function [ica_segments, segmentlabel, segcentroid] = CellsortSegmentation(ica_filters, smwidth, thresh, arealims, plotting)
% [ica_segments, segmentlabel, segcentroid] = CellsortSegmentation(ica_filters, smwidth, thresh, arealims, plotting)
%
%CellsortSegmentation
% Segment spatial filters derived by ICA
%
% Inputs:
%     ica_filters - X x Y x nIC matrix of ICA spatial filters
%     smwidth - standard deviation of Gaussian smoothing kernel (pixels)
%     thresh - threshold for spatial filters (standard deviations)
%     arealims - 2-element vector specifying the minimum and maximum area
%     (in pixels) of segments to be retained; if only one element is
%     specified, use this as the minimum area
%     plotting - [0,1] whether or not to show filters
%
% Outputs:
%     ica_segments - segmented spatial filters
%     segmentabel - indices of the ICA filters from which each segment was derived
%     segcentroid - X,Y centroid, in pixels, of each segment
%
% Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, 2009
% Email: eran@post.harvard.edu, mschnitz@stanford.edu
%

tic
fprintf('-------------- CellsortSegmentation %s -------------- \n', date)

if (nargin<3)||isempty(thresh)
    thresh = 2;
end
if (nargin<4)||isempty(arealims)
    arealims = 200;
end
if (nargin<5)||isempty(plotting)
    plotting = 0;
end
[nic,pixw,pixh] = size(ica_filters);

ica_filtersorig = ica_filters / abs(std(ica_filters(:)));
ica_filters = (ica_filters - mean(ica_filters(:)))/abs(std(ica_filters(:)));
if smwidth>0
    % Smooth mixing filter with a Gaussian of s.d. smwidth pixels
    smrange = max(5,3*smwidth);
    [x,y] = meshgrid([-smrange:smrange]);

    smy = 1; smx = 1;
    ica_filtersfilt = exp(-((x/smx).^2 + (y/smy).^2)/(2*smwidth^2));
    
    ica_filtersfilt = ica_filtersfilt/sum(ica_filtersfilt(:));
    ica_filtersbw = false(pixw,pixh,nic);
    tic
    for j = 1:size(ica_filters,1)
        ica_filtersuse = ica_filters(j,:,:);
        ica_filtersuse = (ica_filtersuse - mean(ica_filtersuse(:)))/abs(std(ica_filtersuse(:)));
        ica_filtersbw(:,:,j) = (imfilter(ica_filtersuse, ica_filtersfilt, 'replicate', 'same') > thresh);
    end
else
    ica_filtersbw = (permute(ica_filters,[2,3,1]) > thresh);
    ica_filtersfilt = 1;
end

tic
if plotting
    clf
    set(gcf,'Color','w')
    colormap(gray)
    
    subplot(223)
    imagesc(squeeze(sum(ica_filters,1)))
    axis image off
    hold on
end
ica_filterslabel = [];
ica_segments = [];
k=0;
L=[];
segmentlabel = [];
segcentroid = [];
[x,y] = meshgrid([1:pixh], [1:pixw]);
for j = 1:nic
    % Label contiguous components
    L = bwlabel(ica_filtersbw(:,:,j), 4);
    Lu = 1:max(L(:));

    % Delete small components
    Larea = struct2array(regionprops(L, 'area'));
    Lcent = regionprops(L, 'Centroid');
    
    if length(arealims)==2
        Lbig = Lu( (Larea >= arealims(1))&(Larea <= arealims(2)));
        Lsmall = Lu((Larea < arealims(1))|(Larea > arealims(2)));
    else
        Lbig = Lu(Larea >= arealims(1));
        Lsmall = Lu(Larea < arealims(1));
    end
    
    L(ismember(L,Lsmall)) = 0;   
        
    for jj = 1:length(Lbig)
        segcentroid(jj+k,:) = Lcent(Lbig(jj)).Centroid;
    end
    
    ica_filtersuse = squeeze(ica_filtersorig(j,:,:));
    for jj = 1:length(Lbig)
        ica_segments(jj+k,:,:) = ica_filtersuse .* ( 0*(L==0) + (L==Lbig(jj)) );  % Exclude background
    end
    
    if plotting && ~isempty(Lbig)
        if smwidth>0
            subplot(2,2,2)
            ica_filtersuse = squeeze(ica_filters(j,:,:));
            ica_filtersuse = (ica_filtersuse - mean(ica_filtersuse(:)))/abs(std(ica_filtersuse(:)));
            imagesc(imfilter((ica_filtersuse), ica_filtersfilt, 'replicate', 'same'),[-1,4])
            hold on
            contour(imfilter((ica_filtersuse), ica_filtersfilt, 'replicate', 'same'), [1,1]*thresh, 'k')
            hold off
            hc = colorbar('Position',[0.9189    0.6331    0.0331    0.2253]);
            ylabel(hc,'Std. dev.')
            title(['IC ',num2str(j),' smoothed'])
            axis image off

            subplot(2,2,1)
        else
            subplot(211)
        end
        imagesc(squeeze(ica_filters(j,:,:)))
        title(['IC ',num2str(j),' original'])
        axis image off

        colord = lines(k+length(Lbig));
        for jj = 1:length(Lbig)
            subplot(223)
            contour(ica_filtersbw(:,:,j), [1,1]*0.5, 'color',colord(jj+k,:),'linewidth',2)
            hold on
            text(segcentroid(jj+k,1), segcentroid(jj+k,2), num2str(jj+k), 'horizontalalignment','c', 'verticalalignment','m')
            set(gca, 'ydir','reverse','tickdir','out')
            axis image
            xlim([0,pixw]); ylim([0,pixh])
            
            subplot(224)
            imagesc(squeeze(ica_segments(jj+k,:,:)))
            hold on
            plot(segcentroid(jj+k,1), segcentroid(jj+k,2), 'bo')
            hold off
            axis image off
            title(['Segment ',num2str(jj+k)])
            drawnow
        end
    end
    k = size(ica_segments,1);
end
toc
