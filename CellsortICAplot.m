function CellsortICAplot(mode, ica_filters, ica_sig, f0, tlims, dt, ratebin, plottype, ICuse, spt, spc)
% CellsortICAplot(mode, ica_filters, ica_sig, f0, tlims, dt, ratebin, plottype, ICuse, spt, spc)
%
% Display the results of ICA analysis in the form of paired spatial filters
% and signal time courses
%
% Inputs:
%     mode - 'series' shows each spatial filter separately; 'contour'
%     displays a single plot with contours for all spatial filters
%     superimposed on the mean fluorescence image
%     ica_filters - nIC x X x Y array of ICA spatial filters
%     ica_sig - nIC x T matrix of ICA temporal signals
%     f0 - mean fluorescence image
%     tlims - 2-element vector specifying the range of times to be displayed
%     dt - time step corresponding to individual movie time frames
%     ratebin - size of time bins for spike rate computation
%     plottype - type of spike plot to use:
%         plottype = 1: plot cellular signals
%         plottype = 2: plot cellular signals together with spikes
%         plottype = 3: plot spikes only
%         plottype = 4: plot spike rate over time
%     ICuse - vector of indices of cells to be plotted
%     spt - vector of spike times (needed for plottype > 1)
%     spc - vector of spike indices (which cell) (needed for plottype > 1)
%
% Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, 2009
% Email: eran@post.harvard.edu, mschnitz@stanford.edu
%

colord=[         0         0    1.0000
    0    0.4000         0
    1.0000         0         0
    0    0.7500    0.7500
    0.7500         0    0.7500
    0.8, 0.5, 0
    0         0    0.5
    0         0.85      0];

% Check input arguments
nIC = size(ica_sig,1);
if nargin<5 || isempty(tlims)
    tlims = [0, size(ica_sig,2)*dt]; % seconds
end
if nargin<8 || isempty(plottype)
    plottype = 1;
end
if nargin<9 || isempty(ICuse)
    ICuse = [1:nIC];
end
if size(ICuse,2)==1
    ICuse = ICuse';
end

% Reshape the filters
[pixw,pixh] = size(f0);
if size(ica_filters,1)==nIC
    ica_filters = reshape(ica_filters,nIC,pixw,pixh);
elseif size(ica_filters,2)==nIC
    ica_filters = reshape(ica_filters,nIC,pixw,pixh);
end

switch mode
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'series'}
        colmax = 20; % Maximum # of ICs in one column
        ncols = ceil(length(ICuse)/colmax);
        if plottype==4
            ncols=1;
        end
        nrows = ceil(length(ICuse)/ncols);

        if size(ica_filters(:,:,1))==size(f0(:,:,1))
            ica_filters = permute(ica_filters,[3,1,2]);
        end

        subplot(1,3*ncols,[2:3])
        tlims(2) = min(tlims(2),size(ica_sig,2)*dt);
        tlims(1) = max(tlims(1),0);

        clf
        f_pos = get(gcf,'Position');
        f_pos(4) = max([f_pos(4),500,50*nrows]);
        f_pos(3) = max(400*ncols,0.9*f_pos(4));

        colormap(hot)
        colord=get(gca,'ColorOrder');
        ll=0;
        filtax = [];
        if ~isempty(ica_filters)
            for k=0:ncols-1
                jj=3*k;
                nrows_curr = min(nrows,length(ICuse)-k*nrows);
                for j=1:nrows_curr
                    filtax= [filtax,subplot(nrows_curr + (plottype==4),3*ncols, jj+1)];
                    jj=jj+3*ncols;
                    ll=ll+1;
                    imagesc(squeeze(ica_filters(ICuse(ll),:,:)))
                    axis image tight off
                end
            end
        end

        ax = [];
        for j=0:ncols-1
            ax(j+1)=subplot(1,3*ncols,3*j+[2:3]);
            ICuseuse = ICuse([1+j*nrows:min(length(ICuse),(j+1)*nrows)]);
            if plottype<=2
                complot(ica_sig, ICuseuse, dt)
            end
            formataxes
            if plottype>=2
                [spcICuse,spcord] = ismember(spc,ICuseuse);
                spc_curr = spcord(spcICuse&(spt<=tlims(2))&(spt>=tlims(1)));
                spt_curr = spt(spcICuse&(spt<=tlims(2))&(spt>=tlims(1)));
                alpha = diff(ylim)/length(ICuseuse);
                switch plottype
                    case 2
                        hold on
                        scatter(spt_curr,-alpha*(spc_curr-1)+2,20,colord(mod(spc_curr-1,size(colord,1))+1,:),'o','filled')
                        hold off
                    case 3
                        cla
                        scatter(spt_curr,-alpha*(spc_curr-1)+0.4*alpha,20,colord(mod(spc_curr-1,size(colord,1))+1,:),'o','filled')
                        yticks=sort(-alpha*([1:length(ICuseuse)]-1)+0.4*alpha);
                        yl = sort(-alpha*[length(ICuseuse)-1,-1]);
                        set(gca,'YTick',yticks)
                        set(gca,'YTicklabel',num2str(fliplr(ICuseuse)'),'YLim',yl)
                        set(gca,'YLim',sort(-alpha*[length(ICuseuse)-1,-0.6]))

                    case 4
                        ax4=[];
                        tvec = [ratebin/2:ratebin:size(ica_sig,2)*dt];
                        binsize = ratebin*ones(size(tvec));
                        binsize(end) = size(ica_sig,2)*dt-tvec(end-1)-ratebin/2;
                        fprintf('IC mean rate\n')
                        sprm = []; sprs = [];
                        for jj=1:length(ICuseuse)
                            sptj = spt_curr(spc_curr==jj);
                            spn = hist(sptj,tvec);
                            spr = spn./binsize;
                            sperr = sqrt(spn)./binsize;

                            ax4(jj) = subplot(length(ICuseuse)+1,3,[3*jj-1,3*jj]);
                            errorbar(tvec,spr,sperr,'Color',colord(mod(jj-1,size(colord,1))+1,:),'LineWidth',1)
                            hold on
                            bar(tvec,spr)
                            hold off
                            axis tight
                            formataxes
                            ylabel({['IC',int2str(ICuseuse(jj))],[num2str(mean(spr),2),'']})
                            sprm(jj) = sum(spn)/sum(binsize);
                            fprintf('%3.3f\n', sprm(jj))
                            box off

                        end
                        histax = subplot(length(ICuseuse)+1,3, 3*length(ICuseuse)+1);
                        hist(sprm,[0:0.1:1])
                        formataxes
                        xlabel('Spike rate (Hz)')
                        ylabel('# Cells')
                        xlim([0,1])

                        spn = hist(spt_curr,tvec);
                        spr = spn./(binsize*length(ICuseuse));
                        sperr = sqrt(spn)./(binsize*length(ICuseuse));  %Standard error for measurement of mean from Poisson process

                        ax4(jj+1) = subplot(length(ICuseuse)+1,3,3*length(ICuseuse)+[2,3]);
                        errorbar(tvec,spr,sperr,'k','LineWidth',3)
                        formataxes
                        set(ax4,'XLim',tlims)
                        set(ax4(1:jj-1),'XTick',[])
                        ylabel({'Pop. mean',[num2str(mean(spr),2),' Hz']},'FontAngle','i')
                        fprintf('Pop. ave\tPop. SD\n')
                        fprintf('%3.3f\t\t%3.3f\n', mean(sprm), std(sprm))
                        xlabel('Time (s)','FontAngle','i')
                        box off
                        ax(j+1) = ax4(1);
                    case 5
                        % Scatter plot of spikes, showing synchrony
                        cla
                        [nsp, nsp_t] = hist(spt_curr, unique(spt_curr));
                        nsp_t = nsp_t(nsp>1);
                        nsp = nsp(nsp>1);
                        for jj = 1:length(nsp_t)
                            plot( nsp_t(jj)*ones(nsp(jj),1), spc_curr(spt_curr==nsp_t(jj)), 'o-', ...
                                'Color', colord(mod(nsp(jj)-2,size(colord,1))+1,:), 'LineWidth', 2, 'MarkerSize',10)
                            hold on
                        end
                        scatter(spt_curr,-spc_curr,20,colord(mod(spc_curr-1,size(colord,1))+1,:),'o','filled')
                        hold off
                        yticks=sort([1:max(spc_curr)]);
                        yl = sort([max(spc_curr)+1,min(spc_curr)-1]);
                        set(gca,'YTick',yticks)
                        set(gca,'YTicklabel',num2str(fliplr(ICuse)'),'YLim',yl)
                end
            end
            formataxes
            xlabel('Time (s)')
            xlim(tlims)
            yl = ylim;
            drawnow
        end
        set(gcf,'Color','w','PaperPositionMode','auto')

        %%%%
        % Resize plots to appropriate size
        if (plottype<4)&(length(ICuse)>=3)
            bigpos = get(ax(1),'Position');
            aheight = 0.9*bigpos(4)/nrows;
            for k=1:length(filtax)
                axpos = get(filtax(k),'Position');
                axpos(3) = aheight;
                axpos(4) = aheight;
                set(filtax(k),'Position',axpos)
            end

            set(gcf,'Units','normalized')
            fpos = get(gcf,'Position');
            for j=1:ncols
                axpos = get(ax(j),'OuterPosition');
                filtpos = get(filtax(1+(j-1)*nrows),'Position');
                axpos(1) = filtpos(1) + filtpos(3)*1.1;
                set(ax(j),'OuterPosition',axpos,'ActivePositionProperty','outerposition')
                axpos = get(ax(j),'Position');
            end
            set(gcf,'Resize','on','Units','characters')
        end

        for j=1:ncols
            ax = [ax,axes('Position',get(ax(j),'Position'),'XAxisLocation','top','Color','none')];
            if plottype==4
                xt = get(ax4(end),'XTick');
            else
                xt = get(ax(j),'XTick');
            end
            xlim(tlims)
            formataxes
            set(gca,'YTick',[],'XTick',xt,'XTickLabel',num2str(xt'/dt, '%15.0f'))
            xlabel('Frame number')
            axes(ax(j))
            box on
        end
        linkaxes(ax,'xy')


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'contour'}
        f_pos = get(gcf,'Position');
        if f_pos(4)>f_pos(3)
            f_pos(4) = 0.5*f_pos(3);
            set(gcf,'Position',f_pos);
        end
        set(gcf,'Renderer','zbuffer','RendererMode','manual')

        subplot(1,2,2)

        clf
        colormap(gray)

        set(gcf,'DefaultAxesColorOrder',colord)
        subplot(1,2,1)

        crange = [min(min(min(ica_filters(ICuse,:,:)))),max(max(max(ica_filters(ICuse,:,:))))];
        contourlevel = crange(2) - diff(crange)*[1,1]*0.8;

        cla
        if ndims(f0)==2
            imagesc(f0)
        else
            image(f0)
        end
        cax = caxis;
        sigsm = 1;
        shading flat
        hold on
        for j=1:length(ICuse)
            ica_filtersuse = gaussblur(squeeze(ica_filters(ICuse(j),:,:)), sigsm);
            contour(ica_filtersuse, [1,1]*(mean(ica_filtersuse(:))+4*std(ica_filtersuse(:))), ...
                'Color',colord(mod(j-1,size(colord,1))+1,:),'LineWidth',2)
        end
        for j=1:length(ICuse)
            ica_filtersuse = gaussblur(squeeze(ica_filters(ICuse(j),:,:)), sigsm);

            % Write the number at the cell center
            [ypeak, xpeak] = find(ica_filtersuse == max(max(ica_filtersuse)),1);
            text(xpeak,ypeak,num2str(j), 'horizontalalignment','c','verticalalignment','m','color','y')
        end
        hold off
        caxis(cax)
        formataxes
        axis image tight off
        title('Avg of movie, with contours of ICs')

        ax = subplot(1,2,2);
        if plottype<=2
            complot(ica_sig, ICuse, dt)
        end
        set(gca,'ColorOrder',colord)
        if plottype>=2
            [spcICuse,spcord] = ismember(spc,ICuse);
            spc = spcord(spcICuse&(spt<=tlims(2))&(spt>=tlims(1)));
            spt = spt(spcICuse&(spt<=tlims(2))&(spt>=tlims(1)));
            alph = diff(ylim)/length(ICuse);
            switch plottype
                case 2
                    % Traces with spikes
                    hold on
                    scatter(spt,-alph*(spc-1)+0.4*alph,20,colord(mod(spc-1,size(colord,1))+1,:),'o','filled')
                    hold off
                case 3
                    % Spike raster only
                    cla
                    scatter(spt,-alph*(spc-1)+0.4*alph,20,colord(mod(spc-1,size(colord,1))+1,:),'o','filled')
                    yticks=sort(-alph*([1:length(ICuse)]-1)+0.4*alph);
                    yl = sort(-alph*[length(ICuse)-1,-1]);
                    set(gca,'YTick',yticks)
                    set(gca,'YTicklabel',num2str(fliplr(ICuse)'),'YLim',yl)
                    set(gca,'YLim',sort(-alph*[length(ICuse)-1,-0.6]))
                case 4
                    % Binned spike rate
                    tvec = [ratebin/2:ratebin:size(ica_sig,2)*dt];
                    binsize = ratebin*ones(size(tvec));
                    binsize(end) = size(ica_sig,2)*dt-tvec(end-1)-ratebin/2;
                    for j=1:length(ICuse)
                        sptj = spt(spc==j);
                        spn = hist(sptj,tvec);
                        spr(j,:) = spn./binsize;
                        sperr(j,:) = sqrt(spn)./binsize;
                    end
                    spn = hist(spt,tvec);
                    mspr = spn./(binsize*length(ICuse));
                    msperr = sqrt(spn)./(binsize*length(ICuse));

                    hold on
                    plot(tvec,mspr,'k','LineWidth',3);
                    patch([tvec(end:-1:1),tvec],[mspr(end:-1:1)+msperr(end:-1:1),mspr-msperr],0.5*[1,1,1],'FaceAlpha',1)
                    hold off
                    axis tight
                    formataxes
            end
        end
        formataxes
        xlim(tlims)
        xlabel('Time (s)','FontAngle','i')
        if plottype<=3
            ylabel('IC #','FontAngle','i')
        else
            ylabel('Spike rate (Hz)','FontAngle','i')
        end
        set(gcf,'Color','w','PaperPositionMode','auto')
        set(gca,'yticklabel',num2str(fliplr([1:length(ICuse)])'))

        axes('Position',get(ax,'Position'),'XAxisLocation','top','Color','none')
        xt = get(ax,'XTick');
        xlim(tlims)
        formataxes
        set(gca,'YTick',[], ...
            'XTick',xt,'XTickLabel',num2str(xt'/dt, '%15.0f'))
        xlabel('Frame number')
        axes(ax)
        box on

end

%%%%%%%%%%%%%%%%%%%%%
function complot(sig, ICuse, dt)

for i = 1:length(ICuse)
    zsig(i, :) = zscore(sig(ICuse(i),:));
end

alpha = mean(max(zsig')-min(zsig'));
if islogical(zsig)
    alpha = 1.5*alpha;
end

zsig2 = zsig;
for i = 1:size(ICuse,2)
    zsig2(i,:) = zsig(i,:) - alpha*(i-1)*ones(size(zsig(1,:)));
end

tvec = [1:size(zsig,2)]*dt;
if islogical(zsig)
    plot(tvec, zsig2','LineWidth',1)
else
    plot(tvec, zsig2','LineWidth',1)
end
axis tight

set(gca,'YTick',(-size(zsig,1)+1)*alpha:alpha:0);
set(gca,'YTicklabel',fliplr(ICuse));


function formataxes

set(gca,'FontSize',12,'FontWeight','bold','FontName','Helvetica','LineWidth',2,'TickLength',[1,1]*.02,'tickdir','out')
set(gcf,'Color','w','PaperPositionMode','auto')

function fout = gaussblur(fin, smpix)
%
% Blur an image with a Gaussian kernel of s.d. smpix
%

if ndims(fin)==2
    [x,y] = meshgrid([-ceil(3*smpix):ceil(3*smpix)]);
    smfilt = exp(-(x.^2+y.^2)/(2*smpix^2));
    smfilt = smfilt/sum(smfilt(:));

    fout = imfilter(fin, smfilt, 'replicate', 'same');
else
    [x,y] = meshgrid([-ceil(smpix):ceil(smpix)]);
    smfilt = exp(-(x.^2+y.^2)/(2*smpix^2));
    smfilt = smfilt/sum(smfilt(:));

    fout = imfilter(fin, smfilt, 'replicate', 'same');
end
