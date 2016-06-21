function [PCuse] = CellsortChoosePCs(fn, mixedfilters)
% [PCuse] = CellsortChoosePCs(fn, mixedfilters)
%
% Allows the user to select which principal components will be kept
% following dimensional reduction.
%
% Inputs:
%   fn - movie file name. Must be in TIFF format.
%   mixedfilters - N x X matrix of N spatial signal mixtures sampled at X
%   spatial points.
%
% Outputs:
%   PCuse - vector of indices of the PCs to be kept for dimensional
%   reduction
%
% Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, 2009
% Email: eran@post.harvard.edu, mschnitz@stanford.edu
%

fprintf('-------------- CellsortChoosePCs %s -------------- \n', date)

[pixw,pixh] = size(imread(fn,1));

npcs = 20; % Number of PCs to display concurrently
currpcs = [1:npcs];
PCf = [];
while isempty(PCf)
    showpcs(currpcs, mixedfilters, pixw, pixh)
    yl = ylim;
    xl = xlim;
    set(gca,'Units','pixels')
    title(['Choose first PC; showing PCs [',num2str(currpcs(1)),':',num2str(currpcs(end)),']'])
    PCf = input('Number of first PC to retain, Klow (''b/f'' to scroll backwards/forwards)): ','s');
    if PCf=='b'
        currpcs = currpcs - min(npcs,currpcs(1)-1);
        PCf = [];
    elseif (PCf=='f')
        currpcs = currpcs+npcs;
        if nnz(currpcs>size(mixedfilters,2))
            currpcs = [-npcs+1:0]+size(mixedfilters,2);
            fprintf('Reached end of stored PCs.\n')
        end
        PCf = [];
    else
        PCf = str2num(PCf);
    end
end
PCl=[];
currpcs = [PCf:PCf+npcs-1];
while isempty(PCl)
    showpcs(currpcs, mixedfilters, pixw, pixh)
    title(['Choose last PC; showing PCs [',num2str(currpcs(1)),':',num2str(currpcs(end)),']'])
    PCl = input('Number of last PC to retain, Khigh (''b/f'' to scroll backwards/forwards): ','s');
    if PCl=='b'
        currpcs = currpcs - min(npcs,currpcs(1)-1);
        PCl = [];
    elseif (PCl=='f')
        currpcs = currpcs+npcs;
        if nnz(currpcs>size(mixedfilters,2))
            currpcs = [-npcs+1:0]+size(mixedfilters,2);
            fprintf('Reached end of stored PCs.\n')
        end
        PCl = [];
    else
        PCl = str2num(PCl);
    end
end
currpcs = [PCf:PCl];
PCbad=[];
showpcs(currpcs, mixedfilters, pixw, pixh)

PCuse = setdiff(currpcs, PCbad);
showpcs(PCuse, mixedfilters, pixw, pixh)

fprintf('  Retaining PCs in the range [Klow - Khigh] = [%d - %d].\n', PCf,PCl)

function showpcs(usepcs, Efull, pixw, pixh)

if nargin<3
    fprintf('Assuming movie frames are square.\n')
    pixw = sqrt(size(Efull,1));
    pixh = sqrt(size(Efull,1));
end
if isempty(usepcs)
    usepcs = [1:size(Efull,2)];
end

if ndims(Efull)>=3
    Efull = reshape(Efull, pixw*pixh, []);
end
for j=usepcs
    Efull(:,j) = zscore(Efull(:,j));
end
pcs = reshape(Efull(:,usepcs), pixw, pixh, []);
pcs = permute(pcs, [1, 2, 4, 3]);
montage(pcs)
colormap(hot)
axis on
xl = xlim;
yl = ylim;
nw = ceil(xl(2)/pixh)-1;
nh = ceil(yl(2)/pixw)-1;
set(gca,'YTick',[pixw:pixw:yl(2)],'YTickLabel',  num2str(usepcs(min([0:nh]*nw+1, length(usepcs)))'), ...
    'XTick',[pixh:pixh:xl(2)], ...
    'XTickLabel',num2str(usepcs([(nh-1)*nw+1:length(usepcs)])'), 'XAxisLocation','bottom','LineWidth',2)
grid on
formataxes
caxis([-1,1]*7)

function formataxes

set(gca,'FontSize',12,'FontWeight','bold','FontName','Helvetica','LineWidth',2,'TickLength',[1,1]*.02,'tickdir','out')
set(gcf,'Color','w','PaperPositionMode','auto')
