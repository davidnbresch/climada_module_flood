function [L,outlets] = drainagebasins(dem,siz,order)

% segment a digital elevation model in drain basins/catchments
%
% Syntax
%
%     L = drainagebasins(dem)
%     L = drainagebasins(M,siz)
%     L = drainagebasins(M,S,order)
%     [L,outlets] = drainagebasins(M,...)
%
%
% Description
%
%     drainagebasins determines drainage basin affiliation of cells in 
%     a digital elevation model (DEM). When drainagebasins is supplied 
%     with a elevation matrix only, the function acts as a simple wrapper 
%     for the function watershed shipped with the Image Processing Toolbox. 
%     drainagebasins thereby allows for NaNs in the DEM. When applying 
%     the watershed algorithm, pixels on divides are returned as zero.
%
%     When drainagebasins is supplied with the flow direction matrix, all
%     cells have a clear affiliation to a drainage basin. Note that only
%     the single flow direction matrix can be supplied. 
%
%     drainagebasins(M,S,order) allows you to delineate drainagebasins for
%     specific (Strahler) order of streams. S is the first output returned
%     by the function streamorder and order is a scalar for the specific
%     stream order. 
%
%     When called with the flow direction matrix as first input argument
%     the second output contains the linear indices of drainage basin
%     outlet cells.
%
% Input
%
%     dem       digital elevation model
%     M         flow direction matrix
%     siz       size of the dem [nrrows nrcols]
%     S         stream order raster (as generated by streamorder)
%     order     stream order (scalar integer)
% 
% Output
%
%     L         output label matrix
%     outlets   linear indices of drainage basin outlet cells
%
% Example 1
%
%     load exampleDEM
%     M = flowdir_single(dem);
%     L = drainagebasins(M,size(dem));
%     imagesc(X(1,:),Y(:,2),L); axis image; axis xy
%     hold on
%     gplot(M,[X(:) Y(:)],'k');
%
% Example 2
%
%     load exampleDEM
%     [A,M] = ezflowacc(X,Y,dem,'type','single');
%     % let's simply assume that channels start where
%     % A is larger than 100;
%     W = A>100;
%     S = streamorder(M,W);
%     DB = drainagebasins(M,S,2);
%     imagesc(X(1,:),Y(:,2),DB); axis image; axis xy
%     subplot(1,2,1)
%     imagesc(X(1,:),Y(:,2),S); axis image; axis xy
%     colorbar
%     title('Strahler stream order')
%     subplot(1,2,2)
%     imagesc(X(1,:),Y(:,2),shufflelabel(DB)); axis image; axis xy
%     colorbar
%     title('order 2 drainage basins')
%
%
% See also: FLOWDIR_SINGLE, STREAMORDER, REGIONPROPS, DRAINAGEDENSITY
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 16. September, 2009



if ~issparse(dem) 
    % use built in watershed function
    I = isnan(dem);
    dem(I) = -inf;
    L = watershed(dem,8);

    L(I) = nan;
    
    if nargout == 2;
        outlets = [];
        warning('TopoToolbox:Output',...
                ['the second output argument is empty. Identification of\n'...
                 'outlets is only supported when the function is called\n'...
                 'with the single flow direction matrix as first input argument']);
    end
    
else
    % use own function 
    M = dem;
    if numel(siz) == 2;
        flagstreamorder = false;
        nrc = prod(siz);
    else
        flagstreamorder = true;
        nrc = numel(siz);
        S   = siz;
        siz = size(S);
        % check stream order
        if ~isscalar(order) || round(order)~=order;
            error('TopoToolbox:incorrectinput',...
                  'order must be a scalar integer')
        end
    end
        
    
    % do we have a multiple flow direction matrix?
    M = spones(M);
    if any(sum(M,2)>1);
        error('TopoToolbox:incorrectinput',...
              'single flow direction must be supplied')
    end
    
    % find indices that don't route water but receive water
    if ~flagstreamorder
        NG = full(sum(M,2) == 0) & (full(sum(M,1)') ~= 0); 
        
        if nargout==2;
            outlets = find(NG);
        end
        
        NG = cumsum(NG).*NG;
    else
        NG = S == order;
        NG = NG(:);
        if any(NG)
        else
            error('TopoToolbox:incorrectinput',...
            'there are no streams of such order.')
        end

        % find basin outlets
        NG = M*NG==0 & NG;
        
        if nargout==2;
            outlets = find(NG);
        end
        
        NG = cumsum(NG).*NG;
        
        
    end
        
    
    L = (speye(nrc)-M)\NG;
    
    % set values in L to zero that occur only once
    I = full(sum(M,2))==0 & full(sum(M,1))'==0;
    L(I) = 0;
    [l,l,l] = unique(L(~I));
    L(~I) = l;
    
    L = reshape(L,siz);
    
end

% 


