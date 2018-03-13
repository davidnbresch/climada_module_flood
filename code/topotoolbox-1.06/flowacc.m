function A = flowacc(M,W0,RR)

% calculate flow accumulation/upslope area from flow direction matrix
%
% Syntax
%
%     A = flowacc(M,siz)
%     A = flowacc(M,W0)
%     A = flowacc(M,[],RR)
%     A = flowacc(M,W0,RR)
%
% Description
%
%     flowacc returns the flow accumulation (or upslope area) based on the
%     flowdirection matrix generated by flowdir or flowdir_single.
%
% Input
%
%     M         sparse flow direction matrix (single and multi)
%     siz       two element vector with size of each dimension of the dem
%               the flow direction matrix is based on
%     W0        inflow rate
%               a matrix same size as the DEM from which the flow direction
%               matrix was derived. Elements in W0 contain the initial
%               amount of water in each cell. By default, W0 is
%               ones(size(dem)).
%     RR        runoff ratio
%               a matrix same size as the DEM from which the flow direction
%               matrix was derived. Elements in RR contain runoff 
%               ratios in each cell, which is defined as the ratio between 
%               water entering and leaving each cell. [RR >= 0 and <=1]
%
% Output
% 
%     A         flow accumulation/upslope area matrix. Values in A
%               indicate the number of cells draining into each cell.
%
%
% Example
%
%     load exampleDEM
%     M = flowdir(X,Y,dem,'exponent',5);
%     A = flowacc(M,size(dem));
%     surf(X,Y,dem,log(A)) 
%
% See also: FLOWDIR, CROSSFLATS, UPSLOPESTATS
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 15. March, 2009


% check input
error(nargchk(2, 3, nargin))

if nargin == 2;

    if numel(W0) == 2;
        siz = W0;
        nrc = prod(siz);
        W0  = ones(nrc,1);
    else    
        nrc = numel(W0);
        siz = size(W0);
        W0  = W0(:);
    end

else
    siz = size(RR);
    nrc = numel(RR);
    if isempty(W0);
        W0 = ones(nrc,1);
    else
        if ~isequal(siz,size(W0));
            error('TopoToolbox:incorrectinput',...
                  'W0 and RR must have same size')
        end
        W0 = W0(:);
    end    
end
    

% does the flow direction matrix correspond to the DEM
if (nrc~=size(M,1)) || (nrc~=size(M,2));
    error('TopoToolbox:incorrectinput',...
          'M must be a square, sparse matrix with size [prod(siz) prod(siz)]')
end

% if runoff ratios are supplied, check whether they range between 0 and 1
% and renormalize M.
if nargin==3;
    if any(RR(:)<0) || any(RR(:)>1)
        error('TopoToolbox:incorrectinput',...
          'values in RR must range between 0 and 1')
    end
    meth = 1; %%% needs to be checked
    % Method 1: Values in A will be inflow rates in each cell
    if meth == 1;
    M = spdiags(RR(:),0,nrc,nrc)*M;
    I = speye(nrc);
    % Method 2: Values in A will be outflow rates of each cell
    else
        RR = 1./RR;
        I  = spdiags(RR(:),0,nrc,nrc);
    end
else
    I = speye(nrc);
    
end
    
% solve flow accumulation equation
A = (I-M')\W0;
% reshape array
A = reshape(A,siz);


