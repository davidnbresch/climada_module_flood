function C = acv(dem)

% Anisotropic coefficient of variation (ACV) 
%
% Syntax
%
%     C = acv(dem)
%
% Description
% 
%     Anisotropic coefficient of variation
%
% Input
%
%     dem       digital elevation model
%
% Output
%
%     C         Anisotropic coefficient of variation (ACV)
%
% Example
%     
%     % highlight ridges and valleys by substracting
%     % the original digital elevation model from
%     % the smoothed one.
%
%     load exampleDEM
%     C  = acv(dem);
%     surf(X,Y,dem,C);
%
% References
%
%     Olaya, V. 2009: Basic land-surface parameters. In: Geomorphometry. 
%     Concepts, Software, Applications, Hengl, T. & Reuter, H. I. (Eds.),
%     Elsevier, 33, 141-169.
% 
% See also: CONV2, FILTER2, BWDIST
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)
% Date: 6. July, 2011


siz = size(dem);

dem     = padarray(dem,[2 2],nan);
inanpad = isnan(dem);
[L,L]   = bwdist(~inanpad); %#ok<ASGLU>
dem     = dem(L);


k       = [ 1  0  1  0  1; ...
            0  0  0  0  0; ...
            1  0  0  0 -1; ...
            0  0  0  0  0; ...
           -1  0 -1  0 -1];
       
dz_AVG  = conv2(dem,k,'valid')/4;

F = {     [ 0  0  0  0  0; ...
            0  0  0  0  0; ...
            1  0  0  0 -1; ...
            0  0  0  0  0; ...
            0  0  0  0  0];...
            ...
          [ 1  0  0  0  0; ...
            0  0  0  0  0; ...
            0  0  0  0  0; ...
            0  0  0  0  0; ...
            0  0  0  0 -1];
            ...
          [ 0  0 -1  0  0; ...
            0  0  0  0  0; ...
            0  0  0  0  0; ...
            0  0  0  0  0; ...
            0  0 -1  0  0];
            ...
          [ 0  0  0  0 -1; ...
            0  0  0  0  0; ...
            0  0  0  0  0; ...
            0  0  0  0  0; ...
            1  0  0  0  0];...
     };
            
ACV = zeros(siz);
 
for r = 1:4;
    ACV = ACV + (conv2(dem,F{r},'valid') - dz_AVG).^2;
end

dem = dem(2:end-1,2:end-1);

F = {      [ 0  0  0;...
             1  0 -1;...
             0  0  0];...
             ...
           [ 1  0  0;...
             0  0  0;...
             0  0 -1];...
             ...
           [ 0  1  0;...
             0  0  0;...
             0 -1  0];...
             ...
           [ 0  0  1;...
             0  0  0;...
            -1  0  0];...
     };

for r = 1:4;
    ACV = ACV + (conv2(dem,F{r},'valid') - dz_AVG).^2;
end

dz_AVG = max(abs(dz_AVG),0.001);

C = log(1 + sqrt(ACV./8)./dz_AVG);
             
             
           
           
    




