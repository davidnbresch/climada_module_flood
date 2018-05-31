function climada_probSet_ls(width,binORlin)

% Script to calibrate flow path parameters 
% MODULE:
%  flood
% NAME:
%  climada_probSet_ls
% PURPOSE:
%    
% CALLING SEQUENCE:
%   
% EXAMPLE:
%   
% INPUTS:
%   
% OPTIONAL INPUT PARAMETERS:
%   
% OUTPUTS:
%  
% MODIFICATION HISTORY:
% Thomas Rölli, thomasroelli@gmail.com, 20180531, init

global climada_global
if ~climada_init_vars, return; end

if ~exist('binORlin') binORlin = 'lin'; end

subS = shaperead('C:\Users\Simon Rölli\Desktop\data\calibration\subData\subS_2x3m.shp');
numrun = 10^6;

rmv = [subS.removed];
reachAngle = [subS.reachAngle];
reachAngle(rmv~=0) = nan;
clear rmv

edges_vmax = [1 4 8 11];
[count_phi,edges_phi] = histcounts(reachAngle,'BinWidth',width,'Normalization', 'probability');
numbin_vmax = numel(edges_vmax)-1;
numbin_phi = numel(edges_phi)-1;

%generate uniform random number within the ranges of the two histograms
%(vmax and phi)
rand_vmax = edges_vmax(1) + rand([1,numrun])*(edges_vmax(end)-edges_vmax(1));
rand_phi = edges_phi(1) + rand([1,numrun])*(edges_phi(end)-edges_phi(1));

%%
%calculate probability for random vmax
prob_vmax = zeros(size(rand_vmax));

%calculate probability
%out of edges = 0;
idx = rand_vmax<edges_vmax(1) | rand_vmax>edges_vmax(end);
prob_vmax(idx) = 0;
%constant probability between 4 and 8 = 1
idx = rand_vmax>=edges_vmax(2) & rand_vmax<=edges_vmax(3);
prob_vmax(idx) = 1;
%linear increasing at left edge (from 1 up to 4)
idx = rand_vmax>=edges_vmax(1) & rand_vmax<edges_vmax(2);
prob_vmax(idx) = (rand_vmax(idx)-edges_vmax(1))/(edges_vmax(2)-edges_vmax(1));
%linear decreasing at right edge (from 8 up to 11)
idx = rand_vmax>edges_vmax(3) & rand_vmax<=edges_vmax(end);
prob_vmax(idx) = 1-(rand_vmax(idx)-edges_vmax(3))/(edges_vmax(end)-edges_vmax(3));

%calculate area under curve to normalize probability
a1 = ((edges_vmax(2)-edges_vmax(1))*1)/2;
a2 = ((edges_vmax(end)-edges_vmax(3))*1)/2;
a3 = (edges_vmax(3)-edges_vmax(2))*1;

prob_vmax = prob_vmax/(a1+a2+a3);
clear a1 a2 a3 idx

%%
%calculate probability for random phi

if binORlin == 'bin'
    prob_phi = zeros(size(rand_phi));

    %set values out of range = 0
    idx = rand_phi<edges_phi(1) | rand_phi>=edges_phi(end);
    prob_phi(idx) = 0;

    for i=1:numel(count_phi)
        idx = rand_phi>=edges_phi(i) & rand_phi<edges_phi(i+1);
        prob_phi(idx) = count_phi(i);
    end
    clear idx i
end

if binORlin == 'lin'
    %for linear approach
    prob_phi = zeros(size(rand_phi));
    edges_phi = [edges_phi(1) edges_phi(1:end-1)+width/2 edges_phi(end)];
    count_phi = [0 count_phi 0];

    idx = rand_phi<edges_phi(1) | rand_phi>=edges_phi(end);
    prob_phi(idx) = 0;
    a = 0;
    for i=1:numel(count_phi)-1
        idx = rand_phi>=edges_phi(i) & rand_phi<edges_phi(i+1);
        prob_phi(idx) = count_phi(i)+((count_phi(i+1)-count_phi(i))...
            /(edges_phi(i+1)-edges_phi(i)))*(rand_phi(idx)-edges_phi(i));
        %a = a+abs((count_phi(i+1)-count_phi(i))*(edges_phi(i+1)-edges_phi(i))/2);
    end
    %normalize
    %prob_phi = prob_phi/a;
    clear idx i
end

%sort random numbers --> better for plot later
figure
[~,idx] = sort(rand_vmax);
plot(rand_vmax(idx),prob_vmax(idx));
figure
[~,idx] = sort(rand_phi);
plot(rand_phi(idx),prob_phi(idx));

sRand_phi = sort(rand_phi);

%%
%generate uniform random numbers (rand_wig) and look for each set of parameters if
%probability product (prob_vmax*prob_phi) greater than rand_wig --> if yes
%choose as a parameter set if not remove
rand_wig = rand([1,numrun]);

idx = find((prob_vmax.*prob_phi) > rand_wig);

para_set(:,1) = rand_vmax(idx);
para_set(:,2) = rand_phi(idx);

figure
h = histogram2(para_set(:,1),para_set(:,2),'Normalization','count');
[histx,histy] = meshgrid(h.YBinEdges(1:end-1),h.XBinEdges(1:end-1));
value = h.Values;

%pcolor(value)
figure
contour(histx,histy,value)
figure
contourf(histx,histy,value)
%hold on
%plot(para_set(:,1),para_set(:,2),'.')
%scatter3(para_set(:,1),para_set(:,2),ones(size(para_set(:,1))),[],ones(size(para_set(:,1))))




end

