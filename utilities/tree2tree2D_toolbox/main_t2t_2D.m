
clear all;
close all;

clc;

%%local threshold for neurons

[file1, path1]=uigetfile('*.tif', 'Select neuron image to load');
inputNeuronFile=strcat(path1,file1);
Image = im2double(imread(inputNeuronFile));

Image = rgb2gray(Image);
% Parameters

winSz = [15 15];
resolution=5;
erosion_count=1;
percent_neuron_size=10;
k=15; % no of nearest neighbors
alpha=0.4;
beta=0.9;
method='Prim'; % change to Kruskal in case Prim does not work

tic

% Convert to binary
BWcell = neuronBW(Image,1-(0.01*percent_neuron_size),winSz,erosion_count);

'after neuronBW'

toc

figure; imshow(BWcell,[],'init','fit');

%%
% pause;



tic

 % make the k-NN graph with nodes as the connected components in BWcell
 
[neuronPrimitiveGraph,connCompNodes]=make_kNN_graph(BWcell,k,resolution);



% create a minimum spanning tree from the k-NN graph

[MST, pred] = graphminspantree(neuronPrimitiveGraph.adjmat,'method','Kruskal');
plot_MST(MST,connCompNodes,neuronPrimitiveGraph.leaf_connectivity,Image);
% return

%%

[filtMST,newConnCompNodes,new_leaf_connectivity]=alphaBetaFilterGraph(...
                                                          MST,connCompNodes,...
                                                          neuronPrimitiveGraph.leaf_connectivity,...
                                                          alpha,beta);

%%


plot_MST(filtMST,newConnCompNodes,new_leaf_connectivity,Image);

% return;

[splinedMedialGraph,spNodePositions] = fitSplineToMST(MST,...
                     connCompNodes,neuronPrimitiveGraph.leaf_connectivity,2);
                 
figure;imshow(Image,[],'init', 'fit');hold on;
% scatter(spNodePositions.x,spNodePositions.y,'r','filled');
gplot(splinedMedialGraph,[spNodePositions.x;spNodePositions.y]','g-');
hold off;



% cl_dist=0.9;
% maxImg=max(Image(:));
% minImg=min(Image(:));
% imageScaled=(Image-minImg)/(maxImg-minImg+0.01);
% grayImage=imadjust(imageScaled);
% 
% [prunedMedialGraph,prunedNodePositions]=...
%                     removeUnlikelyBranches(splinedMedialGraph,...
%                     spNodePositions,cl_dist,grayImage,BWcell,resolution);
%                 
% figure;imshow(Image,[],'init','fit');hold on;
% gplot(prunedMedialGraph,[prunedNodePositions.x'
% prunedNodePositions.y'],'g')
                

                
                