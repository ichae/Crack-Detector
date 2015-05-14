% This program automatically segments neurons from the 3D stack
% The image stacks are loaded by the read3DStack routine
% Author: Suvadip Mukherjee (sm5vp@virginia.edu)
% Acknowledgement: The original code was written by Saurav Basu

%%
clear all;
close all;
clc;

% Read the 3D Stack
skip_interval = 3 ;       % skip intermediate stacks for memory issues
downscale = 0.5 ;         % downsize the image for memory issues
Img = read3DStack(downscale,skip_interval);
fprintf('Neuron Stack is read');

origImg = Img(:,:,2:size(Img,3)-1);
% Write the 3D stack. This is the downsampled and downscaled version
write3DStack(origImg,1,1,'original');

%------------ Code to Run --------------------------------------------
suvadip = 1; % Run Suvadip's code
if suvadip == 1
    threshOption = 2;       % 1: otsu, 2: local adaptive
    pruning_option = 2 ;    % 1: Alpha Beta , 2: Path Cost Based
    fit_option = 2;         % 1 for old method,2 for PathSearch
    Rmin = 1;
    Rmax = 2.5;
else
    threshOption = 1;       % 1: otsu, 2: local adaptive
    pruning_option = 1 ;    % 1: Alpha Beta , 2: Path Cost Based
    fit_option = 1;         % 1 for old method,2 for PathSearch
    Rmin = 2;
    Rmax = 2;
end

% -------------------------------------------------------------------------
% -------------     ENHANCEMENT STAGE      --------------------------------
% -------------------------------------------------------------------------
%%

option = 1 ;    % option = 1 for hessian based, option = 2 for AD, option = 3 for Anisotropic Hessian
avgTimes = 1;   % times to perform the averaging
eImg = imgEnhance(Img,avgTimes,option);
fprintf('Enhancement step is completed\n');
eImg = eImg(:,:,2:size(eImg,3)-1);

write3DStack(eImg,1,1,'enhanced');


%%
% -------------------------------------------------------------------------
% -------------     BINARIZATION STAGE      -------------------------------
% -------------------------------------------------------------------------

threshOption = 2;       % 1: otsu, 2: local adaptive
tholdLevelReduce = .3;  % reduce the threshold limit obtained from Otsu
openSz = 3;             % size of the structuring element to open
minObjSize = 0;         % objectsizes below this will be removed
bndCut = 5;

segI = binarize3D(Img,eImg,tholdLevelReduce,openSz,minObjSize,bndCut,threshOption);
fprintf('\nBinarization Completed\n');
write3DStack(segI,1,1,'threshold'); % just to check

%% Remove small components/dusts
tmp = segI;
areaOpenSz = 80;
segI = bwareaopen(tmp,areaOpenSz,26);
write3DStack(segI,1,1,'threshold');

%% Get the skeleton

% skeleton3D = find3DSkeleton(segI);
% write3DStack(skeleton3D,1,1,'skeleton');

% Get the Markerfile...can be dropped on V3D
% toMarkerfile(skeleton3D);


%% Connect the components
% -------------------------------------------------------------------------
% -------------     Connect the Components  -------------------------------
% -------------------------------------------------------------------------

connComps = get3DConnComps(segI);
num_cc = size(connComps,2);
clear connComps ;

resolution = 2;
k = 10;
[neuronPrimitiveGraph,connCompNodes] = make_3D_kNN_graph(segI,k,resolution);

%% Minimum Spanning Tree -- select the best of the K-MST
% -------------------------------------------------------------------------
% -------------  Create Greedy K-MST   ------------------------------------
% -------------------------------------------------------------------------

clc
xy_offset = 10;
z_offset = xy_offset/5;
idx = 1;
delete_num_nodes = 6;

subGraph = [];
subGraph{idx}.graphMat = neuronPrimitiveGraph.adjmat;
subGraph{idx}.nodeList = connCompNodes;  % contains info about each node -- medial graph/nodePosition/LeafData/NodeValue
method = 'Prim';                         % Change to 'Kruskal' if Prim doesnt work
MST = graphminspantree(neuronPrimitiveGraph.adjmat,'method',method);            % MST after deletion of a node
subGraph{idx}.MST = MST;
[pathcost, pathBetweenLeaves] = pathSearch(Img, MST, connCompNodes , neuronPrimitiveGraph.leaf_connectivity,eImg,xy_offset,z_offset);
subGraph{idx}.pathCost = pathcost;
subGraph{idx}.leafConnectivity = neuronPrimitiveGraph.leaf_connectivity;
subGraph{idx}.pathBetweenLeaves = pathBetweenLeaves;


%%
idx=1;
while idx < delete_num_nodes && idx < num_cc - 1
    
    idx = idx+1;
       
    node2delete = findNode2delete(subGraph{idx-1}); % decide which node to delete
    
    [modifiedGraph,modifiedNodeList,modifiedleafConnectivity,~] = modifyGraph(node2delete, subGraph{idx-1}); 
    
    subGraph{idx}.graphMat = modifiedGraph;
    subGraph{idx}.nodeList = modifiedNodeList;  % contains info about each node -- medial graph/nodePosition/LeafData/NodeValue
    subGraph{idx}.leafConnectivity = modifiedleafConnectivity;
        
    [modifiedMST, pred] = graphminspantree(subGraph{idx}.graphMat,'method',method);
    [modifiedPathcost, modifiedPathBetweenLeaves] = pathSearch(Img, modifiedMST, modifiedNodeList , modifiedleafConnectivity,eImg,xy_offset,z_offset);

    subGraph{idx}.MST = modifiedMST;            % MST after deletion of a node
    subGraph{idx}.pathCost = modifiedPathcost;
    subGraph{idx}.pathBetweenLeaves = modifiedPathBetweenLeaves;
    
end

%% Obtain the optimal grap by solving the optimizer
%   SubGraph -- contains all the k-subMST's
clc
tau1 = 1;
tau2 = 1;

[optimalTree,score] = findOptimalTree(subGraph,tau1,tau2);
stem(score);
filtMST = optimalTree{1}.MST;
newConnCompNodes = optimalTree{1}.nodeList;
new_leaf_connectivity = optimalTree{1}.leafConnectivity;
pathBetweenCC = optimalTree{1}.pathBetweenLeaves;

%% Global Tree
% -------------------------------------------------------------------------
% -----------------------   Global Tree   ---------------------------------
% -------------------------------------------------------------------------

% check if filtMST is all zeros: SUVADIP-- bugfix
num_cc = length(optimalTree{1}.nodeList);
if num_cc == 1
    filtMST = [];
    newConnCompNodes = connCompNodes;
    new_leaf_connectivity = [];
    pathBetweenCC = [];
end

[splinedMedialGraph,spNodePositions] = PathFollow(filtMST,newConnCompNodes,new_leaf_connectivity,pathBetweenCC,0.2);

%% create the .swc file
% -------------------------------------------------------------------------
% -----------------------   Create SWC    ---------------------------------
% -------------------------------------------------------------------------
clc
[row col depth] = size(eImg);
flip = 1;
Rmin = 1;
Rmax = 2;

toSWC(origImg,splinedMedialGraph,spNodePositions,1,row,col,depth,Rmin,Rmax,flip);
delete('*.asv');
