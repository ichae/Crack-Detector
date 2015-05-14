
function skel_3D = find3DSkeleton(binaryVol)

% Implements the Kalman Palagyi algorithm

modified=1;

while modified>0,
    
    modified=0;
    [binaryVol,newMod]=directionalPruning(binaryVol,'U');
    modified=modified+newMod;
%     modified
    [binaryVol,newMod]=directionalPruning(binaryVol,'D');
    modified=modified+newMod;
%     modified
    [binaryVol,newMod]=directionalPruning(binaryVol,'N');
    modified=modified+newMod;
%     modified
    [binaryVol,newMod]=directionalPruning(binaryVol,'S');
    modified=modified+newMod;
%     modified
    [binaryVol,newMod]=directionalPruning(binaryVol,'E');
    modified=modified+newMod;
%     modified
    [binaryVol,newMod]=directionalPruning(binaryVol,'W');
    modified=modified+newMod;
%     modified
    
end

skel_3D=binaryVol;

end

function [afterPruning,howMany]=directionalPruning(binaryVol,direction)

list=[];
[ptList.y,ptList.x,ptList.z]=ind2sub(size(binaryVol),find(binaryVol)); % CHECK ORDER !!: checked
nPts=length(ptList.x);
% nPts

for I=1:nPts,
    
%     if ~isBorderPt(ptList,binaryVol,direction,I), disp('NOT'); pause; end;
    if isBorderPt(ptList,binaryVol,direction,I),
        
%         'borderpoint!'
        Np=collect_26_nbrs(ptList,I,binaryVol);
%         sum(Np)
        if sum(Np)>1 && isSimple(ptList,binaryVol,I),
            
%             'entered'
            list(end+1)=I;
            
        end
        
    end
    
end

howMany=0;
% list
for J=1:length(list),
    
    Np=collect_26_nbrs(ptList,list(J),binaryVol);
    
    if sum(Np)>1 && isSimple(ptList,binaryVol,list(J)),
        
        binaryVol(ptList.y(list(J)),ptList.x(list(J)),ptList.z(list(J)))=0; % CHECK ORDER !! : checked
        howMany=howMany+1;
        
    end

end

afterPruning=binaryVol;

end

function yes=isBorderPt(ptList,binaryVol,direction,ptIndex)

x=ptList.x(ptIndex);
y=ptList.y(ptIndex);
z=ptList.z(ptIndex);
yes=0;

switch direction
    
    case {'U'}
        
        if binaryVol(y,x,z+1)==0, yes=1; end;
        
    case {'D'}
            
        if binaryVol(y,x,z-1)==0, yes=1; end;
        
    case {'N'}
        
        if binaryVol(y+1,x,z)==0, yes=1; end;
        
    case {'S'}
        
        if binaryVol(y-1,x,z)==0, yes=1; end;
        
    case {'E'}
        
        if binaryVol(y,x+1,z)==0, yes=1; end;
        
    case {'W'}
        
        if binaryVol(y,x-1,z)==0, yes=1; end;
        
    otherwise
        
        disp('wrong direction sign');
        
end

end
    
function nbhd=collect_26_nbrs(ptList,ptIndex,binaryVol)

x=ptList.x(ptIndex);
y=ptList.y(ptIndex);
z=ptList.z(ptIndex);

Np.x=zeros(1,26);
Np.y=zeros(1,26);
Np.z=zeros(1,26);

Np.x(1)=x-1;Np.y(1)=y+1;Np.z(1)=z+1; % WNU
Np.x(2)=x;Np.y(2)=y+1;Np.z(2)=z+1; % PNU
Np.x(3)=x+1;Np.y(3)=y+1;Np.z(3)=z+1; % ENU
Np.x(4)=x-1;Np.y(4)=y;Np.z(4)=z+1; % WPU
Np.x(5)=x;Np.y(5)=y;Np.z(5)=z+1; % PPU
Np.x(6)=x+1;Np.y(6)=y;Np.z(6)=z+1; % EPU
Np.x(7)=x-1;Np.y(7)=y-1;Np.z(7)=z+1; % WSU
Np.x(8)=x;Np.y(8)=y-1;Np.z(8)=z+1; % PSU
Np.x(9)=x+1;Np.y(9)=y-1;Np.z(9)=z+1; % ESU
Np.x(10)=x-1;Np.y(10)=y+1;Np.z(10)=z; % WNP
Np.x(11)=x;Np.y(11)=y+1;Np.z(11)=z; % PNP
Np.x(12)=x+1;Np.y(12)=y+1;Np.z(12)=z; % ENP
Np.x(13)=x-1;Np.y(13)=y;Np.z(13)=z; % WPP
% Np{14}.x=x;Np{14}.y=y;Np{14}.z=z; % PPP
Np.x(14)=x+1;Np.y(14)=y;Np.z(14)=z; % EPP
Np.x(15)=x-1;Np.y(15)=y-1;Np.z(15)=z; % WSP
Np.x(16)=x;Np.y(16)=y-1;Np.z(16)=z; % PSP
Np.x(17)=x+1;Np.y(17)=y-1;Np.z(17)=z; % ESP
Np.x(18)=x-1;Np.y(18)=y+1;Np.z(18)=z-1; % WND
Np.x(19)=x;Np.y(19)=y+1;Np.z(19)=z-1; % PND
Np.x(20)=x+1;Np.y(20)=y+1;Np.z(20)=z-1; % END
Np.x(21)=x-1;Np.y(21)=y;Np.z(21)=z-1; % WPD
Np.x(22)=x;Np.y(22)=y;Np.z(22)=z-1; % PPD
Np.x(23)=x+1;Np.y(23)=y;Np.z(23)=z-1; % EPD
Np.x(24)=x-1;Np.y(24)=y-1;Np.z(24)=z-1; % WSD
Np.x(25)=x;Np.y(25)=y-1;Np.z(25)=z-1; % PSD
Np.x(26)=x+1;Np.y(26)=y-1;Np.z(26)=z-1; % ESD

% size(binaryVol)
% Np.y(:)
% Np.x(:)
% Np.z(:)
t = find(Np.z == 0);
Np.z(t) = 1;
linIndex=sub2ind(size(binaryVol),Np.y(:),Np.x(:),Np.z(:));
nbhd=binaryVol(linIndex);

end

function yes=isSimple(ptList,binaryVol,I)

yes=0;
% I
% 'inside isSimple'
% if cond2Satisfied(ptList,binaryVol,I), disp('hurrah!'); end;
% if cond4Satisfied(ptList,binaryVol,I), disp('atleast 4!'); end;

if cond2Satisfied(ptList,binaryVol,I) && cond4Satisfied(ptList,binaryVol,I),
    
    yes=1;
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% function yes=cond2Satisfied(ptList,binaryVol,ptIndex)
% 
% L=zeros(1,26);
% persistent S26
% S26=cell(1,26);
% S26{1}=[];
% S26{2}=1;
% S26{3}=2;
% S26{4}=[1 2];
% S26{5}=[1 2 3 4];
% S26{6}=[2 3 5];
% S26{7}=[4 5];
% S26{8}=[4 5 6 7];
% S26{9}=[5 6 8];
% S26{10}=[1 2 4 5];
% S26{11}=[1 2 3 4 5 6 10];
% S26{12}=[2 3 5 6 11];
% S26{13}=[1 2 4 5 7 8 10 11];
% S26{14}=[2 3 5 6 8 9 11 12];
% S26{15}=[4 5 7 8 13];
% S26{16}=[4 5 6 7 8 9 13 14 15];
% S26{17}=[4 6 8 9 14 16];
% S26{18}=[10 11 13 17];
% S26{19}=[10 11 12 13 14 18];
% S26{20}=[11 12 14 19];
% S26{21}=[10 11 13 15 16 18 19];
% S26{22}=[10 11 12 13 14 15 16 17 18 19 20 21];
% S26{23}=[11 12 14 16 17 19 20 22];
% S26{24}=[13 15 16 21 22];
% S26{25}=[13 14 15 16 17 21 22 23 24];
% S26{26}=[14 16 17 22 23 25];
% 
% Np=collect_26_nbrs(ptList,ptIndex,binaryVol);
% % sum(Np)
% % pause
% yes=1;
% label=0;
% 
% for I=1:26,
%     
%   label=label+Np(I);
%     
%   for J=1:length(S26{I}),
%        
%        if L(S26{I}(J))>0,
%            
%            for K=1:I,
%                
%                if L(K)==L(J),
%                    
%                    L(K)=label;
%                    
%                end
%                
%            end
%            
%        end
%   end
%   
% end
% 
% % Np
% % L
% % label
% % pause
% 
% for I=1:26,
%     
%     if Np(I)==1 && L(I)~=label,
%         
%         yes=0;
%         
%     end  
%     
% end
% 
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function yes=cond4Satisfied(ptList,binaryVol,ptIndex)

yes=1;

x=ptList.x(ptIndex);
y=ptList.y(ptIndex);
z=ptList.z(ptIndex);

node_list.x=[];
node_list.y=[];
node_list.z=[];

dummy_x=x;dummy_y=y;dummy_z=z+1; % U
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
%     'U'
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x;dummy_y=y;dummy_z=z-1; % D
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
%     'D'
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x;dummy_y=y+1;dummy_z=z; % N
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
%     'N'
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x;dummy_y=y-1;dummy_z=z; % S
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
%     'S'
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x+1;dummy_y=y;dummy_z=z; % E
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
%     'E'
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x-1;dummy_y=y;dummy_z=z; % W
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
%     'W'
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

n_crit_nodes=length(node_list.x);

dummy_x=x+1;dummy_y=y+1;dummy_z=z; % NE
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x+1;dummy_y=y-1;dummy_z=z; % SE
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x-1;dummy_y=y+1;dummy_z=z; % NW
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x-1;dummy_y=y-1;dummy_z=z; % SW
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x;dummy_y=y+1;dummy_z=z+1; % UN
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x;dummy_y=y-1;dummy_z=z+1; % US
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x;dummy_y=y+1;dummy_z=z-1; % DN
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x;dummy_y=y-1;dummy_z=z-1; % DS
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x+1;dummy_y=y;dummy_z=z+1; % UE
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x+1;dummy_y=y;dummy_z=z-1; % DE
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x-1;dummy_y=y;dummy_z=z+1; % UW
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x-1;dummy_y=y;dummy_z=z-1; % DW
if binaryVol(dummy_y,dummy_x,dummy_z)==0,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

n_total_nodes=length(node_list.x);

% n_crit_nodes
% n_total_nodes
% node_list
% pause

conn_graph=zeros(n_total_nodes);

for I=1:n_total_nodes,
    
    for J=1:I-1,
        
        dist=(node_list.x(I)-node_list.x(J)).^2+ ...
             (node_list.y(I)-node_list.y(J)).^2+ ...
             (node_list.z(I)-node_list.z(J)).^2;
                    
        if dist<=1, conn_graph(I,J)=1; end; % Check 6-adjacency
                    
    end
    
end

conn_graph=conn_graph+conn_graph'+diag(ones(1,n_total_nodes));

% sparse(conn_graph)
% pause

distMat=graphallshortestpaths(sparse(conn_graph),'directed',false);

conn_6_graph=distMat(1:n_crit_nodes,1:n_crit_nodes);
% conn_6_graph
% isinf(conn_6_graph)
% sum(sum(isinf(conn_6_graph)))
% pause

if sum(sum(isinf(conn_6_graph)))>0, yes=0; end; % not 6-connected

end

function yes=cond2Satisfied(ptList,binaryVol,ptIndex)

yes=1;

x=ptList.x(ptIndex);
y=ptList.y(ptIndex);
z=ptList.z(ptIndex);

node_list.x=[];
node_list.y=[];
node_list.z=[];

dummy_x=x;dummy_y=y;dummy_z=z+1; % U
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
%     'U'
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x;dummy_y=y;dummy_z=z-1; % D
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
%     'D'
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x;dummy_y=y+1;dummy_z=z; % N
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
%     'N'
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x;dummy_y=y-1;dummy_z=z; % S
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
%     'S'
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x+1;dummy_y=y;dummy_z=z; % E
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
%     'E'
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x-1;dummy_y=y;dummy_z=z; % W
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
%     'W'
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

% n_crit_nodes=length(node_list.x);

dummy_x=x+1;dummy_y=y+1;dummy_z=z; % NE
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x+1;dummy_y=y-1;dummy_z=z; % SE
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x-1;dummy_y=y+1;dummy_z=z; % NW
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x-1;dummy_y=y-1;dummy_z=z; % SW
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x;dummy_y=y+1;dummy_z=z+1; % UN
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x;dummy_y=y-1;dummy_z=z+1; % US
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x;dummy_y=y+1;dummy_z=z-1; % DN
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x;dummy_y=y-1;dummy_z=z-1; % DS
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x+1;dummy_y=y;dummy_z=z+1; % UE
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x+1;dummy_y=y;dummy_z=z-1; % DE
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x-1;dummy_y=y;dummy_z=z+1; % UW
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x-1;dummy_y=y;dummy_z=z-1; % DW
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x+1;dummy_y=y+1;dummy_z=z+1; % UNE
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x+1;dummy_y=y+1;dummy_z=z-1; % DNE
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x-1;dummy_y=y+1;dummy_z=z+1; % UNW
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x-1;dummy_y=y+1;dummy_z=z-1; % DNW
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x-1;dummy_y=y-1;dummy_z=z+1; % USW
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x+1;dummy_y=y-1;dummy_z=z+1; % USE
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x-1;dummy_y=y-1;dummy_z=z-1; % DSW
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

dummy_x=x+1;dummy_y=y-1;dummy_z=z-1; % DSE
if binaryVol(dummy_y,dummy_x,dummy_z)==1,
   
    node_list.x(end+1)=dummy_x;
    node_list.y(end+1)=dummy_y;
    node_list.z(end+1)=dummy_z;
    
end

n_total_nodes=length(node_list.x);

conn_graph=zeros(n_total_nodes);

for I=1:n_total_nodes,
    
    for J=1:I-1,
        
        dist=(node_list.x(I)-node_list.x(J)).^2+ ...
             (node_list.y(I)-node_list.y(J)).^2+ ...
             (node_list.z(I)-node_list.z(J)).^2;
                    
        if dist<=3, conn_graph(I,J)=1; end; % Check 26-adjacency
                    
    end
    
end

conn_graph=conn_graph+conn_graph'+diag(ones(1,n_total_nodes));

distMat=graphallshortestpaths(sparse(conn_graph),'directed',false);

if sum(sum(isinf(distMat)))>0, yes=0; end; % not 26-connected

end

