
function [RG, node_positions, leaf_data]=giveReebGraph(connCompImage,topLeftOffset,resolution)

[r,c]=size(connCompImage);

[y,x]=find(connCompImage==1);

x_bar=round(mean(x));
y_bar=round(mean(y));

if x_bar<1, x_bar=1; end % boundary check
if x_bar>c, x_bar=c; end
if y_bar<1, y_bar=1; end
if y_bar>r, x_bar=r; end

node_row=1; % first node at centroid
node_col=1;
vertex_x=x_bar;
vertex_y=y_bar;

circle_image=zeros(r,c);
circle_image(y_bar,x_bar)=1; % the center point is 1;
circle_image=bwdist(circle_image); % now we have a circle function

radius=resolution; % current radius
dummy_image=zeros(r,c); % will use this for conn comp of level sets
dummy_image(find(circle_image<=radius))=1; % conn comp within the current radius
prev_circle=bwperim(dummy_image); % levet set at current radius
level_set_at_prev_circle=zeros(r,c);
level_set_at_prev_circle(find(prev_circle==1 && connCompImage==1))=1;

% find how many conn comps and process them

[conn_comp_at_prev_circle,num]=bwlabel(level_set_at_prev_circle);

for count=1:num, % record the vertex positions of first level set
    
    [i_cc,j_cc]=find(conn_comp_at_prev_circle==count);
    vertex_x(end+1)=mean(j_cc);
    vertex_y(end+1)=mean(i_cc);
    
end

compute=1;

% now go on adding vertices to the Reeb graph

while compute==1, % graph can be extended
    
   radius=radius+resolution; % we will check next level set
   dummy_image=zeros(r,c); % will use this for conn comp of level sets
   dummy_image(find(circle_image<=radius))=1; % conn comp within the current radius
   curr_circle=bwperim(dummy_image); % levet set at current radius
   level_set_at_curr_circle=zeros(r,c);
   level_set_at_curr_circle(find(curr_circle==1 && connCompImage==1))=1;
   [conn_comp_at_curr_circle,num]=bwlabel(level_set_at_curr_circle);
   
   if num==0, % no more possible nodes
       
       compute=0;
       
   elseif sum(sum(curr_circle-level_set_at_curr_circle))==0, % no new morphology
       
       % do nothing
       
   else % there are more than 1 conn comps
       
       interior_of_curr=zeros(r,c);
       interior_of_curr(find(circle_image<=(radius-resolution)))=1;
       exterior_of_prev=zeros(r,c);
       exterior_of_prev(find(circle_image>=radius))=1;
       intercepted_area=interior_of_curr & exterior_of_prev & connCompImage;
       
       for count=1:num, % create proper vertex connectivity
           
               [i_cc,j_cc]=find(conn_comp_at_curr_circle==count);
               vertex_x(count)=mean(j_cc);
               vertex_y(count)=mean(i_cc);
           
           
           
       end
   end
   
    
end
    
    

