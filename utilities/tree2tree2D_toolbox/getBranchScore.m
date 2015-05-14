
function score=getBranchScore(branch_buffer,branch_x,branch_y,...
               grayImage,bwImage,nbhd_image,res)

branch_length=length(branch_buffer);
[r,c]=size(grayImage);
score=0;

% branch_buffer
% pause

% first handle the case when branch length > 2
% make an ellipse with minor axis along branch tanget
% and major axis perpendicular to that
% with major axis = 2* minor axis

L=2;
% t=linspace(0,2*pi,51);
% ellipse_x1=2*L*res*cos(t(1:end-1));
% ellipse_y1=L*res*sin(t(1:end-1));
box_x0=[-(L/2)*res (L/2)*res (L/2)*res -(L/2)*res];
box_y0=[L*res L*res -L*res -L*res];

% figure;imshow(nbhd_image,[]);hold on;pause;

for I=2:branch_length-1,
    
    if bwImage(branch_y(I),branch_x(I))~=1, % was 0 from local threshold
        
        minor_axis_x=branch_x(I+1)-branch_x(I);
        minor_axis_y=branch_y(I+1)-branch_y(I);
        angle=givePrincipalArgumentAngle(minor_axis_x,minor_axis_y);
        
        % put the major axis of ellipse
        % perpendicular to the branch direction
        offsetAngle=angle+pi/2;
        
%         offsetAngle
        
        center_x=branch_x(I);
        center_y=branch_y(I);
        
        % rotate by offsetAngle
        
        box=[cos(offsetAngle) -sin(offsetAngle);...
                 sin(offsetAngle) cos(offsetAngle)]*...
                 [box_x0;box_y0];
             
        box_x=box(1,:)+center_x;
        box_y=box(2,:)+center_y;
        
        % boundary check
        box_x=max(1,box_x);
        box_x=min(box_x,c);
        box_y=max(1,box_y);
        box_y=min(box_y,r);
        
        bw_box=roipoly(grayImage,box_x,box_y);
        grayBox=bw_box.*grayImage; % take a window centered on the neuron
        tubular_region=bw_box.*nbhd_image;
        bkgrnd_image=(~tubular_region).*bw_box;
        grayBox0=imrotate(grayBox,-offsetAngle);
        bw_box0=imrotate(bw_box,-offsetAngle);
        ii=find(bw_box0==1);
        gray_cropped=grayBox0(ii);
        B=floor(length(ii)/L);
        gray_cropped=gray_cropped(1:L*B); % making sure that there are exactly L*B elements
        gray_cropped=reshape(gray_cropped,B,L);
        level=graythresh(gray_cropped);
        fgrnd=im2bw(grayBox,level).*bw_box;
        tubular_mean=sum(sum(bw_box.*fgrnd.*nbhd_image))...
                       /sum(sum(bw_box.*nbhd_image));
        outside_mean=sum(sum(bw_box.*fgrnd.*bkgrnd_image))...
                       /sum(sum(bw_box.*bkgrnd_image));
  
        score=score+abs(outside_mean-tubular_mean);
        
%         score  
%         ellipse'
%         imshow((bw_ellipse | nbhd_image),[]);pause;
%         imshow((bw_ellipse & nbhd_image),[]);pause;
        
    
    else % if already detected as neuron
        
        score=score+1; % give full weight
        
    end
    
end

if branch_length==2,
    
    score=1;
    
end

score=1-(score/(branch_length+0.1));
