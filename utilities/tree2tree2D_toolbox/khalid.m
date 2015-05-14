
% clear all;
% close all;
% clc;

[file1, path1]=uigetfile('*.tif', 'Select neuron image to load');
inputNeuronFile=strcat(path1,file1);
Image = double(imread(inputNeuronFile));
expand=60;
[r,c]=size(Image);
imbig=zeros(r+(2*expand),c+(2*expand));
imbig(expand+1:expand+r,expand+1:expand+c)=Image;

figure; imshow(Image,[],'init','fit'); hold on;

sp=ginput(1);
spx=sp(1);
spy=sp(2);

scatter(spx,spy,20,'y+');

LPD=[-1 -2 0 2 1]';
K=5;
step=3;
filterL=repmat(LPD,1,K);
filterR=repmat(flipud(LPD),1,K);

N=32;
angle=linspace(0,180,N+1);

width=1;
middle=zeros(width*2+1,K);

filter=[filterL;middle;filterR];
cx=spx;
cy=spy;
x=[];
y=[];

for I=1:100,
  
%   I
  x(end+1)=cx;
  y(end+1)=cy;
  
  buffer=imbig(round(cy)-2*(width+K)+expand:round(cy)+2*(width+K)+expand,...
               round(cx)-2*(width+K)+expand:round(cx)+2*(width+K)+expand);
  best_val=-99999999;
  best_dir=0;
  
  for J=1:N,
      
%       J
      buffer_r=imrotate(buffer,angle(J));
      [r1,c1]=size(buffer_r);
      midy=round(r1/2);
      midx=round(c1/2);
      to_match=buffer_r(midy-K-width:midy+K+width,...
                        midx-2:midx+2);
      matchval=sum(sum(to_match.*filter));
%       
%       sum(buffer_r(:))
%       sum(to_match(:))
%       matchval
%       
      if matchval>best_val,
          
          best_val=matchval;
          best_dir=J;
          
      end
      
  end
%   
%   best_val
%   best_dir
%   
  prevx=cx;
  prevy=cy;
  
  cx=prevx+step*cos(angle(best_dir));
  cy=prevy+step*sin(angle(best_dir));
    
  plot([prevx,cx],[prevy,cy],'y','Linewidth',3);drawnow;
  
%   pause;
  
end


      
      
      