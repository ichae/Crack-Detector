
function write3DStack(Img,sliceReduce,stackSkip,newDir)

% write some portion in case hardware is in adequate

% Img=(Img-min(Img(:)))/(max(Img(:))-min(Img(:)));

% Prepare for writing
message = ['Write The Stacks (' newDir ')'];
dname = uigetdir('.\temporary',message);
[dirmade,msg,msgID] = mkdir(dname,newDir);
p_write = strcat(dname,'\',newDir,'\');

if ~isempty(msg), % directory already exists
    
    delete(strcat(p_write,'*.tif')); % remove the exixting tiffs

end

fname = 'slice';

for f = 1:stackSkip:size(Img,3)
    
   buffer = imresize(Img(:,:,f),sliceReduce);
   imwrite(uint8(buffer.*255),[p_write fname num2str(f) '.tif']);
   
end
