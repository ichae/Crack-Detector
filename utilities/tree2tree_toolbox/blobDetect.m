function [bImg] = blobDetect(Img)
%Detect the blobs

[row col depth] = size(Img);

fprintf('\n Blob Detection \n');
    dlg_title = 'Enter Parameters for Blob Detection';
    prompt = {'Min. Scale','Max. Scale','Scale Interval'};
    num_lines = 1;
    def = {'4','6','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def,'on');
     
    smin = str2num(answer{1}); smax = str2num(answer{2}); s_skip = str2num(answer{3});
    scaleRange = smin:s_skip:smax;
%     filterRad = str2num(answer{4});
    ImStack = [];
    
    wait_cnt = 0;
    h = waitbar(0,'Blobness detection in Progress');
    gama = 2;
    
    I_w = (fftn(Img));

    parfor I=1:length(scaleRange),

        s = scaleRange(I);
        xRad = round(row/2);
        yRad = round(col/2);
        zRad = round(depth/2);
        
        [x1,x2,x3] = ndgrid(-xRad:1:xRad, -yRad:1:yRad, -zRad:1:zRad);
        gaussFilt = exp((-x1.^2-x2.^2-x3.^2)/(2*s));
        G_w = (fftn(gaussFilt,size(Img)));
        filtImg = abs(ifftn(G_w.*I_w));
        
        filtImg = log(1+filtImg);
        filtImg = (filtImg - min(filtImg(:)))/(max(filtImg(:))-min(filtImg(:)));

        filtImg = ifftshift(filtImg);

        ImStack = [ImStack s^gama*blobness3D(filtImg)];
        disp(strcat('scale = ',num2str(s),' in blobDetect computed '));
    end
    
    ImStack = ImStack';
    bImg = max(ImStack);
    bImg = reshape(bImg,size(Img));
    fprintf('\n Blob Detection completed');

    close(h);

end

