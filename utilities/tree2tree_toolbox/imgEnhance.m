
function eImg = imgEnhance(Img,avgTimes,option)
% Enhance the Image
% USAGE: eImg = imgEnhance(Img,scaleRange,filterRad,avgTimes,option)
% Img: 3D image , scaleRange: a row vector of the scales to look for (for Hessian based)
% avgTimes: How many times we want to do mean filtering 
% option : '1' for hessian based , '2 for AD' , '3' for AD based Hessian
% eImg : the image after enhancement and contrast stretch

if (option == 1) % For Hessian Based Filtering
    fprintf('\n Hessian Based Enhancement \n');
    dlg_title = 'Enter Parameters for Hessian based enhancement';
    prompt = {'Min. Scale','Max. Scale','Scale Interval','Filter Radius'};
    num_lines = 1;
    def = {'1','4','1','8'};
    answer = inputdlg(prompt,dlg_title,num_lines,def,'on');
     
    smin = str2num(answer{1}); smax = str2num(answer{2}); s_skip = str2num(answer{3});
    scaleRange = smin:s_skip:smax;
    filterRad = str2num(answer{4});
    ImStack = [];
    
    wait_cnt = 0;
%     h = waitbar(0,'Enhancement in Progress');
    gama = 1.5;
    parfor I=1:length(scaleRange),
%         wait_cnt = wait_cnt+1;
%         h = waitbar(wait_cnt/length(scaleRange));
        s = scaleRange(I);
%         filterRad = 3*s;
        [x1,x2,x3] = ndgrid(-filterRad:1:filterRad, -filterRad:1:filterRad, -filterRad:1:filterRad);
        gaussFilt = exp((-x1.^2-x2.^2-x3.^2)/(2*s));
        gaussFilt = gaussFilt/sum(gaussFilt(:));
        filtImg = imfilter(Img,gaussFilt,'same');
        ImStack = [ImStack s^gama*vesselness3D(filtImg)];
        disp(strcat('scale = ',num2str(s),' in imEnhance computed '));
    end
    
    ImStack = ImStack';
    if length(scaleRange) > 1
        eImg = max(ImStack);
    else
        eImg = ImStack;
    end
    eImg = reshape(eImg,size(Img));
%     close(h);
    
elseif (option == 2)    % Anisotropic Diffusion
    
    fprintf('\nEnhancing by Anisotropic Diffusion \n');
    dlg_title = 'Enter Parameters for Anisotropic Diffusion based enhancement';
    prompt = {'Delta','Kappa','# Iterations','smoothing function(1/2)'};
    num_lines = 1;
    def = {'3/44','2','6','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def,'on');
    
    delta = str2num(answer{1}); 
    kappa = str2num(answer{2}); 
    num_iter = str2num(answer{3});
    kernel = str2num(answer{4}) ;
    voxel_spacing = ones(3,1);
    eImg =  anisodiff3D(Img,num_iter, delta, kappa, kernel, voxel_spacing);
    
elseif(option == 3)     % Anisotropic Hessian
    
    fprintf('\nAnisotropic Hessian based enhancement\n');
    dlg_title = 'Enter Parameters for Anisotropic Hessian based enhancement';
    prompt = {'Delta','Kappa','# Iterations','smoothing function(1/2)'};
    num_lines = 1;
    def = {'3/44','2','8','1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def,'on');
    
    delta = str2num(answer{1}); 
    kappa = str2num(answer{2}); 
    num_iter = str2num(answer{3});
    kernel = str2num(answer{4}) ;
    voxel_spacing = ones(3,1);
    eImg = anisotropicHessian(Img,num_iter, delta, kappa, kernel, voxel_spacing);
       
end

    avgFilter=ones(3,3,3);
    
    for I=1:avgTimes,
        eImg = imfilter(eImg,avgFilter,'same');
    end

eImg = (eImg-min(eImg(:)))/(max(eImg(:))-min(eImg(:)));


