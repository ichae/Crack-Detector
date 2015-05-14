classdef StdIP
    %STDIP Class of standard image processing functions
    %   Includes various general purpose IP routine
    %   Author: Suvadip Mukherjee (sm5vp@virginia.edu)
    
    properties
    end
    
    methods(Static)

        function [rgbImg, grayImg, fname, pname] = readImg2D(res)
            [fname,pname] = uigetfile('*.*','Input Data');
            rgbImg = imread(strcat(pname,fname));
            dim = ndims(rgbImg);
            rgbImg = imresize(rgbImg,res);
            if(dim == 3)
                grayImg = rgb2gray(rgbImg);
            else
                grayImg = rgbImg;
            end
            grayImg = im2double(grayImg);
            
        end
        
        function h = imageOverlay(im1,im2,mag,col)
           %-- Overlay im2 over im1. Generally, im2 is binary 
           figure; imshow(im1,'InitialMagnification',mag); hold on;
           R = col(1)*im2; G = col(2)*im2; B = col(3)*im2;
           rgbim2 = cat(3,R,G,B);
           h = imshow(rgbim2,'InitialMagnification',mag); hold off;
           alpha_data = im2;
           set(h,'AlphaData',alpha_data);
        end
       
    end
    
end

