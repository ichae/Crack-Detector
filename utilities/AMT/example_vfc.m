% Vector field convolution (VFC) external force field example.
% 
%     See also AMT, EXAMPLE_PIG, AM_VFC, AM_VFK, AC_DISPLAY.
% 
%     Reference
%     [1] Bing Li and Scott T. Acton, "Active contour external force using
%     vector field convolution for image segmentation," Image Processing,
%     IEEE Trans. on, vol. 16, pp. 2096-2106, 2007.  
%     [2] Bing Li and Scott T. Acton, "Automatic Active Model
%     Initialization via Poisson Inverse Gradient," Image Processing,
%     IEEE Trans. on, vol. 17, pp. 1406-1420, 2008.   
% 
% (c) Copyright Bing Li 2005 - 2009.

clear all
disp('======================================')
disp('Vector field convolution (VFC) example')

%% parameter settings
disp('Initializing parameters ...')
SAVE_AVI = 0;           % set it to 1 if you want to save the process as .avi movie
DISPLAY_STREAMLINE = 0; % set it to 1 if you want to plot streamlines, note that it takes a while
mu = .2;
GVF_ITER = 100;
normalize = 1;
alpha = .5;
beta = 0;
tau = .5;
SNAKE_ITER = 5;
SNAKE_ITER1 = 60;
RES = .5;
clr = {'b' 'b' 'r'};

%% Read images
disp('Reading images ...')
U = imread('im_U.bmp');
noisyU=imread('im_Unoisy.bmp');
figure(1)

%% compare 3 different cases
for cs = 1:3,
    %% compute external force fields
    switch cs,
        case 1, % traditional GVF with Gaussian filter
            disp('--------------------------------------------------')
            disp('Case 1: GVF snake with initial circle close to FOI')
            disp('Computing the external force field ...')
            h = fspecial('gaussian',[5 5],5);
            f = imfilter(double(noisyU),h);
            titl = 'GVF';
            Fext = AM_GVF(f, mu, GVF_ITER, normalize);
            R = 20;
        case 2, % traditional GVF with Gaussian filter
            disp('--------------------------------------------------')
            disp('Case 2: GVF snake with initial circle far away from FOI')
            disp('Computing the external force field ...')
            h = fspecial('gaussian',[5 5],5);
            f = imfilter(double(noisyU),h);
            titl = 'GVF';
            Fext = AM_GVF(f, mu, GVF_ITER, normalize);
            R = 28;
        case 3, % VFC
            disp('--------------------------------------------------')
            disp('Case 3: VFC snake with initial circle far away from FOI')
            disp('Computing the external force field ...')
            f = noisyU;
            K = AM_VFK(2, 32, 'power',1.8);
            Fext = AM_VFC(f, K, 1);
            R = 28;
            titl = 'VFC';
    end
  
    %% display
    I = (1-noisyU)*.3+.7;   % for display
    subplot(2,3,cs)
    disp('Displaying the external force field ...')

    if DISPLAY_STREAMLINE,   
        cla
        [x y] = meshgrid(.5:64,.5:64);
        vt = [x(:) y(:)];   % seeds
        VT = zeros([size(vt) 40]);
        VT(:,:,1) = vt;
        for i=1:39, % moving these seeds
            vt = AC_deform(vt,0,0,tau,Fext,1);
            VT(:,:,i+1) = vt;
        end

        [Ty Tx] = find(~U);  % ground truth  
        hold on
        for i=1:size(vt,1),
            if min(abs(VT(i,1,end)-Tx)+abs(VT(i,2,end)-Ty))<=2,  
                % converge to U-shape
                plot(squeeze(VT(i,1,:)), squeeze(VT(i,2,:)),'r','linewidth',1)
            else
                plot(squeeze(VT(i,1,:)), squeeze(VT(i,2,:)),'k','linewidth',1)
            end
        end
        hold off
        axis equal; axis 'ij';
        axis([1 64 1 64])
    else
        AC_quiver(Fext, I);
        title(['normalized ',titl,' field']);
    end
    %% uncomment these 2 lines to save the display 
    %     F = getframe(gca);
    %     imwrite(F.cdata,['vector_field',num2str(cs),'.bmp']);

    %% initialize a circle at (32 32) with radius R
    disp('Initializing the snake ...')
    vert  = AC_initial(RES, 'circle', [32 32 R]);    
    vert0 = vert;

    subplot(2,3,3+cs)
    imshow(I)
    AC_display(vert,'close',clr{cs});
    drawnow, pause(.5)

    if SAVE_AVI
        mov = avifile(['example_vfc_',num2str(cs),'.avi'],'fps',4,'quality',100,'compression','None');
        frame = getframe(gca);
        mov = addframe(mov,frame);
    end

    disp('Deforming the snake ...')
    for i=1:SNAKE_ITER1,
        vert = AC_deform(vert,alpha,beta,tau,Fext,SNAKE_ITER);
        vert = AC_remesh(vert,.5);

        if mod(i,2)==0,
            imshow(I)
            h=AC_display(vert,'close',clr{cs});
            h=AC_display(vert0,'close',[clr{cs} '--']);
            title([titl ' iteration ' num2str(i)])
            drawnow, pause(.5)
        end

        if SAVE_AVI
            frame = getframe(gca);
            mov = addframe(mov,frame);
        end
    end
    disp('Done!')

    if SAVE_AVI
        h=close(mov);
    end
end