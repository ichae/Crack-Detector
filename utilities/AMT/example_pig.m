% Poisson inverse gradient (PIG) automatic initialization example.
% 
%     See also AMT, EXAMPLE_VFC, AM_PIG, AC_ISOLINE, AM_VFC.
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
disp('Poisson inverse gradient (PIG) example')

%% parameter settings
disp('Initializing parameters ...')

SAVE_AVI = 0;   % set it to 1 if you want to see the deformation and save it as .avi movie
alpha = .5;
beta = .1;
tau = .5;
clr = {'r','m','b'};

%% Read images
disp('Reading images ...')
I0 = imread('im_lung.png');

%% Filter noise
disp('Filtering images ...')
I = 255-imfill(255-I0,'holes');
I = imopen(I, strel('disk',4));

%% Compute edge map and external force field
disp('Computing edge map and force field...')
f = edge(I,'canny',.2,1.5);
omega = ~f;

K = AM_VFK(2, 64, 'power',3);
Fext = AM_VFC(f, K, 1);

%% compare 3 different methods
for method = 1:3
    %% Initialize
    switch method
        case 1 % CoG
            disp('--------------------------------------------------')
            disp('Initializer 1: center of divergence (CoD)')
            disp('Initializing ...')
            titl = 'CoD';
            tic;
            vert = AM_CoD(Fext,6);
            TimeInitial(method) = toc;
        case 2 % FFS
            disp('--------------------------------------------------')
            disp('Initializer 2: force field segmentation (FFS)')
            disp('Initializing ...')
            titl = 'FFS';
            tic
            vert = AM_FFS(Fext,imdilate(f,strel('disk',1)));
            TimeInitial(method) = toc;
        case 3 % PIG
            disp('--------------------------------------------------')
            disp('Initializer 3: Poisson inverse gradient (PIG)')
            disp('Initializing ...')
            titl = 'PIG';
            lambda = 5:6:50;    % isovalues
            eta = 2;            % # of models to initialize

            tic
            E = AM_PIG(Fext, omega, double(f));
            [vert isoValue] = AC_isoLine(E, lambda, eta);  
            TimeInitial(method) = toc;

            figure(3) % display E with isolines
            imagesc(E-max(E(:))+min(E(:))-1);
            hold on
            contour(E,lambda,'linewidth',2.5)
            hold off
            
            title(titl), axis image, axis off
            set(gca,'Clim',[min(E(:))*2-max(E(:))-1,max(E(:))])
            colormap([gray(128); jet(32) ;jet(96)]);
    end
    
    %% display
    figure(2)
    K = length(vert);
    ModelNum(method) = K;
    if SAVE_AVI
        c = clr{method};
    else
        c = 'g--';
    end
    subplot(1,3,method)
    imshow(I0)
    title(titl)
    for k = 1:K
        h = AC_display(vert{k},'close',c);
        set(h,'LineWidth',3)
    end
    drawnow, pause(.5)

    if SAVE_AVI
        title([titl ' iteration 0'])
        mov = avifile(['example_pig_',num2str(method),'.avi'],'fps',4,'quality',100,'compression','None');
        frame = getframe(gca);
        mov = addframe(mov,frame);
    end
    
    %% snake deformation
    disp('Deforming the snake ...')
    if ~SAVE_AVI
        disp('iteration   0')
    end
    area_previous = zeros(1,K);
    for k = 1:K,
        area_previous(k) = polyarea(vert{k}(:,1),vert{k}(:,2));
    end
    area_diff = area_previous;  % compute areas for converge condition
    iter = 0;
    flagsConverged = zeros(1,K);
    
    tic
    while ~all(flagsConverged)
        iter = iter+3;
        flagsConverged = abs(area_diff)<1;  % converge condition: area change less than 1
        for k = 1:K,
            if flagsConverged(k)
                continue
            end
            vert{k} = AC_deform(vert{k},alpha,beta,tau,Fext,3);
            area = polyarea(vert{k}(:,1),vert{k}(:,2));
            area_diff(k) = area - area_previous(k);
            area_previous(k) = area;
            
            if mod(iter,9)==0,
                vert{k} = AC_remesh(vert{k},1);
            end
        end
        if SAVE_AVI && mod(iter,9)==0,
            imshow(I0)
            title([titl ' iteration ' num2str(iter)])
            clear h
            for k = 1:K,
                h(k) = AC_display(vert{k},'close',c);
            end
            set(h(~isnan(h)),'linewidth',3)
            drawnow
            frame = getframe(gca);
            mov = addframe(mov,frame);
        else
            disp(sprintf('\b\b\b\b%3d', iter));
        end
    end
    TimeDeform(method) = toc;
    
    if SAVE_AVI
        mov = close(mov);
    else
        for k = 1:K
            try
                h = AC_display(vert{k},'close','y');
                set(h,'linewidth',3)
            end
        end
        drawnow, pause(.5)
    end
    
    %% report
    disp([num2str(K) ' models are initialized for ' num2str(TimeInitial(method)) ...
        ' seconds, and converge after ' num2str(iter) ' iterations and ' ...
        num2str(TimeDeform(method)) ' seconds.' ])
end