% ACTIVE MODEL TOOLBOX (AMT)
% 
% Acronyms
%   AC  = Active Contour 2D
%   AM  = Active Models
%   CoD = Center of Divergence
%   FFS = Force Field Segmentation
%   GVF = Gradient Vector Flow
%   PIG = Poisson Inverse Gradient 
%   VFC = Vector Field Convolution
%   VFK = Vector Field Kernel
%
% Examples
%   example_vfc  - Vector field convolution (VFC) external force field example.
%   example_pig  - Poisson inverse gradient (PIG) automatic initialization example.
% 
% Functions
%   AC_deform    - Deform an active contour (AC), also known as snake.
%   AC_display   - Display the contour on the current axis.
%   AC_initial   - Initialize an active contour (AC), also known as snake.
%   AC_isoLine   - Select the optimal close isoline(s).
%   AC_quiver    - Display the 2D external force field in arrows.
%   AC_remesh    - Remesh an active contour (AC) to ensure the resolution.
%   AM_CoD       - Initialize active contours using the CoD algorithm [3].
%   AM_FFS       - Initialize active contours using the FFS algorithm [4].
%   AM_gradient  - 2D/3D gradient, a faster implementation than build-in GRADIENT().
%   AM_GVF       - Compute the gradient vector flow (GVF) force field [5].
%   AM_laplacian - 2D/3D discrete Laplacian, a faster implementation than build-in DEL2().
%   AM_PIG       - Estimate the energy via Poisson inverse gradient (PIG) [2].
%   AM_VFC       - Compute the vector field convolution (VFC) force field [1].
%   AM_VFK       - Compute the vector field kernel (VFK).
% 
% References
%   [1] Bing Li and Scott T. Acton, "Active contour external force using
%   vector field convolution for image segmentation," Image Processing,
%   IEEE Trans. on, vol. 16, pp. 2096-2106, 2007.    
%   [2] Bing Li and Scott T. Acton, "Automatic Active Model Initialization
%   via Poisson Inverse Gradient," Image Processing, IEEE Trans. on, vol.
%   17, pp. 1406-1420, 2008.     
%   [3] Xingfei Ge and Jie Tian, "An automatic active contour model for
%   multiple objects," in Pattern Recognition, Proceedings of the
%   International Conference on, 2002.      
%   [4] Chunming Li, Jundong Liu, and Martin D. Fox, "Segmentation of
%   external force field for automatic initialization and splitting of
%   snakes," Pattern Recognition, vol.38, pp. 1947?960, 2005.        
%   [5] Chenyang Xu and Jerry L. Prince, "Snakes, Shapes, and Gradient
%   Vector Flow," Image Processing, IEEE Trans. on, vol. 7, pp. 359-369, 1998.  
% 
% (c) Copyright Bing Li 2005 - 2009.
