%% BE

wr = -5:0.01:3;
wi = -4:0.01:4;
[Wr,Wi] = meshgrid(wr,wi);
BE_sc = abs(1./(1-(Wr+1i*Wi)));
BE_sc(BE_sc<1) = 1; % stable
BE_sc(BE_sc>1) = 2;
cmap = [ 1 0.5 0 ; 1 1 1 ];
contourf(Wr,Wi,BE_sc,[1 2])
colormap(cmap), axis equal
xlabel( '$\mathcal{R}(\Delta t \lambda_i$)', 'interpreter', 'latex', 'fontsize', 16)
ylabel( '$\mathcal{I}(\Delta t \lambda_i$) ', 'interpreter', 'latex', 'fontsize', 16)
title( 'Backward Euler Region of Absolute Stability', 'interpreter', 'latex', 'fontsize', 16)
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 14 )

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [15 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15])
set(gcf, 'PaperPosition', [0 0 15 15])
svnm = 'BE_stab';
print( '-dpng', svnm, '-r200' );

%% FE

wr = -5:0.01:3;
wi = -4:0.01:4;
[Wr,Wi] = meshgrid(wr,wi);
FE_sc = abs(1-(Wr+1i*Wi));
FE_sc(FE_sc<1) = 1; % stable
FE_sc(FE_sc>1) = 2;
cmap = [ 1 0.5 0 ; 1 1 1 ];
contourf(Wr,Wi,FE_sc,[1 2])
colormap(cmap), axis equal
xlabel( '$\mathcal{R}(\Delta t \lambda_i$)', 'interpreter', 'latex', 'fontsize', 16)
ylabel( '$\mathcal{I}(\Delta t \lambda_i$) ', 'interpreter', 'latex', 'fontsize', 16)
title( 'Forward Euler Region of Absolute Stability', 'interpreter', 'latex', 'fontsize', 16)
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 14 )

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [15 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15])
set(gcf, 'PaperPosition', [0 0 15 15])
svnm = 'FE_stab';
print( '-dpng', svnm, '-r200' );

%% CN

wr = -5:0.01:3;
wi = -4:0.01:4;
[Wr,Wi] = meshgrid(wr,wi);
CN_sc = abs((1+0.5*(Wr+1i*Wi))./(1-0.5*(Wr+1i*Wi)));
CN_sc(CN_sc<1) = 1; % stable
CN_sc(CN_sc>1) = 2;
cmap = [ 1 0.5 0 ; 1 1 1 ];
contourf(Wr,Wi,CN_sc,[1 2])
colormap(cmap), axis equal
xlabel( '$\mathcal{R}(\Delta t \lambda_i$)', 'interpreter', 'latex', 'fontsize', 16)
ylabel( '$\mathcal{I}(\Delta t \lambda_i$) ', 'interpreter', 'latex', 'fontsize', 16)
title( 'Crank-Nicolson Region of Absolute Stability', 'interpreter', 'latex', 'fontsize', 16)
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 14 )

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [15 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15])
set(gcf, 'PaperPosition', [0 0 15 15])
svnm = 'CN_stab';
print( '-dpng', svnm, '-r200' );