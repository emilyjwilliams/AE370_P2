%% AE 370 Project 2
% Emily Williams
% Acoustic wave equation

%% Spatial Convergence Test

clear all
close all
clc

L = 20;

a = 0;
b = L;
c = 343;
ln = b - a;
T = 20;
dt = 0.1;
u_m = 50;
w = 1;
k = 1;

uex = @(x,t) u_m.*sin(w.*t-k.*x);

% g(x)
g = @(x,t) u_m.*(k.^2.*c.^2 - w.^2).*sin(w.*t-k.*x);

% u(x=a,t) = g_a(t)
g_a = @(t) u_m.*sin(w.*t-k.*a);

% u(x=b,t) = g_b(t)
g_b = @(t) u_m.*sin(w.*t-k.*b);

% initial condition
ueta = @(x) uex(x,0);
veta = @(x) u_m.*w.*cos(-k.*x);

%# of n points to use
nvect = [20; 40; 80; 100];

%initialize error vect
err_cn = zeros( size(nvect) );
err_be = zeros( size(nvect) );

for j = 1 : length( nvect )

        n = nvect(j);

        % discretize
        xj = (a:ln/n:b)';

        dx = ln/n;
        
        % A construction
        A = (c^2/dx^2)*(diag(-2*ones(n-1,1),0) + diag(ones(n-2,1),1) + diag(ones(n-2,1),-1));
        
        B = [ zeros(size(A)) eye(size(A)) ; A zeros(size(A)) ];
        
        p = @(t) [ zeros(size(g(xj(2:end-1),t))) ; g(xj(2),t)+(c^2/dx^2).*g_a(t) ; g(xj(3:end-2),t) ; g(xj(end-1),t)++(c^2/dx^2).*g_b(t) ];

        % f(z,t)
        f = @(z,t) [ v ; A*u + g(t) ];
    %---

    %--- initialize
        zk_cn = double([ (ueta(xj(2:end-1))) ; veta(xj(2:end-1)) ]);
        zk_be = double([ (ueta(xj(2:end-1))) ; veta(xj(2:end-1)) ]);
        tk = 0;
        tvect = dt : dt : T;
        
        nsnps = 100;
        ind = max( 1, round(length(tvect)/nsnps) );
        tsv = tvect( 1 : ind : end );

        z_cn = zeros( 2*n-2, length(tsv));
        z_be = zeros( 2*n-2, length(tsv));
        cnt = 1;
    %---
    
    %--- time stepping
    for jj = 1 : length( tvect )

        stat = [ num2str(j), ' out of ', num2str(length(nvect)), ': ', num2str(jj), ' out of ', num2str(length(tvect)) ];
        disp(stat)
        
        tkp1 = tk + dt;
        
        % kp1
        zkp1_cn = (eye(size(B)) - 0.5*dt*B)\(zk_cn + 0.5*dt*B*zk_cn + 0.5*dt*(p(tk)+p(tkp1)));
        zkp1_be = (eye(size(B)) - dt*B)\(zk_be + dt*p(tkp1));
  
        % update
        zk_cn = zkp1_cn;
        zk_be = zkp1_be;
        tk = tkp1;

        if min(abs( tkp1-tsv ) ) < 1e-8
            z_cn(:,cnt) = zk_cn;
            z_be(:,cnt) = zk_be;
            cnt = cnt + 1;
        end

    end
    %---
    
    err_cn(j) = norm( vpa(zkp1_cn(1:n-1,:)) - vpa(uex(xj(2:end-1),tk)) ) / norm( vpa(uex(xj(2:end-1),tk)) );
    err_be(j) = norm( vpa(zkp1_be(1:n-1,:)) - vpa(uex(xj(2:end-1),tk)) ) / norm( vpa(uex(xj(2:end-1),tk)) );

end

figure(100)
F(nsnps) = struct('cdata',[],'colormap',[]);
anim = VideoWriter('animation','MPEG-4');
anim.FrameRate = 10;
anim.Quality = 97;
open(anim)
for i = 1:nsnps
   plot(20,zkp1_cn(i,:),'b.','MarkerSize',20);
   F(i) = getframe;
   writeVideo(anim,F(i));
end
close(anim) 

figure(1000)
for i = 1:nsnps-1
    scatter(xj(i),zkp1_cn(i,:),'filled'), hold on
end
set( gca, 'fontsize', 15, 'ticklabelinterpreter', 'latex' )
title('Crank-Nicolson Spatial Propagation', 'fontsize', 20, 'interpreter', 'latex')
xlabel('$x$', 'fontsize', 15, 'interpreter', 'latex')
ylabel('$u$', 'fontsize', 15, 'interpreter', 'latex')
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [20 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 20 15])
set(gcf, 'PaperPosition', [0 0 20 15])
svnm = 'animation';
print( '-dpng', svnm, '-r200' );


%-- waterfall
    [X,T] = meshgrid( xj(2:end-1), tsv );

    figure(1), subplot(1,2,1)
    waterfall( X,T, uex(X,T) ), hold on
    set( gca, 'fontsize', 15, 'ticklabelinterpreter', 'latex' )
    title('$u_{exact}$', 'fontsize', 20, 'interpreter', 'latex')
    xlabel('$x$', 'fontsize', 15, 'interpreter', 'latex')
    ylabel('$t$', 'fontsize', 15, 'interpreter', 'latex')
    zlim([-80 80])

    subplot(1,2,2)
    y1 = z_cn';
    waterfall( X,T, y1(:,1:n-1) ), hold on
    set( gca, 'fontsize', 15, 'ticklabelinterpreter', 'latex' )
    title({'Crank-Nicolson','$u_{FD}$'}, 'fontsize', 20, 'interpreter', 'latex')
    xlabel('$x$', 'fontsize', 15, 'interpreter', 'latex')
    ylabel('$t$', 'fontsize', 15, 'interpreter', 'latex')
    zlim([-80 80])
    
    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'Color', [1 1 1])
    set(gca, 'Color', [1 1 1])
    set(gcf, 'PaperUnits', 'centimeters')
    set(gcf, 'PaperSize', [20 15])
    set(gcf, 'Units', 'centimeters' )
    set(gcf, 'Position', [0 0 20 15])
    set(gcf, 'PaperPosition', [0 0 20 15])
    svnm = 'waterfall';
    print( '-dpng', svnm, '-r200' );

    figure(2), subplot(1,2,1)
    waterfall( X,T, uex(X,T) ), hold on
    set( gca, 'fontsize', 15, 'ticklabelinterpreter', 'latex' )
    title('$u_{exact}$', 'fontsize', 20, 'interpreter', 'latex')
    xlabel('$x$', 'fontsize', 15, 'interpreter', 'latex')
    ylabel('$t$', 'fontsize', 15, 'interpreter', 'latex')
    zlim([-80 80])

    subplot(1,2,2)
    y2 = z_be';
    waterfall( X,T, y2(:,1:n-1) ), hold on
    set( gca, 'fontsize', 15, 'ticklabelinterpreter', 'latex' )
    title({'Backward Euler','$u_{FD}$'}, 'fontsize', 20, 'interpreter', 'latex')
    xlabel('$x$', 'fontsize', 15, 'interpreter', 'latex')
    ylabel('$t$', 'fontsize', 15, 'interpreter', 'latex')
    zlim([-80 80])
    
    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'Color', [1 1 1])
    set(gca, 'Color', [1 1 1])
    set(gcf, 'PaperUnits', 'centimeters')
    set(gcf, 'PaperSize', [20 15])
    set(gcf, 'Units', 'centimeters' )
    set(gcf, 'Position', [0 0 20 15])
    set(gcf, 'PaperPosition', [0 0 20 15])
%--

%-- error
    figure(4)
    loglog( ln./nvect, err_cn , 'b.-', 'markersize', 26, 'linewidth', 2 ), hold on
    loglog( ln./nvect, err_be , 'r.-', 'markersize', 26, 'linewidth', 2 )
    h = legend('CN','BE');
    set(h, 'Interpreter','latex', 'fontsize', 16, 'Location', 'NorthWest' )
    
    xlabel( '$\Delta x$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel( '$||\textbf{e}||/||\textbf{u}_e||$ ', 'interpreter', 'latex', 'fontsize', 16)

    set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'Color', [1 1 1])
    set(gca, 'Color', [1 1 1])
    set(gcf, 'PaperUnits', 'centimeters')
    set(gcf, 'PaperSize', [15 15])
    set(gcf, 'Units', 'centimeters' )
    set(gcf, 'Position', [0 0 15 15])
    set(gcf, 'PaperPosition', [0 0 15 15])
    svnm = 'spat';
    print( '-dpng', svnm, '-r200' );
%--

%% Temporal Convergence Test CN

clear all
close all
clc

L = 10;

a = 0;
b = L;
c = 5;
ln = b - a;
T = 20;
u_m = 50;
w = 1;
k = 1;

uex = @(x,t) u_m.*sin(w.*t-k.*x);

n = 3000; 
dx = (b-a)/n;  

% g(x)
g = @(x,t) u_m.*(k.^2.*c.^2 - w.^2).*sin(w.*t-k.*x);

% u(x=a,t) = g_a(t)
g_a = @(t) u_m.*sin(w.*t-k.*a);

% u(x=b,t) = g_b(t)
g_b = @(t) u_m.*sin(w.*t-k.*b);

% initial condition
ueta = @(x) uex(x,0);
veta = @(x) u_m.*w.*cos(-k.*x);

% dt's
dtvect = [0.8e-1; 0.62e-1; 0.68e-1; 0.5e-1];

% error vect
err_cn = zeros( size( dtvect ) );

for j = 1 : length( dtvect )

        dt = dtvect(j);

        % discretize
        xj = (a:ln/n:b)';

        % A construction
        A = (c^2/dx^2)*(diag(-2*ones(n-1,1),0) + diag(ones(n-2,1),1) + diag(ones(n-2,1),-1));
        
        B = [ zeros(size(A)) eye(size(A)) ; A zeros(size(A)) ];
        
        p = @(t) [ zeros(size(g(xj(2:end-1),t))) ; g(xj(2),t)+(c^2/dx^2).*g_a(t) ; g(xj(3:end-2),t) ; g(xj(end-1),t)++(c^2/dx^2).*g_b(t) ];

        % f(z,t)
        f = @(z,t) [ v ; A*u + g(t) ];
    %---

    %--- initialize
        zk_cn = double([ ueta(xj(2:end-1)) ; veta(xj(2:end-1)) ]);
        tk = 0;
        tvect = dt : dt : T;
    %---
    
        M = inv(eye(size(B)) - 0.5*dt*B);
    
    %--- time stepping
    for jj = 1 : length( tvect )
        
        tkp1 = tk + dt;
        
        % kp1
        zkp1_cn = M*(zk_cn + 0.5*dt*B*zk_cn + 0.5*dt*(p(tk)+p(tkp1)));
  
        % update
        zk_cn = double(zkp1_cn);
        tk = tkp1;

    end
    %---
    
    stat = [ num2str(j), ' out of ', num2str(length(dtvect)) ];
    disp(stat)
    
    err_cn(j) = norm( vpa(zkp1_cn(1:n-1,:)) - vpa(uex(xj(2:end-1),tk)) ) / norm( vpa(uex(xj(2:end-1),tk)) );

end


%-- error plot
    figure(2)
    m = err_cn(end)*(1./dt^2);
    loglog( dtvect, m*(dtvect).^2, 'k--', 'linewidth', 2 ), hold on
    loglog( dtvect, err_cn , 'b.', 'markersize', 26 )
    ylim([1e-5 1.5e-3])
    xlim([0.049 0.081])
    yticks(1e-7:0.5e-1:1e-4)
    
    % formatting
    h = legend('$O(\Delta t^2)$', 'CN');
    set(h, 'Interpreter','latex', 'fontsize', 16, 'Location', 'NorthWest' )
    
    xlabel( '$\Delta t$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel( '$||\textbf{e}||/||\textbf{u}_e||$ ', 'interpreter', 'latex', 'fontsize', 16)

    set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'Color', [1 1 1])
    set(gca, 'Color', [1 1 1])
    set(gcf, 'PaperUnits', 'centimeters')
    set(gcf, 'PaperSize', [15 15])
    set(gcf, 'Units', 'centimeters' )
    set(gcf, 'Position', [0 0 15 15])
    set(gcf, 'PaperPosition', [0 0 15 15])
    svnm = 'temp';
    print( '-dpng', svnm, '-r200' );
%--

%% Temporal Convergence Test BE

clear all
close all
clc

L = 10;

a = 0;
b = L;
c = 5;
ln = b - a;
T = 20;
u_m = 50;
w = 1;
k = 1;

uex = @(x,t) u_m.*sin(w.*t-k.*x);

n = 3000; 
dx = (b-a)/n;  

% g(x)
g = @(x,t) u_m.*(k.^2.*c.^2 - w.^2).*sin(w.*t-k.*x);

% u(x=a,t) = g_a(t)
g_a = @(t) u_m.*sin(w.*t-k.*a);

% u(x=b,t) = g_b(t)
g_b = @(t) u_m.*sin(w.*t-k.*b);

% initial condition
ueta = @(x) uex(x,0);
veta = @(x) u_m.*w.*cos(-k.*x);

% dt's
dtvect = [0.8e-1; 0.62e-1; 0.68e-1; 0.5e-1];

% error vect
err_be = zeros( size( dtvect ) );

for j = 1 : length( dtvect )

        dt = dtvect(j);

        % discretize
        xj = (a:ln/n:b)';

        % A construction
        A = (c^2/dx^2)*(diag(-2*ones(n-1,1),0) + diag(ones(n-2,1),1) + diag(ones(n-2,1),-1));
        
        B = [ zeros(size(A)) eye(size(A)) ; A zeros(size(A)) ];
        
        p = @(t) [ zeros(size(g(xj(2:end-1),t))) ; g(xj(2),t)+(c^2/dx^2).*g_a(t) ; g(xj(3:end-2),t) ; g(xj(end-1),t)++(c^2/dx^2).*g_b(t) ];

        % f(z,t)
        f = @(z,t) [ v ; A*u + g(t) ];
    %---

    %--- initialize
        zk_be = double([ ueta(xj(2:end-1)) ; veta(xj(2:end-1)) ]);
        tk = 0;
        tvect = dt : dt : T;
    %---
    
        M = inv(eye(size(B)) - dt*B);
    
    %--- time stepping
    for jj = 1 : length( tvect )
        
        tkp1 = tk + dt;
        
        % kp1
        zkp1_be = M*(zk_be + dt*p(tkp1));;
  
        % update
        zk_be = double(zkp1_be);
        tk = tkp1;

    end
    %---
    
    stat = [ num2str(j), ' out of ', num2str(length(dtvect)) ];
    disp(stat)
    
    err_be(j) = norm( vpa(zkp1_be(1:n-1,:)) - vpa(uex(xj(2:end-1),tk)) ) / norm( vpa(uex(xj(2:end-1),tk)) );

end


%-- error plot
    figure(2)
    m = err_be(end)*(1./dt);
    loglog( dtvect, m*(dtvect), 'k--', 'linewidth', 2 ), hold on
    loglog( dtvect, err_be , 'r.', 'markersize', 26 )
    yticks(1e-7:0.5e-1:1e-4)
    
    % formatting
    h = legend('$O(\Delta t)$', 'BE');
    set(h, 'Interpreter','latex', 'fontsize', 16, 'Location', 'NorthWest' )
    
    xlabel( '$\Delta t$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel( '$||\textbf{e}||/||\textbf{u}_e||$ ', 'interpreter', 'latex', 'fontsize', 16)

    set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'Color', [1 1 1])
    set(gca, 'Color', [1 1 1])
    set(gcf, 'PaperUnits', 'centimeters')
    set(gcf, 'PaperSize', [15 15])
    set(gcf, 'Units', 'centimeters' )
    set(gcf, 'Position', [0 0 15 15])
    set(gcf, 'PaperPosition', [0 0 15 15])
    svnm = 'temp_be';
    print( '-dpng', svnm, '-r200' );
%--


%% Frequency Analysis

clear all
close all
clc

L = 10;

a = 0;
b = L;
c = 343;
ln = b - a;
T = 10;
u_m = 50; % dB
freq = [ 261.6256; 523.2511; 391.9954; 739.9888; 1108.731; 7040; 3520; 1760; 880; 440; 220; 110; 55 ];
n = 1000; 
dx = (b-a)/n;  
dt = 0.1;

for j = 1:length(freq)
    
    w = 2*pi*freq(j);
    lambda = c/freq(j);
    k = 2*pi/lambda;

    uex = @(x,t) u_m.*sin(w.*t-k.*x);

    % g(x)
    g = @(x,t) u_m.*(k.^2.*c.^2 - w.^2).*sin(w.*t-k.*x);

    % u(x=a,t) = g_a(t)
    g_a = @(t) u_m.*sin(w.*t-k.*a);

    % u(x=b,t) = g_b(t)
    g_b = @(t) u_m.*sin(w.*t-k.*b);

    % initial condition
    ueta = @(x) uex(x,0);
    veta = @(x) u_m.*w.*cos(-k.*x);

    % discretize
    xj = (a:ln/n:b)';

    % A construction
    A = (c^2/dx^2)*(diag(-2*ones(n-1,1),0) + diag(ones(n-2,1),1) + diag(ones(n-2,1),-1));

    B = [ zeros(size(A)) eye(size(A)) ; A zeros(size(A)) ];

    p = @(t) [ zeros(size(g(xj(2:end-1),t))) ; g(xj(2),t)+(c^2/dx^2).*g_a(t) ; g(xj(3:end-2),t) ; g(xj(end-1),t)++(c^2/dx^2).*g_b(t) ];

    % f(z,t)
    f = @(z,t) [ v ; A*u + g(t) ];
    %---

    %--- initialize
    zk_cn = double([ ueta(xj(2:end-1)) ; veta(xj(2:end-1)) ]);
    tk = 0;
    tvect = dt : dt : T;
    %---

    M = inv(eye(size(B)) - 0.5*dt*B);

    %--- time stepping
    for jj = 1 : length( tvect )

        tkp1 = tk + dt;

        % kp1
        zkp1_cn = M*(zk_cn + 0.5*dt*B*zk_cn + 0.5*dt*(p(tk)+p(tkp1)));

        % update
        zk_cn = double(zkp1_cn);
        tk = tkp1;

    end
    %---
    
    z_cn(:,j) = zk_cn(1:n-1);
    
end


%-- plot
    figure(2)
    plot( z_cn(:,2), '-','color',[0.8500 0.3250 0.0980], 'linewidth', 2 ), hold on
    plot( z_cn(:,4), 'b-', 'linewidth', 2 ), hold on
    plot( z_cn(:,5), 'k-', 'linewidth', 2 )
    
    % formatting
    h = legend('A$_7$ (3520 Hz)','A$_5$ (880 Hz)','A$_4$ (440 Hz)','interpreter','latex');
    set(h, 'Interpreter','latex', 'fontsize', 16, 'Location', 'NorthWest' )
    
    xlabel( '$x$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel( '$u$ ', 'interpreter', 'latex', 'fontsize', 16)

    set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'Color', [1 1 1])
    set(gca, 'Color', [1 1 1])
    set(gcf, 'PaperUnits', 'centimeters')
    set(gcf, 'PaperSize', [25 15])
    set(gcf, 'Units', 'centimeters' )
    set(gcf, 'Position', [0 0 25 15])
    set(gcf, 'PaperPosition', [0 0 25 15])
%     svnm = 'a2';
%     print( '-dpng', svnm, '-r200' );

    figure(6)
    plot( z_cn(:,1), '-','color',[0.8500 0.3250 0.0980],'linewidth', 2 ), hold on
    plot( z_cn(:,3), 'b-', 'linewidth', 2 )
    
    % formatting
    h = legend('A$_8$ (7040 Hz)','A$_6$ (1760 Hz)','interpreter','latex');
    set(h, 'Interpreter','latex', 'fontsize', 16, 'Location', 'SouthWest' )
    
    xlabel( '$x$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel( '$u$ ', 'interpreter', 'latex', 'fontsize', 16)

    set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'Color', [1 1 1])
    set(gca, 'Color', [1 1 1])
    set(gcf, 'PaperUnits', 'centimeters')
    set(gcf, 'PaperSize', [25 15])
    set(gcf, 'Units', 'centimeters' )
    set(gcf, 'Position', [0 0 25 15])
    set(gcf, 'PaperPosition', [0 0 25 15])
%     svnm = 'a3';
%     print( '-dpng', svnm, '-r200' );

    figure(3)
    plot( z_cn(:,3), 'b-', 'linewidth', 2 ), hold on
    plot( z_cn(:,4), 'r-', 'linewidth', 2 )
    
    % formatting
    h = legend('C$_4$ (262 Hz)','C$_5$ (523 Hz)','interpreter','latex');
    set(h, 'Interpreter','latex', 'fontsize', 16, 'Location', 'NorthWest' )
    
    xlabel( '$x$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel( '$u$ ', 'interpreter', 'latex', 'fontsize', 16)

    set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'Color', [1 1 1])
    set(gca, 'Color', [1 1 1])
    set(gcf, 'PaperUnits', 'centimeters')
    set(gcf, 'PaperSize', [15 15])
    set(gcf, 'Units', 'centimeters' )
    set(gcf, 'Position', [0 0 15 15])
    set(gcf, 'PaperPosition', [0 0 15 15])
%     svnm = 'C';
%     print( '-dpng', svnm, '-r200' );

    figure(4)
    plot( z_cn(:,5), 'b-', 'linewidth', 2 ), hold on
    plot( z_cn(:,4), 'r-', 'linewidth', 2 )
    
    % formatting
    h = legend('G$_4$ (391.9954 Hz)','interpreter','latex');
    set(h, 'Interpreter','latex', 'fontsize', 16, 'Location', 'NorthWest' )
    
    xlabel( '$x$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel( '$u$ ', 'interpreter', 'latex', 'fontsize', 16)

    set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'Color', [1 1 1])
    set(gca, 'Color', [1 1 1])
    set(gcf, 'PaperUnits', 'centimeters')
    set(gcf, 'PaperSize', [15 15])
    set(gcf, 'Units', 'centimeters' )
    set(gcf, 'Position', [0 0 15 15])
    set(gcf, 'PaperPosition', [0 0 15 15])
%     svnm = 'G4';
%     print( '-dpng', svnm, '-r200' );

    figure(5)
    plot( z_cn(:,6), 'b-', 'linewidth', 2 ), hold on
    plot( z_cn(:,7), 'r-', 'linewidth', 2 )
    
    % formatting
    h = legend('F$\sharp_5$ (739.9888 Hz)','D$\flat_6$ (1108.731 Hz)','interpreter','latex');
    set(h, 'Interpreter','latex', 'fontsize', 16, 'Location', 'NorthWest' )
    
    xlabel( '$x$', 'interpreter', 'latex', 'fontsize', 16)
    ylabel( '$u$ ', 'interpreter', 'latex', 'fontsize', 16)

    set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'Color', [1 1 1])
    set(gca, 'Color', [1 1 1])
    set(gcf, 'PaperUnits', 'centimeters')
    set(gcf, 'PaperSize', [15 15])
    set(gcf, 'Units', 'centimeters' )
    set(gcf, 'Position', [0 0 15 15])
    set(gcf, 'PaperPosition', [0 0 15 15])
%     svnm = 'fs5';
%     print( '-dpng', svnm, '-r200' );
%--



