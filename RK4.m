%% Spatial Convergence Test RK4

clear all
close all
clc

L = 20;

a = 0;
b = L;
c = 343;
ln = b - a;
T = 0.25;
dt = 0.0001;
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
nvect = [80; 100; 150; 200];

%initialize error vect
err_rk4 = zeros( size(nvect) );

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
        zk_rk4 = double([ (ueta(xj(2:end-1))) ; veta(xj(2:end-1)) ]);
        tk = 0;
        tvect = dt : dt : T;
        
        nsnps = 100;
        ind = max( 1, round(length(tvect)/nsnps) );
        tsv = tvect( 1 : ind : end );

        z_rk4 = zeros( 2*n-2, length(tsv));
        cnt = 1;
    %---
    
    %--- time stepping
    for jj = 1 : length( tvect )

        stat = [ num2str(j), ' out of ', num2str(length(nvect)), ': ', num2str(jj), ' out of ', num2str(length(tvect)) ];
        disp(stat)
        
        tkp1 = tk + dt;
        
        % kp1
        y1 = B*zk_rk4 + p(tk);
        y2 = B*zk_rk4 + 0.5*dt*B*y1 + p(tk+0.5*dt);
        y3 = B*zk_rk4 + 0.5*dt*B*y2 + p(tk+0.5*dt);
        y4 = B*zk_rk4 + dt*B*y3 + p(tk+dt);
        zkp1_rk4 = zk_rk4 + (1/6)*dt*(y1+2*y2+2*y3+y4);
  
        % update
        zk_rk4 = zkp1_rk4;
        tk = tkp1;

        if min(abs( tkp1-tsv ) ) < 1e-8
            z_rk4(:,cnt) = zk_rk4;
            cnt = cnt + 1;
        end

    end
    %---
    
    err_rk4(j) = norm( vpa(zkp1_rk4(1:n-1,:)) - vpa(uex(xj(2:end-1),tk)) ) / norm( vpa(uex(xj(2:end-1),tk)) );

end

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
    y1 = z_rk4';
    waterfall( X,T, y1(:,1:n-1) ), hold on
    set( gca, 'fontsize', 15, 'ticklabelinterpreter', 'latex' )
    title({'RK4','$u_{FD}$'}, 'fontsize', 20, 'interpreter', 'latex')
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
    m = err_rk4(end)/(dx^2);
    loglog( ln./nvect, m*(ln./nvect).^2, 'k--', 'linewidth', 2 ), hold on
    loglog( ln./nvect, err_rk4 , 'g.', 'markersize', 26, 'linewidth', 2 ), hold on
    h = legend('$O(\Delta x^2)$', 'RK4');
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
    svnm = 'spat_rk4';
    print( '-dpng', svnm, '-r200' );
%--

