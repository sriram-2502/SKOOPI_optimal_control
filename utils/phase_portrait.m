%% phase portrait

plot_phase_portrait = true;

n = 2;
x = sym('x',[n,1],'real');
f_x = [x(1)-x(1).^3 - x(2); -x(2) - x(2).^3];
f = matlabFunction(f_x,'vars',{x}); 

if(plot_phase_portrait)
    % linearization at (0,0) saddle
    Dom = [-5 5];
    [X, Y] = meshgrid(Dom(1):0.1:Dom(2), Dom(1):0.1:Dom(2));
    
    % Initialize components of the vector field
    u = zeros(size(X));
    v = zeros(size(Y));
    
    % Evaluate f at each grid point
    for i = 1:numel(X)
        xy = [X(i); Y(i)];   % 2x1 input
        result = f(xy);      % 2x1 output
        u(i) = result(1);    % x-component
        v(i) = result(2);    % y-component
    end
    
    figure(1)
    l = streamslice(X,Y,u,v); hold on;
    set(l,'LineWidth',1)
    set(l,'Color','k');
    
    % tspan = [0,5]; ic_pts = 1;
    % xl = Dom(1); xh = Dom(2);
    % yl = Dom(1); yh = Dom(2);
    % for x0 = linspace(Dom(1), Dom(2), ic_pts)
    %     for y0 = linspace(Dom(1), Dom(2), ic_pts)
    %         [ts,xs] = ode45(@(t,x)f(x),tspan,[x0 y0]');
    %         plot(xs(:,1),xs(:,2),'k','LineWidth',1); hold on;
    %     end
    % end
    xlim([-5,5])
    ylim([-5,5])
    axes = gca;
    axis square
    set(axes,'FontSize',15);
    xlabel('$x_1$','FontSize',20, 'Interpreter','latex')
    ylabel('$x_2$','FontSize',20, 'Interpreter','latex')
    box on
    axes.LineWidth=2;
end