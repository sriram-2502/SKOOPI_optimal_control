function plot_figure(robot_params,navigation_params,euler_params,joint_control,joints)
%% get params
a1 = robot_params.l1 ; a2 = robot_params.l2;
w = 0.05; % width of robotic arm
% joint_angles = joint_control.angles;
% joint_vel = joint_control.vel;
% control = joint_control.control;
time = joint_control.t;
time = time(1:end-2);
x_goal = navigation_params.x_goal_f;
colors = colororder;
blue = colors(1,:);
red = colors(2,:);
yellow = colors(3,:);
green = colors(5,:);
grayColor = [.7 .7 .7];
redColor = [1 0 0];
obsColor = [.7 .7 .7]; % Obstacle color -> Grey
figure(111)
hold on;
%% ------------------------ plot task space solution --------------------------------------------
p1 = subplot(2,4,2);
% get links config
q1 = joints(1,1); q2 = joints(1,2);
x1 = a1*sin(q1); y1 = -a1*cos(q1);
x2 = x1 + a2*sin(q1+q2); y2 = y1 - a2*cos(q1+q2);
idx = [5:75:500];
% idx = [5:50:200];
% idx = [200:50:300];
% idx = [300:50:400];
% idx = [400:50:498];
for angles = 1:length(idx)
    q1 = joints(idx(angles),1); q2 = joints(idx(angles),2);
    x1 = a1*sin(q1); y1 = -a1*cos(q1);
    x2 = x1 + a2*sin(q1+q2); y2 = y1 - a2*cos(q1+q2);
    
    link1 = polyshape([x1+w*cos(q1) x1-w*cos(q1) -w*cos(q1) +w*cos(q1)], ...
    [y1+w*sin(q1) y1-w*sin(q1) -w*sin(q1) w*sin(q1)]);
    plot(link1,'FaceColor',red, 'FaceAlpha',0.8^(length(idx)-angles), 'EdgeColor','black'); hold on
        
    link2 = polyshape([x2+w*cos(q1+q2) x2-w*cos(q1+q2) x1-w*cos(q1+q2) x1+w*cos(q1+q2)], ...
    [y2+w*sin(q1+q2) y2-w*sin(q1+q2) y1-w*sin(q1+q2) y1+w*sin(q1+q2)]);
    plot(link2,'FaceColor',red, 'FaceAlpha',0.8^(length(idx)-angles), 'EdgeColor','black'); hold on
end
%% plot goal position on the end effector
N = length(time)-1; 
skip_rate=50;
dT = euler_params.step_size;
state_index = idx(end);
trail_length = 400; 
tail_grad = linspace(0,skip_rate, trail_length);
if trail_length == 1
    trail_skip_rate = 1;
else
    trail_skip_rate = round(length(tail_grad)/10);
end
xc = 0.0;  yc = 0.0;    % circle center (m)
r  = 1.5;               % radius (m)
omega = 0.2;               % rad/s
x_goal = [xc + r*cos(omega*time); yc + r*sin(omega*time)]';
scatter(x_goal(state_index,1), x_goal(state_index,2), 80,'filled','o','MarkerFaceColor',green, 'LineWidth', 2); hold on;
s_plot = scatter(x_goal(state_index-trail_length:skip_rate:state_index,1),...
    x_goal(state_index-trail_length:skip_rate:state_index,2), 80,'filled','o','MarkerFaceColor',green, 'LineWidth', 2); hold on;
s_plot.AlphaData = tail_grad(1:skip_rate:trail_length)./45;
s_plot.MarkerFaceAlpha = 'flat';
plot(x_goal(state_index,1), x_goal(state_index,2), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor',green); hold on;
hold on;
xlabel('position, $x_1$','interpreter','latex', 'FontSize', 20);
ylabel('position, $x_2$','interpreter','latex', 'FontSize', 20);
% p1.XLim = [-0.5, 2.1]; p1.YLim = [-1.0, 1.6];
p1.XLim = [-1.0, 2.0]; p1.YLim = [-2.0, 0.5];
% plot optionss
axes1 = gca;
box(axes1,'on');
axis(axes1,'square');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontSize',15,'LineWidth',1.5);
% legend
subplot(2,4,3)
dummy_robot = plot(1,1,'Color',red,'LineWidth',4); hold on;
dummy_goal = plot(NaN,NaN,'o','MarkerSize',10,'MarkerEdgeColor', green, 'MarkerFaceColor',green); hold on;
lgd = legend('Robot','Goal', ...
        'Location', 'southwest','Interpreter','Latex');