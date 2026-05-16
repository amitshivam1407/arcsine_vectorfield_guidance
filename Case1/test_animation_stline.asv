%-------------------------------------------------------------------------%
%------------                    14th April 2022            --------------%
%------------                   Autopilot design            --------------%
%------------            Straight line path following    -----------------%
%-------------------------------------------------------------------------%

close all;clear all; clc

%--------------- declaration of constant parameters ----------------------%
global  Vgd kv k_kai k_kaidot k2 figure_type_good
Vgd = 10;
kv = 5;
k_kai = 30;
k_kaidot = sqrt(k_kai);  % 20 for better curvature profile  %%
figure_type_good = 1;

%-------------------------------------------------------------------------%

%------------------ initial conditions   ---------------------------------%

tspan = [0 25];

x0 = 100;
y0 = -100;
k2 = 0.002;
if x0 < 0
    kaid0 =  asin(1./(1+ k2 *(x0).^2)) ;
    %    kaidot_des(i,1) = 2*k2*Vgd.*(x(i,1))./((1+k2*(x(i,1)).^2).^2);
else
    kaid0 = pi - asin(1./(1 + k2 *(x0).^2)) ;
    %   kaidot_des(i,1) = - 2*k2*Vgd.*(x(i,1))./((1 + k2*(x(i,1)).^2).^2);
end
kai0 = kaid0;
kaidot0 = -(2*k2*Vgd.*x0)./((1+k2*(x0).^2).^2);

% k2 = abs(14.92597/x0^2);

x_initial = [x0;y0;kai0;kaidot0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%  vector field construction  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
range  = -110:10:110;
[X,Y] = meshgrid(range);
[kaid,kaid_dot] = vf_proposed(X);
xdot  = Vgd*cos(kaid);
ydot  = Vgd*sin(kaid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = odeset('RelTol',1e-8,'AbsTol',1e-8);

%%%%%%%%%%%%%%%%%  ode solver   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,x] = ode45(@(t,x)fun_stline_upd(t,x) ,tspan, x_initial,options);

%%%%%%%%%%%%%%%%%%%%%%%  parameters to plot   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psi = x(:,3) + asin((wx*sin(x(:,3)) - wy*cos(x(:,3)))./ x(:,5) );
% Vg_act = sqrt(x(:,5).^2 + wx^2 + wy^2 + 2*wx*x(:,5).*cos(psi) + 2*wy*x(:,5).*sin(psi));
x_ini = x(1,1) ;
y_ini = x(1,2) ;
x_end = x(end,1) ;
y_end = x(end,2) ;
for i = 1:length(t)
    if x(i,1) < 0
        kai_des(i,1) =  asin(1./(1+ k2*(x(i,1)).^2)) ;
        kaidot_des(i,1) = 2*k2*Vgd.*(x(i,1))./((1+k2*(x(i,1)).^2).^2);
    else
        kai_des(i,1) = pi - asin(1./(1 + k2*(x(i,1)).^2)) ;
        kaidot_des(i,1) = - 2*k2*Vgd.*(x(i,1))./((1 + k2*(x(i,1)).^2).^2);
    end
end
kaidot_actual = x(:,4);
kappa_des = kaidot_des./Vgd;
kappa_actual = x(:,4)./Vgd ;
% [kai_des,kaidot_des] = vf_proposed(x(:,1));

% Va_des = sqrt((Vgd*cos(kai_des) - wx).^2 + (Vgd*sin(kai_des) - wy).^2) ;
% psi_des = kai_des + asin((wx*sin(kai_des) - wy*cos(kai_des))./ Va_des );
% % psiddot_des = kaidot_des.*((Va_des.*cos(psi_des - kai_des) + (wx*cos(kai_des) + wy*sin(kai_des)))./(Va_des.*cos(psi_des - kai_des)));
% c1 = (1./x(:,5)).*kaidot_des.*(wx*cos(kai_des)+wy*sin(kai_des));
% c2 = (wx*sin(kai_des)-wy*cos(kai_des)).*((1./x(:,5)).^2).*(Vgd./Va_des).*kaidot_des.*(-wx*sin(kai_des)+wy*cos(kai_des));
% psiddot_des = kaidot_des + (1./sqrt(1 - ((1./x(:,5)).*(wx*sin(kai_des)-wy*cos(kai_des))).^2)).*(c1 + c2);
% psi_actual = x(:,3);
% psidot_actual = x(:,4);
% Vg_actual = sqrt((x(:,5).*cos(x(:,3)) + wx).^2 + (x(:,5).*sin(x(:,3)) + wy).^2) ;
% % Vg_actual = (wx.*cos(x(:,3)) + wy.*sin(x(:,3))) + sqrt(x(:,5).^2 - (wx.*sin(x(:,3)) - wy.*cos(x(:,3))).^2);
% kai_actual = x(:,3) + asin((-wx*sin(x(:,3)) + wy*cos(x(:,3)))./ Vg_actual );
% kaidot_actual = x(:,4).*((Vg_actual.*cos(kai_actual - x(:,3)) - (wx*cos(x(:,3)) + wy*sin(x(:,3))))./(Vg_actual.*cos(kai_actual - x(:,3))));
% kappa_actual = kaidot_actual./Vg_actual;
%---------------------- plotting figures ---------------------------------%
%
% figure_plot = fun_stlineplot(t,X,Y,xdot,ydot,x(:,1),x(:,2),x_ini,y_ini,x_end,y_end,...
%     kaidot_des,kaidot_actual,kappa_des,kappa_actual,kai_des,x(:,3));

% %% Animate Trajectory (used variables: x, y, si)
%
% fig = gcf;
% fig.WindowState = 'maximized';
% %pause(1)
%
% isACmarker = false;
% markerFileName = 'plane.jpeg';
% delayTime = 0.08;
%
% noLOS = min(60,length(x));
% jumpStep = floor(length(x)/10);
%
% %% Setup Video creation
% vidFlnm = 'AnimatePlot';
%
% if(ischar(vidFlnm))
%
% 	vid = VideoWriter(vidFlnm,'MPEG-4');
% 	vid.Quality = 100;
% 	vid.FrameRate = 15;
%
% 	open(vid)
%
% end
%
% %% Plot the initial coordinates
% if(isACmarker)
%     % Find the initial coordinates of the marker WITHOUT orientation
%     [xm0, ym0] = setupPlotMarker(markerFileName, 1/4, 0.9);
%
%     % Find the initial coordinates of the marker WITH orientation
%     [xmM, ymM] = calcMarkerCoord(xm0,ym0,x(1,3));
%     %[xmT, ymT] = calcMarkerCoord(xm0,ym0,alt(1));
%
%     h1 = plot(xmM+x(1,1),ymM+x(1,2),'s','MarkerSize',0.8,'Color','blue','MarkerFaceColor','blue');
%     %h2 = plot(xmT+XT(1),ymT+YT(1),'o','MarkerSize',1.3,'Color','red','MarkerFaceColor','red');
% else
%     h1 = plot(x(1,1),x(1,2),'ob','MarkerSize',8,'MarkerFaceColor','blue');
%    % h2 = plot(XT(1),YT(1),'or','MarkerSize',8,'MarkerFaceColor','red');
% end
%
%
% %% Animate the marker
%
% for i = 1+jumpStep:jumpStep:length(x)
%
%     pause(delayTime)
%
%     plot(x(i-jumpStep:i,1),x(i-jumpStep:i,2),'-b','LineWidth',3)
%     %plot(XT(i-jumpStep:i),YT(i-jumpStep:i),'-r','LineWidth',1.5)
%
%     % Marker Position Update
%     if(isACmarker)
%         [xmM, ymM] = calcMarkerCoord(xm0,ym0,x(i,3));
%         %[xmT, ymT] = calcMarkerCoord(xm0,ym0,alt(i));
%         set(h1,'XData',xmM+x(i,1),'YData',ymM+x(i,2))
%         %set(h2,'XData',xmT+XT(i),'YData',ymT+YT(i))
%     else
%         set(h1,'XData',x(i,1),'YData',x(i,2))
%         %set(h2,'XData',XT(i),'YData',YT(i))
%     end
%
%     drawnow
%
%     if(ischar(vidFlnm))
%
% 		F = getframe(gcf);
% 		writeVideo(vid,F)
%
%     end
%
% end
%
% plot(x(i-jumpStep:end,1),x(i-jumpStep:end,2),'-b','LineWidth',3)
% %plot(XT(i-jumpTime:end),YT(i-jumpTime:end),'-r','LineWidth',1.5)
%
% if(isACmarker)
%     [xmM, ymM] = calcMarkerCoord(xm0,ym0,x(end,3));
%     %[xmT, ymT] = calcMarkerCoord(xm0,ym0,alt(end));
%     set(h1,'XData',xmM+x(end,1),'YData',ymM+x(end,2))
%     %set(h2,'XData',xmT+XT(end),'YData',ymT+YT(end))
% else
%     set(h1,'XData',x(end,1),'YData',x(end,2))
%     %set(h2,'XData',XT(end),'YData',YT(end))
% end
%
% F = getframe(gcf);
% writeVideo(vid,F)
%
%  close(vid)

%  myVideo = VideoWriter('poly1','MPEG-4'); myVideo.FrameRate =5;
%  open(myVideo);

%% Animate Trajectory (used variables: x, y, si)
% plot parameter
ax_fnt = 23;
lbl_fnt = 25;
ax_wdth = 3;
lgd_fnt = 17;

%  fig = gcf;
% f1 = figure
% fig.WindowState = 'maximized';

%  f1 = figure;
%  f1.Position = [10 50 1000 950];
fig1 = figure('Position', [0 0 1920 1080]);
% figure
% subplot(2,2,1)
ax1 = gca;
% ax1.Position(1) = [10 50 1000 950];
h = quiver(ax1,X,Y,xdot,ydot,'color',[0.75  0.75   0.75],'linewidth',1);hold(ax1,'on');
xl = xline(ax1,0,'k','linewidth',3);hold(ax1,'on');
h1 = plot(ax1,x(1,1),x(1,2),'ob','MarkerSize',10,'MarkerFaceColor','blue');
% h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h11 = plot(ax1,x(:,1),x(:,2),'-r','LineWidth',4);hold on; grid on;
ax1.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax1.XColor = 'black';         % Box horizontal lines' color
ax1.YColor = 'black';         % Box vertical lines' color
set(ax1,'linewidth',5)        % Axis linewidth (box and grid)
xlabel(ax1,' $ x, $ m','Fontsize',lbl_fnt);
ylabel(ax1,'$ y, $ m','Fontsize',lbl_fnt);
% l= legend(ax1,[h xl h11],'Vector field','Desired path','UAV trajectory');hold on;
hold(ax1,'on');
axis equal;
ax1.XLim = [-100 100];
ax1.YLim = [-110 100];

% f2 = figure;
% f2.Position = [1010 570 900 430];
fig2 = figure('Position', [0 0 1920 1080]);
% subplot(2,2,2)
ax2 = gca;
% ax2.Position(2) = [1010 550 900 450];
% h21 = plot(ax2,t,x(:,1),'-r','LineWidth',4);hold on; grid on;
ax2.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax2.XColor = 'black';         % Box horizontal lines' color
ax2.YColor = 'black';         % Box vertical lines' color
set(ax2,'linewidth',5)        % Axis linewidth (box and grid)
xlabel(ax2,' $ t, $ s','Fontsize',lbl_fnt);
ylabel(ax2,'$ x, $ m','Fontsize',lbl_fnt);
hold(ax2,'on');grid(ax2,'on');
ax2.XLim = [0 t(end,1)];
ax2.YLim = [x(end,1) x(1,1)];
axis(ax2,'tight');

% f3 = figure;
% f3.Position = [1010 50 900 430];
fig3 = figure('Position', [20 20 1920 1080]);
% subplot(2,2,3)
% f3.Position = [1010 50 900 450];
ax3 = gca;
% ax3.Position(3) = [1010 50 900 450];
% h31 = plot(ax3,t,x(:,3),'-r','LineWidth',4);hold on; grid on;
ax3.FontSize = ax_fnt;
box on                        % Switch on the box around the axis
ax3.XColor = 'black';         % Box horizontal lines' color
ax3.YColor = 'black';         % Box vertical lines' color
set(ax3,'linewidth',5)        % Axis linewidth (box and grid)
xlabel(ax3,' $$ t, $$ s','Fontsize',lbl_fnt);
ylabel(ax3,'$$ \chi, $$ rad.','Fontsize',lbl_fnt);
hold(ax3,'on');grid(ax3,'on');
ax3.XLim = [0 t(end,1)];
ax3.YLim = [x(end,3) x(1,3)];
axis(ax3,'tight');




% % jumpstep = 5;
% %  figure(5)
% % delayTime = 0.1;
%
% % l = legend(ax1,[h xl h11] ,'Vector field','Desired path','UAV trajectory');hold on;
%
% % % plot(x(:,1),x(:,2))
% for i = 1+jumpStep:jumpStep:length(t)
%     %      set(h1,'XData',x(i,1),'YData',x(i,2))
%     frame = getframe(gcf); %get frame
%      writeVideo(myVideo,frame);
%     h11 = plot(ax1,x(i-jumpStep:i,1),x(i-jumpStep:i,2),'-r','LineWidth',4);hold on; grid on;
% %      h21 = plot(ax2,t(i-jumpStep:i,1),x(i-jumpStep:i,1),'-r','LineWidth',4);hold on; grid on;
% %       h31 = plot(ax3,t(i-jumpStep:i,1),x(i-jumpStep:i,3),'-r','LineWidth',4);hold on; grid on;
%
%
% %      set(h1,'XData',x(i,1),'YData',x(i,2))
%      hold off ;
%     pause(0.01)
%
%
%     %         delete(p1)
% end
% l= legend(ax1,[h xl h11],'Vector field','Desired path','UAV trajectory');hold on;
% % subplot(1,1,1)
% % subplot(2,1,1)
%
%   close(myVideo);
%  hFig1 = figure;
% ax1 = gca;
%  hPlot1 = plot(NaN,NaN);

%  hFig2 = figure;
%  ax2 = gca;
%  hPlot2 = plot(NaN,NaN);

h11 = animatedline(ax1,'Color', 'r', 'LineStyle', '-', 'LineWidth', 4) ;
h21 = animatedline(ax2,'Color', 'r', 'LineStyle', '-', 'LineWidth', 4) ;
h31 = animatedline(ax3,'Color', 'r', 'LineStyle', '-', 'LineWidth', 4) ;

% h1 = animatedline('ob','MarkerSize',10,'MarkerFaceColor','blue') ;
% set(gca, 'XLim', [0 1*pi]);
% ax1.XLim = [-100 100];
% ax1.YLim = [-110 100];
% set(gca, 'XLim', [-100 100]);
% set(gca, 'YLim', [-100 100]);

l = legend(ax1,[h xl h11],'Vector field','Desired path','UAV trajectory');hold on;
box on
f = cell(length(x),1) ;
g = cell(length(x),1) ;
gg = cell(length(x),1) ;
for i = 1 : length(x)

    set(gca, 'XLim', [-100 100],'YLim', [-100 100]);
    addpoints(h11, x(i,1), x(i,2));hold on;
    set(h1,'XData',x(i,1),'YData',x(i,2));
    drawnow
    f{i} = getframe(gcf) ;
end 
for i = 1 : length(x)
    set(gca, 'XLim', [0 25],'YLim', [x(end,1) x(1,1)]);
    addpoints(h21, t(i,1), x(i,1));hold on;
     ax2.XLim = [0 t(end,1)];
    ax2.YLim = [0 x(1,1)];
    drawnow
    gg{i} = getframe(gcf) ;
end  
for i = 1 : length(x)
    set(gca, 'XLim',  [0 25],'YLim', [x(end,3) x(1,3)]);
    addpoints(h31, t(i,1), x(i,3));hold on;
    ax3.XLim = [0 t(end,1)];
    ax3.YLim = [x(end,3) x(1,3)];

    drawnow
    g{i} = getframe(gcf) ;
end

% xlabel(ax1,' $$ x, $$ m','Fontsize',lbl_fnt);
% ylabel(ax1,'$$ y, $$ m','Fontsize',lbl_fnt);

obj = VideoWriter('Trajectory','MPEG-4');
obj.Quality = 100;
obj.FrameRate = 50;
% obj1 = VideoWriter('Error','MPEG-4');
% obj1.Quality = 100;
% obj1.FrameRate = 50;
% obj2 = VideoWriter('Course','MPEG-4');
% obj2.Quality = 100;
% obj2.FrameRate = 50;

open(obj);
% open(obj1);
% open(obj2);
for i = 1:length(x)
    writeVideo(obj, f{i}) ;
%     writeVideo(obj1, gg{i}) ;
%     writeVideo(obj2, g{i}) ;
end
obj.close();
% obj1.close();
% obj1.close();

% open(obj1);
% for i = 1:length(x)
%     writeVideo(obj1, gg{i}) ;
% end
% obj1.close();
% 
% 
% 
% open(obj2);
% for i = 1:length(x)
%     writeVideo(obj2, g{i}) ;
% end
% obj1.close();


function out = fun_stline_upd(t,x)
global Vgd k2 wx wy  kv k_kai k_kaidot


if x(1) < 0
    kaid =  asin(1./(1+ k2*(x(1)).^2)) ;
    kaid_dot = 2*k2*Vgd.*(x(1))./((1+k2*(x(1)).^2).^2);
else
    kaid = pi - asin(1./(1 + k2*(x(1)).^2)) ;
    kaid_dot = - 2*k2*Vgd.*(x(1))./((1 + k2*(x(1)).^2).^2);
end
%
% Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
%
% psid = kaid + asin((wx*sin(kaid)-wy*cos(kaid))/Vad);
% psi = x(3) + asin((wx*sin(x(3))-wy*cos(x(3)))/x(5));
%
% Vg = sqrt(x(5)^2  + wx^2 +wy^2 + 2*x(5)*(wx*cos(psi) + wy*sin(psi)));


out(1,1) = Vgd*cos(x(3));
out(2,1) = Vgd*sin(x(3));
out(3,1) = x(4) ;
out(4,1) = k_kai*(kaid - x(3)) + k_kaidot*(kaid_dot - x(4) )  ;
% out(5,1) = kv*(Vad - x(5));
end

% function out = fun_stline_upd(t,x)
% global Vgd k2 wx wy  kv k_kai k_kaidot
%
%
% if x(1) < 0
%    kaid =  asin(1./(1+ k2*(x(1)).^2)) ;
%    kaid_dot = 2*k2*Vgd.*(x(1))./((1+k2*(x(1)).^2).^2);
% else
%    kaid = pi - asin(1./(1 + k2*(x(1)).^2)) ;
%    kaid_dot = - 2*k2*Vgd.*(x(1))./((1 + k2*(x(1)).^2).^2);
% end
%
% Vad = sqrt(Vgd^2  + wx^2 +wy^2 - 2*(Vgd*wx*cos(kaid) + Vgd*wy*sin(kaid)));
%
% psid = kaid + asin((wx*sin(kaid)-wy*cos(kaid))/Vad);
% % psiddot = kaid_dot.*((Vad.*cos(psid - kaid) + (wx*cos(kaid) + wy*sin(kaid)))./Vad.*cos(psid - kaid));
% c1 = (1/x(5))*kaid_dot*(wx*cos(kaid)+wy*sin(kaid));
% c2 = (wx*sin(kaid)-wy*cos(kaid))*((1/x(5))^2)*(Vgd/Vad)*kaid_dot*(-wx*sin(kaid)+wy*cos(kaid));
% psiddot = kaid_dot + (1/sqrt(1 - ((1/x(5))*(wx*sin(kaid)-wy*cos(kaid)))^2))*(c1 + c2);
% % psi = x(3) + asin((wx*sin(x(3))-wy*cos(x(3)))/x(5));
%
% % Vg = sqrt(x(5)^2  + wx^2 +wy^2 + 2*x(5)*(wx*cos(psi) + wy*sin(psi)));
%
%
% out(1,1) = x(5)*cos(x(3))  + wx ;
% out(2,1) = x(5)*sin(x(3))  + wy ;
% out(3,1) = x(4) ;
% out(4,1) = k_kai*(psid - x(3)) + k_kaidot*(psiddot - x(4))  ;
% out(5,1) = kv*(Vad - x(5));
% end

function [kaid,kaid_dot] = vf_proposed(x)
global Vgd  k2
for i = 1:length(x)
    for j = 1:length(x)
        if x(i,j) < 0
            kaid(i,j) = asin(1./(1+k2*(x(i,j)).^2)) ;
            kaid_dot(i,j) =  2*k2*Vgd.*(x(i,j))./((1+k2*(x(i,j)).^2).^2);
        else
            kaid(i,j) = pi - asin(1./(1+k2*(x(i,j)).^2)) ;
            kaid_dot(i,j) = - 2*k2*Vgd.*(x(i,j))./((1+k2*(x(i,j)).^2).^2);
        end
    end
end
end


function [xm0, ym0] = setupPlotMarker(imgfile, factor, ms)

im = resizeImage(imgfile,factor);
%imshow(im)

im = flipud(im);
[m, n, ~] = size(im);
x = zeros(m*n,1);
y = x;

k = 1;
for i = 1:m
    for j = 1:n

        R = double(im(i,j,1))/255;
        G = double(im(i,j,2))/255;
        B = double(im(i,j,3))/255;

        if(R==1 && G==1 && B==1)
            continue;
        end

        x(k) = ((j-1)/(n-1))*2-1;
        y(k) = ((i-1)/(m-1))*2-1;

        k = k+1;
    end
end

if(k<m*n)
    x(k:end) = [];
    y(k:end) = [];
end

pos = get(gcf,'Position');
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
axlim = max(ylim(2)-ylim(1),xlim(2)-xlim(1));
figlim = max(pos(3:4));

%     xm0 = x*(axlim/20)/(figlim/560);
%     ym0 = y*(axlim/20)/(figlim/560);

xm0 = ms*x*(axlim/20)/(figlim/400);
ym0 = ms*y*(axlim/20)/(figlim/400);

end

%% LOCAL SUB-FUNCTIONS
function [imr] = resizeImage(imgfile, factor)

im = imread(imgfile);

[m, n, ~] = size(im);
imr = im;

%     dm = round((m-1)/(m*factor-1));
%     dn = round((n-1)/(n*factor-1));

dm = round(1/factor-1);
dn = round(1/factor-1);

for j = 2:(n-1)/dn+1
    [~, n1, ~] = size(imr);
    if(j<=n1)
        if(j+dn-1<=n1)
            imr(:,j:j+dn-1,:) = [];
        else
            imr(:,j:end,:) = [];
        end
    else
        break;
    end
end

for i = 2:(m-1)/dm+1
    [m1, ~ , ~] = size(imr);
    if(i<=m1)
        if(i+dm-1<=m1)
            imr(i:i+dm-1,:,:) = [];
        else
            imr(i:end-1,:,:) = [];
        end
    else
        break;
    end
end
% imshow(imr)
end


function [xm, ym] = calcMarkerCoord(xm0, ym0, ang)

rotM = [cos(ang) -sin(ang);sin(ang) cos(ang)];

P = rotM*[xm0';ym0'];
xm = P(1,:)';
ym = P(2,:)';

end




