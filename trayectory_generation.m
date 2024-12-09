function [waypoints] = ...
trayectory_generation(waypoints,type_int,type_yaw,r,t)

%% WAYPOINTS COORDINATES, YAW AND SPEED
 x = waypoints(:,1);
 y = waypoints(:,2);
 z = waypoints(:,3);
 yaw = waypoints(:,4);
 v = waypoints(:,5);

 % Curvilinear coordinate distance in horizontal plane
 s = ...
 cumsum(sqrt(sum((waypoints(2:end,1:2)-waypoints(1:end-1,1:2)).^2,2)));

 % If interpolation function is lineal: Intermediate lineal points are
 % calculated and saved in vector sq: Curvilinear distances
 % If 1: lineal interpolation, If other: curved interpolation
 if type_int == 1
 npoints = 100;
 sq = zeros(length(s),npoints);
 lows = [0;s+r];
 upps = [s-r;0];
 for i=1:length(s)
 sq(i,:) = linspace(lows(i),upps(i),npoints);
 end
 sq = reshape(sq',1,[])';
 else
 sq = [0;s];
 end
 s = [0;s];

 % Coordinates of intermediate points
 xq = interp1(s,x,sq,'linear');
 yq = interp1(s,y,sq,'linear');
 zq = interp1(s,z,sq,'linear');
 yawq = interp1(s,yaw,sq,'linear');
 vq = interp1(s,v,sq,'linear');

 % Coordinates of points in the trajectory
 %sq2 = 0:0.01:s(end);%====== ORIGINAL =======
 sq2 = 0:0.001:s(end);
 xq2 = interp1(sq,xq,sq2,'spline');
 yq2 = interp1(sq,yq,sq2,'spline');
 zq2 = interp1(sq,zq,sq2,'spline');
 yawq2 = interp1(sq,yawq,sq2,'spline');
 vq2 = interp1(sq,vq,sq2,'spline');
 % Check type of yaw: If 1 yaw parallel to velocity
 if type_yaw==1
 yawq2 = atan2(yq2(2:end)-yq2(1:end-1),xq2(2:end)-xq2(1:end-1));
 yawq2 = [yawq2,yawq2(end)];
 end

 % References as function of time considering speed
 s_check = 0;
 waypoints = zeros(length(t),5);
 for i = 1:length(t)-1
 [~,ind] = min(abs(sq2-s_check));

 waypoints(i+1,:) = [xq2(ind),yq2(ind),zq2(ind),yawq2(ind),vq2(ind)];
 s_check = s_check+vq2(ind)*(t(i+1)-t(i));
 end
 waypoints(1,:) = waypoints(2,:);

 end
