clc
clear all

limtime = 10;
alltime = 0;
delta_t = 1e-4;

% container porperties
centiapply = true;
n = 20;
container_radius = 10;
container_radiusofswing = 1500;
omega_c = 0;
container_frequency = 1;

radius = 1;

% container runtime parameter
totalrotation = 0;

% physics constants
worldSize = [10 10]; % [m]
eulermethod_scale = 0.5;
F_adh = [0,0,0];
F_env = 6*pi * obj.air_viscosity * obj.radius * (obj.local_average_velocity_of_external_air_at_contact - obj.velocity);
F_g = [0,0,0];
%F_g = obj.mass * obj.gravity_acc;

totalextraforce = F_adh + F_env + F_g;

% render setting
isrenderline = true;
isrendercircle = true;
renderhertz = 2;

% runder runtime parameter
pauseiteration = 0;

% ball_contactforce_constant = 4/3 * sqrt(radius/2) * young_modulus/(2*(1- poisson_cof * poisson_cof));

ballarr = ballgeneration(n);

% set up GUI
fig = figure('Name', 'Example 1: Ball stack', 'Units', 'normalized', 'Position', [0.2 0.2 0.6 0.6]);

% axes
set(gca, 'OuterPosition', [0.1 0.2 0.8 0.8])
axis([-10 worldSize(1) -10 worldSize(2)])
axis square
axis off

% outer render
diameter = container_radius*2;
outerline = rectangle('Position',[-container_radius,-container_radius,diameter,diameter],'Curvature',[1,1],'EdgeColor', 'black');

%balltemp = line([ballarr(1).position(1),ballarr(1).position(1)] , [-ballarr(1).position(2),ballarr(1).position(2)], 'color', 'black');

% circles
graphics = gobjects(n);
linegraphic = gobjects(n);

for circle = 1:n % skip edges
	diameter = 2*radius;
	temp = ballarr(circle).position(1,1:2) - radius;

	position = [temp(1) temp(2) diameter diameter];
    graphics(circle) = rectangle('Position', position, 'Curvature', [1 1], 'EdgeColor', 'black', 'FaceColor', [0.75 + rand()*0.25 0.0 0.0]); % circles are drawn as squares with curved corners
end



while alltime < limtime
    
    % FL
    forcevector_FL = zeros(n,3) + totalextraforce;
    % ML
    torqueforcevector_FL = zeros(n,3);

    % O(n^2) solution; naive compare
    for idx = 1:n

        if centiapply == true
        
            % centifugal force calculation
            % centifugal normal direction
            centifugal_direcvec = ball_iteration_arr(idx).position - (container_radiusofswing * [cos(totalrotation),sin(totalrotation),0]);
            centifugal_direcvec = centifugal_direcvec/magnitude(centifugal_direcvec);
    
            totalrotation = totalrotation + 2 * pi * delta_t * container_frequency;
            
            if totalrotation > 2 * pi
                totalrotation = totalrotation - 2*pi;
            end

            centifugalforce = ball_iteration_arr(idx).mass * (2 * pi * container_frequency)^2 * container_radiusofswing;
            
            % apply centifugal
            forcevector_FL(idx,:) = forcevector_FL(idx,:) + centifugalforce * centifugal_direcvec;

        end

        if magnitude(ball_iteration_arr(idx).position) + radius > container_radius
            
            [Fi,fi] = contactforce_wall(ball_iteration_arr(idx),container_radius,delta_t);

            % apply force
            forcevector_FL(idx,:) = forcevector_FL(idx,:) + Fi;
            
            % friction angular momentum calculation
            torqueforcevector_FL(idx,:) = torqueforcevector_FL(idx,:) + fi;

        end
            
        for idy = (idx+1):n
            
            % distance size
            dis = magnitude(ball_iteration_arr(idx).position - ball_iteration_arr(idy).position);

            if dis < 2*radius

                [Fi,Fj,fi,fj] = contactforce_ball(ball_iteration_arr(idx),ball_iteration_arr(idy),delta_t);

                % apply force
                forcevector_FL(idx,:) = forcevector_FL(idx,:) + Fi;
                forcevector_FL(idy,:) = forcevector_FL(idy,:) + Fj;
                
                % friction angular momentum calculation
                torqueforcevector_FL(idx,:) = torqueforcevector_FL(idx,:) + fi;
                torqueforcevector_FL(idy,:) = torqueforcevector_FL(idy,:) + fj;
                
            end
        end
    end

    iteration = 0;

    ball_iteration_arr(:) = ballarr(:);

    if eulermethod_scale == 0
        continue
    end
    
    % iteration
    centirotated_temp = totalrotation + 2 * pi * delta_t * container_frequency;

    while

        forcevector = zeros(n,3) + totalextraforce;
        torqueforcevector = zeros(n,3);
        
        % O(n^2) solution; naive compare
        for idx = 1:n
    
            if centiapply == true
            
                % centifugal force calculation
                % centifugal normal direction
                centifugal_direcvec = ball_iteration_arr(idx).position - (container_radiusofswing * [cos(centirotated_temp),sin(centirotated_temp),0]);
                centifugal_direcvec = centifugal_direcvec/magnitude(centifugal_direcvec);

                centifugalforce = ball_iteration_arr(idx).mass * (2 * pi * container_frequency)^2 * container_radiusofswing;
                
                % apply centifugal
                forcevector(idx,:) = forcevector(idx,:) + centifugalforce * centifugal_direcvec;
    
            end
    
            if magnitude(ball_iteration_arr(idx).position) + radius > container_radius
                
                [Fi,fi] = contactforce_wall(ball_iteration_arr(idx),container_radius,delta_t);
    
                % apply force
                forcevector(idx,:) = forcevector(idx,:) + Fi;
                
                % friction angular momentum calculation
                torqueforcevector(idx,:) = torqueforcevector(idx,:) + fi;
    
            end
                
            for idy = (idx+1):n
                
                % distance size
                dis = magnitude(ball_iteration_arr(idx).position - ball_iteration_arr(idy).position);
    
                if dis < 2*radius
    
                    [Fi,Fj,fi,fj] = contactforce_ball(ball_iteration_arr(idx),ball_iteration_arr(idy),delta_t);
    
                    % apply force
                    forcevector(idx,:) = forcevector(idx,:) + Fi;
                    forcevector(idy,:) = forcevector(idy,:) + Fj;
                    
                    % friction angular momentum calculation
                    torqueforcevector(idx,:) = torqueforcevector(idx,:) + fi;
                    torqueforcevector(idy,:) = torqueforcevector(idy,:) + fj;
                    
                end
            end
        end
        
        ball_iteration_arr(:) = ballarr(:);
        
        
        for iter = 1:n
            ball_iteration_arr(iter) = applyforce(ball_iteration_arr(iter),forcevector_FL(iter,:),forcevector(iter,:), ...
                torqueforcevector_FL(iter,:),torqueforcevector(iter,:),delta_t);
        end

    end
    
    alltime = alltime + delta_t;

    if pauseiteration == renderhertz
        pauseiteration = 0;
        % ball render
        for circle = 1:n % skip edges
            %radius = ballarr(circle).radius;
            
            diameter = 2 * radius;
	        temp = ballarr(circle).position(1,1:2) - radius;
        
	        position = [temp(1) temp(2) 2 2];
            
            set(graphics(circle), 'Position', position);
        end
        
        % line render
        if isrenderline == true
    
            delete(linegraphic);
            linegraphic = gobjects(n);
    
            for circle = 1:n 
                %radius = ballarr(circle).radius;

                rotated = ballarr(circle).rotation;
                position = ballarr(circle).position;
                
                xconst = radius*cos(rotated(3));
                yconst = radius*sin(rotated(3));
                linegraphic(circle) = line(position(1) + [-xconst, xconst], position(2) + [-yconst,yconst], 'color', 'black');
    
            end
    
        end
        pause(0.01);

    else
        pauseiteration = pauseiteration + 1;

    end
       
    %delete(balltemp)
    %balltemp = line([ballarr(1).position(1) - radius,ballarr(1).position(1) + radius] , [ballarr(1).position(2),ballarr(1).position(2)], 'color', 'black');
    
    
    
end
    