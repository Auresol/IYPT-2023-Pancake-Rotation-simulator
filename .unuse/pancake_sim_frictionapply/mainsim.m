clc
clear all

% simulating time
limtime = 20;
alltime = 0;
delta_t = 1e-5;

% result setting
calculatingomega = true;
calculateomega_after = 10;

% result runtime parameter
totalomega = 0;

% container porperties
centiapply = true;
n = 10;
container_radius = 5;
container_radiusofswing = 4
container_frequency = [0,0,1];

radius = 0.7;

% container runtime parameter
totalrotation = 0;
container_omega = 2*pi*container_frequency(3);

% render setting
isrender = true;
isrenderline = true;
isrendercircle = true;
renderhertz = 10;
timefordisplay = 1;

% render runtime parameter
timeforrenderpass = 0;
pauseiteration = 0;

%ball_contactforce_constant = 4/3 * sqrt(radius/2) * young_modulus/(2*(1- poisson_cof * poisson_cof));

ballarr = ballgeneration(n,radius);

if isrender == true

    % set up GUI
    fig = figure('Name', 'Example 1: Ball stack', 'Units', 'normalized', 'Position', [0.2 0.2 0.6 0.6]);
    
    % axes
    set(gca, 'OuterPosition', [0.1 0.2 0.8 0.8])
    axis([-container_radius,container_radius,-container_radius,container_radius])
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

end



while alltime < limtime

    forcevector = zeros(n,3);
    torqueforcevector = zeros(n,3);
    
    % O(n^2) solution; naive compare
    for idx = 1:n
        
        containerswing_direc = (ballarr(idx).position) - container_radiusofswing * [cos(totalrotation),sin(totalrotation),0];
        containerswing_direc(3) = 0;
        
        container_velocity = cross(2*pi*container_frequency,containerswing_direc);

        totalrotation = totalrotation + container_omega * delta_t;

        if totalrotation > 2*pi
            totalrotation = totalrotation - 2*pi;

        end
        
        isincontactwitlwall = false;

        if magnitude(ballarr(idx).position) + radius > container_radius
            
            isincontactwitlwall = true;

        end

        [Fi,fi] = contactforce_wall(ballarr(idx),container_radius,container_velocity,isincontactwitlwall,delta_t);

        %coriolisforce = 2 * ballarr(idx).mass * magnitude(velo) * container_omega;
        %coriolis_direc = cross(velo/magnitude(velo),[0,0,1]);
        %coriolis_direc(3) = 0;
        
        %if velo ~= 0
            %Fi = Fi + coriolis_direc * coriolisforce;
        %end

        % apply force
        forcevector(idx,:) = forcevector(idx,:) + Fi;
        
        % friction angular momentum calculation
        torqueforcevector(idx,:) = torqueforcevector(idx,:) + fi;

        
            
        for idy = (idx+1):n
            
            % distance size
            dis = magnitude(ballarr(idx).position - ballarr(idy).position);

            if dis < 2*radius

                [Fi,Fj,fi,fj] = contactforce_ball(ballarr(idx),ballarr(idy),delta_t);

                % apply force
                forcevector(idx,:) = forcevector(idx,:) + Fi;
                forcevector(idy,:) = forcevector(idy,:) + Fj;
                
                % friction angular momentum calculation
                torqueforcevector(idx,:) = torqueforcevector(idx,:) + fi;
                torqueforcevector(idy,:) = torqueforcevector(idy,:) + fj;
                
            end
        end
    end
     
    for iter = 1:n
        ballarr(iter) = applyforce(ballarr(iter),forcevector(iter,:),torqueforcevector(iter,:),delta_t);
    end
    
    alltime = alltime + delta_t;
    timeforrenderpass = timeforrenderpass + delta_t;
  
    if timeforrenderpass > timefordisplay
        timeforrenderpass = timeforrenderpass - timefordisplay;
        disp("TIME PASS " + alltime);
    end

    if calculatingomega == true && alltime > calculateomega_after
        totalomega = totalomega + calculate_omega(ballarr,n,delta_t);
    end

    if pauseiteration == renderhertz && isrender == true
        pauseiteration = 0;
        % ball render
        for circle = 1:n % skip edges
            %radius = ballarr(circle).radius;
            
            diameter = 2 * radius;
	        temp = ballarr(circle).position(1,1:2) - radius;
        
	        position = [temp(1) temp(2) diameter diameter];
            
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

disp(totalomega/ pauseiteration)
disp(totalomega/ (2*pi*pauseiteration))
    