function [resultarr,outomega] = autorunner(n)
    
    % simulating time
    limtime = 5;
    alltime = 0;
    delta_t = 1e-5;
    
    % result setting
    calculatingomega = true;
    calculateomega_after = 2;
    plodisplayseperate = 1000;
    
    % result runtime parameter
    totalomega = 0;
    currentframe = 0;
    resultarr = [0,0];
    currentindicator = 1;
    
    % container porperties
    centiapply = true;
    container_radius = 0.05;
    container_radiusofswing = 0.1;
    container_frequency = [0,0,0.5];
    
    radius = 0.007;
    
    % container runtime parameter
    totalrotation = 0;
    container_omega = 2*pi*container_frequency(3);
    
    % render setting
    isrender = true;
    isrenderline = true;
    isrendercircle = true;
    isrendervector = false;
    isrendercenterofmass = true;
    renderhertz = 50;
    timefordisplay = 1;
    
    rendermarker = false;
    marker_radius = 0.002;
    
    centerofmassmarker_radius = 0.002;
    
    % render runtime parameter
    timeforrenderpass = 0;
    centerofmass = [0,0,0];
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
        coriolisvector = gobjects(n);
        centifugalvector = gobjects(n);
        velocityvector = gobjects(n);
        
        for circle = 1:n % skip edges
	        diameter = 2*radius;
	        temp = ballarr(circle).position(1,1:2) - radius;
        
	        position = [temp(1) temp(2) diameter diameter];
            graphics(circle) = rectangle('Position', position, 'Curvature', [1 1], 'EdgeColor', 'black', 'FaceColor', [0.75 + rand()*0.25 0.0 0.0]); % circles are drawn as squares with curved corners
        end
        
        if rendermarker == true
            position = [container_radius - marker_radius,-marker_radius,marker_radius,marker_radius];
            ballmarker = rectangle('Position', position, 'Curvature', [1 1], 'EdgeColor', 'black', 'FaceColor', [0.0 0.0 0.0]);
        end
        
        centerofmass = [0,0,0];
        for i = 1:n
            centerofmass = centerofmass + ballarr(i).position;
        end
        centerofmass = centerofmass/n;
        position = [centerofmass(1) - centerofmassmarker_radius,centerofmass(2) - centerofmassmarker_radius,centerofmassmarker_radius,centerofmassmarker_radius];
        centerofmass_marker = rectangle('Position', position, 'Curvature', [1 1], 'EdgeColor', 'black', 'FaceColor', [0.0 0.0 0.0]);
    
    
    end
    
    for i = 1:n
    
        containerswing_direc = (ballarr(i).position) - container_radiusofswing * [cos(totalrotation),sin(totalrotation),0];
        containerswing_direc(3) = 0;
        
        container_velocity = cross(2*pi*container_frequency,containerswing_direc);
    
    end
    
    while alltime < limtime
    
        forcevector = zeros(n,3);
        torqueforcevector = zeros(n,3);
        coriolisvectorarr = zeros(n,3);
        centifugalvectorarr = zeros(n,3);
        
        % O(n^2) solution; naive compare
        for idx = 1:n
            
            containerswing_direc = ballarr(idx).position - container_radiusofswing * [0,1,0];
            
            container_velocity = cross(2*pi*container_frequency,containerswing_direc);
            
            isincontactwitlwall = false;

            if magnitude(ballarr(idx).position) + radius > container_radius
                
                isincontactwitlwall = true;
                container_contact_velocity = cross(2*pi*container_frequency,containerswing_direc/magnitude(containerswing_direc) * container_radius);
    
            end
    
            [Fi,fi] = contactforce_wall(ballarr(idx),container_radius,container_velocity,isincontactwitlwall,delta_t);
            
            centifugalforce = ballarr(idx).mass * container_omega^2 * containerswing_direc * 100;
            centifugalvectorarr(idx,:) = centifugalforce(:);
    
            Fi = Fi + centifugalforce;
            %ballarr(idx).velocity = ballarr(idx).velocity + container_velocity;
    
            coriolisforce = 2 * ballarr(idx).mass * cross(2*pi*container_frequency,ballarr(idx).velocity);
            coriolisforce(3) = 0;
            coriolisvectorarr(idx,:) = coriolisforce(:);
         
            Fi = Fi + coriolisforce;
    
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
            temp = calculate_omega(ballarr,n,delta_t);

            if currentframe == plodisplayseperate
                resultarr(currentindicator) = temp;
                currentframe = 0;
                currentindicator = currentindicator + 1;
            else
                currentframe = currentframe + 1;
            end

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
    
            if isrendervector == true
        
                delete(coriolisvector);
                delete(centifugalvector);
                delete(velocityvector);
                coriolisvector = gobjects(n);
                centifugalvector = gobjects(n);
                velocityvector = gobjects(n);
        
                for circle = 1:n 
                    %radius = ballarr(circle).radius;
    
                    velomag = magnitude(ballarr(idx).velocity);
                    position = ballarr(circle).position;
    
                    coriolisvector(circle) = line(position(1) + [0, coriolisvectorarr(circle,1)], position(2) + [0,coriolisvectorarr(circle,2)], 'color', 'green');
                    centifugalvector(circle) = line(position(1) + [0, centifugalvectorarr(circle,1)], position(2) + [0,centifugalvectorarr(circle,2)], 'color', 'blue');
                    velocityvector(circle) = line(position(1) + [0, ballarr(circle).velocity(1)/velomag * container_radius/5], position(2) + [0,ballarr(circle).velocity(2)/velomag * container_radius/5], 'color', 'black');
        
                end
        
            end
    
            if rendermarker == true
                pos = (container_radius * [cos(totalrotation),sin(totalrotation),0,0]) + [-marker_radius,-marker_radius,marker_radius,marker_radius];
                set(ballmarker,'Position',pos)
            end
            
            if isrendercenterofmass == true
                
                centerofmass = 0;
                for i = 1:n
                    centerofmass = centerofmass + ballarr(i).position;
                end
                centerofmass = centerofmass/n;
                pos = [centerofmass(1) - centerofmassmarker_radius,centerofmass(2) - centerofmassmarker_radius,centerofmassmarker_radius,centerofmassmarker_radius];
                set(centerofmass_marker,'Position',pos)
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
    
    %plot(resultarr);

    outomega = totalomega/ (2*pi*pauseiteration);
        
end

