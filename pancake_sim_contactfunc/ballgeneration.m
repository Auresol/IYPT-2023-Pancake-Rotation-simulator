function ballarr = ballgeneration(n,radius)
    
    % ball data
    centershift = [0,0,0];
    mass = 5;

    initial_position = closepacking_initial(n);
    
    ballarr = ball(initial_position(1) + centershift,radius,mass);

    for i = 2:n
        ballarr(i) = ball(initial_position(i,:)*2*radius + centershift,radius,mass);
    end
    
end

