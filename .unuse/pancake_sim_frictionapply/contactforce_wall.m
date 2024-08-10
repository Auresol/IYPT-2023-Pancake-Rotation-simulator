function [Fi,fi] = contactforce_wall(ball1,container_radius,container_velocity,incontactwithwall,delta_t)
    
    % --------------------------- %
    wall_young_modulus = 1e9;
    wall_poisson_cof = 0.25;
    wall_hertz_daming_cof = 0.4;
    
    % friction force parameter
    friction_compliment = 1e15;
    static_mue_w = 1;
    kinetic_mue_w = 1;
    % ---------------------------- %
    
    Fi = [0,0,0];
    fi = [0,0,0];

    % direction vector (x -> wall)
    direcvec_ground = [0,0,1];

    % overlap parameter
    overlap_ground = -(ball1.position(3));
    if overlap_ground < 0
        overlap_ground = 0;
    end

    %overlap_ground = 0.01;
    
    % R_effective
    radius_effective = ball1.radius;

    % E_effective
    young_modulus_effective = 1/( (1 - ball1.poisson_cof^2)/ball1.young_modulus + (1 - wall_poisson_cof ^ 2)/wall_young_modulus );

    % contact area parameter
    contactarea_ground = pi * radius_effective/2 * overlap_ground;
    
    contact_const = 4/3 * sqrt(radius_effective/2) * young_modulus_effective;

    % contact force (scalar)
    contactforce_ground = contact_const * sqrt(overlap_ground)^3;
    
    % vRel = -vPi
    vRel_ground = container_velocity - (ball1.velocity + cross(ball1.angular_velocity, ball1.radius * -direcvec_ground));

    vRelt_ground = vRel_ground - (dot(vRel_ground, direcvec_ground)) * direcvec_ground;
    

    if magnitude(vRelt_ground) ~= 0
        
        % friction force calculation
        % friction torque direction (x -> y)
        tDirection_ground = vRelt_ground/magnitude(vRelt_ground);
        
        friction_condition = friction_compliment * magnitude(vRelt_ground) * contactarea_ground * delta_t;
        
        frictionforce_ground = 0;
        if friction_condition < static_mue_w * contactforce_ground
            frictionforce_ground = friction_condition;
        else
            frictionforce_ground = kinetic_mue_w * contactforce_ground;
        end 

        fi = fi + cross(-direcvec_ground * ball1.radius, frictionforce_ground * tDirection_ground);
        Fi = Fi + frictionforce_ground * tDirection_ground;

    end

    % damping calculation
    % relative speed at surface
    vReln_ground = dot(vRel_ground,direcvec_ground) * direcvec_ground;

     % mass_effective
    mass_effective = ball1.mass;

    % damping_coefficient_effective
    damping_effective = (ball1.young_modulus * ball1.damping_cof + wall_young_modulus * wall_hertz_daming_cof)/(ball1.young_modulus + wall_young_modulus);
    
    % damping force
    dampingforce_ground = 2*damping_effective * sqrt(2 * young_modulus_effective * mass_effective) * (radius_effective * overlap_ground)^(1/4);
    
    % apply force
    Fi = Fi + contactforce_ground * direcvec_ground + dampingforce_ground * vReln_ground;

    if incontactwithwall == true

        % direction vector (x -> wall)
        direcvec = -ball1.position/magnitude(ball1.position);
    
        % overlap parameter
        overlap = magnitude(ball1.position) + ball1.radius - container_radius;
    
        % contact area parameter
        contactarea = pi * radius_effective/2 * overlap;
    
        % contact force (scalar)
        contactforce =  contact_const * sqrt(overlap^3);
        
        % vRel = -vPi
        vRel = -(ball1.velocity + cross(ball1.angular_velocity, ball1.radius * -direcvec));
    
        vRelt = vRel - (dot(vRel, direcvec)) * direcvec;
        
        if magnitude(vRelt) ~= 0
            
            % friction force calculation
            % friction torque direction (x -> y)
            tDirection = vRelt/magnitude(vRelt);
            
            friction_condition = friction_compliment * magnitude(vRelt) * contactarea * delta_t;
            
            frictionforce = 0;
            if friction_condition < static_mue_w * contactforce
                frictionforce = friction_condition;
            else
                frictionforce = kinetic_mue_w * contactforce;
            end 
    
            fi = fi + cross(-direcvec * ball1.radius, frictionforce * tDirection);
            Fi = Fi + frictionforce * tDirection;
    
        end
    
        % damping calculation
        % relative speed at surface
        vReln = dot(vRel,direcvec) * direcvec;
        
        % damping force
        dampingforce = 2*damping_effective * sqrt(2 * young_modulus_effective * mass_effective) * (radius_effective * overlap)^(1/4);
        
        % apply force
        Fi = Fi + contactforce * direcvec + dampingforce * vReln;
    end

end

