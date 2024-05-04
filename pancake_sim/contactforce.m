function [Fi,Fj] = contactforce(ball1,ball2)
    % direction vector (x -> wall)
            direcvec = -ballarr(idx).position/magnitude(ballarr(idx).position);

            % overlap parameter
            overlap = magnitude(ballarr(idx).position) + radius - container_radius;
            
            % young modulus effectiveness
            young_modulus_effectiveness = 1/(((1-poisson_cof^2)/young_modulus + (1-wall_poisson_cof^2)/wall_young_modulus));
            
            % contact force
            coctactforce = 4/3 * sqrt(radius) * young_modulus_effectiveness * sqrt(overlap^3) * direcvec;

            contactarea = pi * R/2 * overlap;
                
            vRel = -(ballarr(idx).velocity + cross(ballarr(idx).angular_velocity, radius * direcvec));

            vRelt = vPi - (dot(vPi, direcvec)) * -direcvec;  

            tDirection = vRelt/magnitude(vRelt);

            % friction torque direction (x -> y)
            frictiontorquevec = vRelt/magnitude(vRelt);
            
            friction_condition = friction_compliment * magnitude(vRelt) * contactarea * delta_t;

            if friction_condition < static_mue_w * magnitude(contactforce)
                frictionforce = friction_condition;
            else
                frictionforce = kinetic_mue_w * contactforce;
            end 

            % damping calculation
            vReln = dot(vRel,direcvec) * direcvec;

            dampingforce = 2*hertz_daming_cof * sqrt(young_modulus/(1- poisson_cof * poisson_cof) * mass/2) * (radius/2 * overlap)^(1/4);
            
            % apply force
            forcevector(idx,:) = forcevector(idx,:) + contactforce * direcvec + dampingforce * vReln;

            torqueforcevector(idx,:) = torqueforcevector(idx,:) + frictionforce * cross(tDirection, (radius - overlap/2) * direcvec);

end

