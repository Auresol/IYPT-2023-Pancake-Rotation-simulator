classdef ball
    %BALL Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = private)

        % F_env calculation (external damping force)
        air_viscosity = 0.00265;
        local_average_velocity_of_external_air_at_contact = [0,0,0];

        % F_g calculation
        gravity_acc = [0,-9.81,0];

        past_velocity (1,3) = [0,0,0];
        past_position (1,3) = [0,0,0];
        past_angular_velocity (1,3) = [0,0,0];
        past_rotation (1,3) = [0,0,0];

    end
    
    properties (Access = public)

        id = 0;

        velocity (1,3) = [0,0,0];
        position (1,3) = [0,0,0];
        angular_velocity (1,3) = [0,0,0];
        rotation (1,3) = [0,0,0];

        radius (1,1) = 1;
        mass (1,1) = 0.01;
        inertia (1,1) = 0;

        young_modulus (1,1) = 1e9;
        poisson_cof (1,1) = 0.25;
        damping_cof (1,1) = 0.4;

    end
    
    methods

        function obj = ball(position,radius,mass)

            %BALL Construct an instance of this class
            %   Detailed explanation goes here
            obj.radius = radius;
            obj.mass = mass;
            obj.position = position;
            obj.inertia = 2/5 * mass * radius^2;

        end
        
        function obj = applyforce(obj,forcevec_FL,forcevec_FLL,rotated_momentum_FL,rotated_momentum_FLL,euler_scale,delta_t)
            
            obj = obj.pastAssign();
            obj.position = obj.past_position + obj.past_velocity * delta_t + euler_scale * delta_t^2/obj.mass * (euler_scale * forcevec_FLL + (1-euler_scale) * forcevec_FL);
            obj.velocity = obj.past_velocity + delta_t / obj.mass * (euler_scale * forcevec_FLL + (1-euler_scale) * forcevec_FL);
            obj.angular_velocity = obj.past_angular_velocity + delta_t / obj.inertia * (euler_scale * rotated_momentum_FLL + (1-euler_scale) * rotated_momentum_FL);
            obj.rotation = obj.past_rotation + obj.past_angular_velocity * delta_t;
            
        end

        function obj = pastAssign(obj)

            obj.past_velocity = obj.velocity;
            obj.past_position = obj.position;
            obj.past_angular_velocity = obj.angular_velocity;
            obj.past_rotation = obj.rotation;

        end

        function obj = past_iteration_assign(obj)

            obj.past_iter_velocity = obj.iter_velocity;
            obj.past_iter_position = obj.iter_position;
            obj.past_iter_angular_velocity = obj.iter_angular_velocity;

        end

    end
end

