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

        radius = 0;
        mass = 0.01;
        inertia = 1;

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
        
        function obj = applyforce(obj,forcevec,rotated_momentum,delta_t)
            
            F_adh = [0,0,0];
            F_env = 6*pi * obj.air_viscosity * obj.radius * (obj.local_average_velocity_of_external_air_at_contact - obj.velocity);
            %F_g = [0,0,0];
            F_g = obj.mass * obj.gravity_acc;

            forcevec = forcevec + F_adh + F_env + F_g;
            
            obj = obj.pastAssign();
            obj.position = obj.past_position + obj.past_velocity * delta_t ;
            obj.velocity = obj.past_velocity + forcevec * delta_t / obj.mass;
            obj.angular_velocity = obj.past_angular_velocity + rotated_momentum * delta_t / obj.inertia;
            obj.rotation = obj.past_rotation + obj.past_angular_velocity * delta_t;
            
        end

        function obj = pastAssign(obj)

            obj.past_velocity = obj.velocity;
            obj.past_position = obj.position;
            obj.past_angular_velocity = obj.angular_velocity;
            obj.past_rotation = obj.rotation;

        end

    end
end

