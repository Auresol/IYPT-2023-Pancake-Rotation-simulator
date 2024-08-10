function totalomega = calculate_omega(ballarr,n,delta_t)
    %CALCULATE_OMEGA Summary of this function goes here
    %   Detailed explanation goes here
    totalomega = 0;
    center_of_mass = [0,0,0];

    for i = 1:n
        center_of_mass = center_of_mass + ballarr(i).position;
    end

    center_of_mass = center_of_mass/n;

    for i = 1:n
        pastvec = center_of_mass - ballarr(i).past_position;
        curvec = center_of_mass - ballarr(i).position;

        angle = acos(dot(pastvec,curvec)/(magnitude(pastvec) * magnitude(curvec)));

        omega = angle/delta_t;
        totalomega = totalomega + omega;

    end

    totalomega = totalomega/n;



end

