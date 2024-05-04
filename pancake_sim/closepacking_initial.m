function output = closepacking_initial(n)
%CLOSEPACKING_INITIAL Summary of this function goes here
%   Detailed explanation goes here
    n = n+10;
    cycle = 0;
    items = 2;

    z = 0;

    output = zeros(n,3);

    while items < n
        cycle = cycle + 1;
        
        pos1 = [0,cycle,z];
        
        for k = 1:6
            output(items,:) = pos1(1,:);
            items = items + 1;

            if items > n
                break
            end

            pos2 = pos1;
            pos1 = [cycle * sin(k * pi/3), cycle * cos(k * pi/3),z];

            for inbet = 1:cycle-1
                t1 = (pos1 * inbet / cycle + pos2 * (cycle-inbet)/cycle);
                
                output(items,:) = t1(1,:);
                items = items + 1;

                if items > n
                    break
                end

            end

            if items > n
                break
            end
        end
    end

end

