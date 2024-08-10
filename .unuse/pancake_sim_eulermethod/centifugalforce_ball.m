function Fi = centifugalforce_ball(position,rotated,container_radius)

    if rotated > 2 * pi
        rotated = rotated - 2*pi;
    end

    % centifugal force calculation
    % centifugal normal direction
    centifugal_direcvec = position - (container_radiusofswing * [cos(rotated),sin(rotated),0]);
    centifugal_direcvec = centifugal_direcvec/magnitude(centifugal_direcvec);

    
end

