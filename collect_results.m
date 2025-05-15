function collect_results(file, run, coord, world, topo, output_dir, timesteps)

    % this function collects results from the PIV analysis and saves them to a specified output directory.

    tic_start = tic;
    fprintf('--------------------------------------------------\n');
    fprintf('Running collect_results.m \n');
    fprintf('--------------------------------------------------\n');

    javaaddpath('HydrolabPIV-master/src/measures');
    javaaddpath('HydrolabPIV-master/src/interp');
    addpath(genpath('HydrolabPIV-master/'));
    
    load(coord, 'pixel');

    [tform, ~, ~] = createcoordsystem(pixel, world, 'cubic');

    %MODEL
    load(topo, 'transformedTriangulation');
    Points = [transformedTriangulation.Points];


    xGrid = linspace(0, 1080, 1000);
    yGrid = linspace(0, 1080, 1000);
    [xq,yq] = meshgrid(xGrid, yGrid);

    
    

    zOriginal = Points(:,3);
    zMin = min(zOriginal);
    zMax = max(zOriginal);
    zNormalized = 5.5 * (zOriginal - zMin) / (zMax - zMin);
    Points(:, 3) = zNormalized;
    F_topo = scatteredInterpolant(Points(:, 1), Points(:, 2), Points(:, 3), 'natural');
    zq = F_topo(xq, yq);

    [dzx, dzy] = gradient(zq, xGrid, yGrid);
    [UU,VV,xx,yy] = pixel2world(tform, dzx, dzy, xq, yq, 1);

    F_dzx = scatteredInterpolant(xx(:), yy(:), UU(:), 'natural');
    F_dzy = scatteredInterpolant(xx(:), yy(:), VV(:), 'natural');

    U = cell(size(timesteps));
    V = cell(size(timesteps));
    x = cell(size(timesteps));
    y = cell(size(timesteps));
    cosTheta_masked_all = cell(size(timesteps));
    interpolated_dzx_all = cell(size(timesteps));
    interpolated_dzy_all = cell(size(timesteps));
    normInterpolated_dzx_all = cell(size(timesteps));
    normInterpolated_dzy_all = cell(size(timesteps));
    normU_all = cell(size(timesteps));
    normV_all = cell(size(timesteps));
    valid_mask_all = cell(size(timesteps));
    masked_normInterpolated_dzx_all = cell(size(timesteps));
    masked_normInterpolated_dzy_all = cell(size(timesteps));
    masked_normU_all = cell(size(timesteps));
    masked_normV_all = cell(size(timesteps));
    x_masked_all = cell(size(timesteps));
    y_masked_all = cell(size(timesteps));

    load(file, 'U2', 'V2', 'x2', 'y2');

    for timestep = timesteps


        if mod(timestep, 100) == 0
            fprintf('  Processing timestep %d of %d (%.1f%%)\n', timestep, timesteps(end), timestep/length(timesteps)*100);
        end

        U_t = U2{timestep};
        V_t = V2{timestep};
        x_t = x2;
        y_t = y2;

        [U1_t, V1_t, x1_t, y1_t] = pixel2world(tform, U_t, V_t, x_t, y_t, 1/5);

        U_current = U1_t;
        V_current = V1_t;
        x_current = x1_t;
        y_current = y1_t;

        interpolated_dzx = F_dzx(x_current, y_current);
        interpolated_dzy = F_dzy(x_current, y_current);


        normInterpolated_dzx = interpolated_dzx ./ sqrt(interpolated_dzx.^2 + interpolated_dzy.^2);
        normInterpolated_dzy = interpolated_dzy ./ sqrt(interpolated_dzx.^2 + interpolated_dzy.^2);


        normU = U_current ./ sqrt(U_current.^2 + V_current.^2);
        normV = V_current ./ sqrt(U_current.^2 + V_current.^2);


        threshold = 1e-4;
        valid_mask = (sqrt(interpolated_dzx.^2 + interpolated_dzy.^2) > threshold);


        masked_normInterpolated_dzx = normInterpolated_dzx(valid_mask);
        masked_normInterpolated_dzy = normInterpolated_dzy(valid_mask);
        masked_normU = normU(valid_mask);
        masked_normV = normV(valid_mask);

        cosTheta_masked = -masked_normInterpolated_dzy .* masked_normU + masked_normInterpolated_dzx .* masked_normV;
        cosTheta_masked = abs(cosTheta_masked);
        

        x_masked = x_current(valid_mask);
        y_masked = y_current(valid_mask);

        if timestep == 1
            x_min = min(x_current(:));
            x_max = max(x_current(:));
            y_min = min(y_current(:));
            y_max = max(y_current(:));
            numPoints = 1000;
            xContourGrid = linspace(x_min, x_max, numPoints);
            yContourGrid = linspace(y_min, y_max, numPoints);
            [xContour, yContour] = meshgrid(xContourGrid, yContourGrid);

            F_world = scatteredInterpolant(xx(:), yy(:), zq(:), 'natural');
            z_world = F_world(xContour, yContour);
        end

        cosTheta_masked_all{timestep} = cosTheta_masked;
        interpolated_dzx_all{timestep} = interpolated_dzx;
        interpolated_dzy_all{timestep} = interpolated_dzy;
        normInterpolated_dzx_all{timestep} = normInterpolated_dzx;
        normInterpolated_dzy_all{timestep} = normInterpolated_dzy;
        normU_all{timestep} = normU;
        normV_all{timestep} = normV;
        valid_mask_all{timestep} = valid_mask;
        masked_normInterpolated_dzx_all{timestep} = masked_normInterpolated_dzx;
        masked_normInterpolated_dzy_all{timestep} = masked_normInterpolated_dzy;
        masked_normU_all{timestep} = masked_normU;
        masked_normV_all{timestep} = masked_normV;
        x_masked_all{timestep} = x_masked;
        y_masked_all{timestep} = y_masked;
        U{timestep} = U_current;
        V{timestep} = V_current;
        x{timestep} = x_current;
        y{timestep} = y_current;
    end

    save(output_dir, 'timesteps', 'U', 'V', 'x', 'y', 'interpolated_dzx_all', 'interpolated_dzy_all',...
    'normInterpolated_dzx_all', 'normInterpolated_dzy_all', 'normU_all', 'normV_all', 'valid_mask_all', ...
    'masked_normInterpolated_dzx_all', 'masked_normInterpolated_dzy_all', 'masked_normU_all', 'masked_normV_all',...
    'cosTheta_masked_all', 'x_masked_all', 'y_masked_all',...
    'xContour','yContour','z_world');
    

    total_time = toc(tic_start);
    minutes = floor(total_time / 60);
    seconds = total_time - minutes * 60;
 
    fprintf('--------------------------------------------------\n');
    fprintf('EXECUTION TIME SUMMARY for Run %d\n', run);
    fprintf('--------------------------------------------------\n');
    fprintf('Total time:      %d min %.1f sec (%.2f sec)\n', minutes, seconds, total_time);
    fprintf('Timesteps:       %d\n', length(timesteps));
    fprintf('Time per step:   %.2f sec\n', total_time/length(timesteps));
    fprintf('--------------------------------------------------\n');
end
