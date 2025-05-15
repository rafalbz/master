% Aligns the model with the reference image
function simple_model_alignment()
    % Load STL model
    fv = stlread('C:/Users/Rafal/Documents/Master/piv_source/3_1_Hill.stl');
    swappedPoints = fv.Points(:, [1, 3, 2]); % Swap Y and Z axes
    
    % Filter out triangles below threshold
    zThreshold = -1;
    keepTriangles = false(size(fv.ConnectivityList, 1), 1);
    for i = 1:size(fv.ConnectivityList, 1)
        triangleVertices = swappedPoints(fv.ConnectivityList(i, :), :);
        if all(triangleVertices(:, 3) >= zThreshold)
            keepTriangles(i) = true;
        end
    end
    
    % Get filtered connectivity list and points
    filteredTriangles = fv.ConnectivityList(keepTriangles, :);
    
    if isempty(filteredTriangles)
        filteredTriangles = fv.ConnectivityList;
        originalSwappedPoints = swappedPoints;
    else
        uniqueVertices = unique(filteredTriangles(:));
        originalSwappedPoints = swappedPoints(uniqueVertices, :);
        [~, filteredTriangles] = ismember(filteredTriangles, uniqueVertices);
    end
    
    % Load reference image
    img = imread('C:\Users\Rafal\Documents\Master\piv_source\final\28.03\scaling\ridge\ridge2.png');
    [imgHeight, imgWidth, ~] = size(img);
    
    % Create figure and UI
    fig = figure('Position', [100, 100, 1200, 700], 'Name', 'Model Alignment');
    ax = axes('Position', [0.05, 0.1, 0.65, 0.85]);
    imshow(img);
    hold on;
    
    % Initial values and slider ranges
    scaleXInit = 4.91;
    scaleYInit = 5.66;
    scaleZInit = 5;
    rotZInit = 91.7;
    transXInit = 597;
    transYInit = -1485.8;
    
    % Slider ranges
    scaleXRange = [1, 15];
    scaleYRange = [1, 15];
    scaleZRange = [1, 50];
    rotZRange = [0, 180];
    transXRange = [-1000, 1000];
    transYRange = [-6000, 3000];
    
    % Calculate slider steps
    function steps = calculateSliderStep(minVal, maxVal, increment)
        range = maxVal - minVal;
        smallStep = increment / range;
        largeStep = increment * 5 / range;
        steps = [smallStep, largeStep];
    end
    
    % UI Controls Panel
    uipanel('Position', [0.75, 0.1, 0.22, 0.85], 'Title', 'Control Parameters');
    

    baseY = 650;
    spacing = 50;
    

    uicontrol('Style', 'text', 'Position', [950, baseY, 120, 20], 'String', 'X Scale:');
    hScaleX = uicontrol('Style', 'slider', 'Position', [950, baseY-20, 200, 20], ...
        'Min', scaleXRange(1), 'Max', scaleXRange(2), 'Value', scaleXInit, ...
        'SliderStep', calculateSliderStep(scaleXRange(1), scaleXRange(2), 0.05), ...
        'Callback', @updateVisualization);
    hScaleXText = uicontrol('Style', 'text', 'Position', [1160, baseY-20, 50, 20], ...
        'String', num2str(scaleXInit, '%.2f'));
    

    baseY = baseY - spacing;
    uicontrol('Style', 'text', 'Position', [950, baseY, 120, 20], 'String', 'Y Scale:');
    hScaleY = uicontrol('Style', 'slider', 'Position', [950, baseY-20, 200, 20], ...
        'Min', scaleYRange(1), 'Max', scaleYRange(2), 'Value', scaleYInit, ...
        'SliderStep', calculateSliderStep(scaleYRange(1), scaleYRange(2), 0.05), ...
        'Callback', @updateVisualization);
    hScaleYText = uicontrol('Style', 'text', 'Position', [1160, baseY-20, 50, 20], ...
        'String', num2str(scaleYInit, '%.2f'));
    

    baseY = baseY - spacing;
    uicontrol('Style', 'text', 'Position', [950, baseY, 120, 20], 'String', 'Z Scale:');
    hScaleZ = uicontrol('Style', 'slider', 'Position', [950, baseY-20, 200, 20], ...
        'Min', scaleZRange(1), 'Max', scaleZRange(2), 'Value', scaleZInit, ...
        'SliderStep', calculateSliderStep(scaleZRange(1), scaleZRange(2), 0.05), ...
        'Callback', @updateVisualization);
    hScaleZText = uicontrol('Style', 'text', 'Position', [1160, baseY-20, 50, 20], ...
        'String', num2str(scaleZInit, '%.2f'));
    

    baseY = baseY - spacing;
    uicontrol('Style', 'text', 'Position', [950, baseY, 120, 20], 'String', 'Rotation Z:');
    hRotZ = uicontrol('Style', 'slider', 'Position', [950, baseY-20, 200, 20], ...
        'Min', rotZRange(1), 'Max', rotZRange(2), 'Value', rotZInit, ...
        'SliderStep', calculateSliderStep(rotZRange(1), rotZRange(2), 0.1), ...
        'Callback', @updateVisualization);
    hRotZText = uicontrol('Style', 'text', 'Position', [1160, baseY-20, 50, 20], ...
        'String', num2str(rotZInit, '%.2f'));
    

    baseY = baseY - spacing;
    uicontrol('Style', 'text', 'Position', [950, baseY, 120, 20], 'String', 'Translation X:');
    hTransX = uicontrol('Style', 'slider', 'Position', [950, baseY-20, 200, 20], ...
        'Min', transXRange(1), 'Max', transXRange(2), 'Value', transXInit, ...
        'SliderStep', calculateSliderStep(transXRange(1), transXRange(2), 1.0), ...
        'Callback', @updateVisualization);
    hTransXText = uicontrol('Style', 'text', 'Position', [1160, baseY-20, 50, 20], ...
        'String', num2str(transXInit, '%.1f'));
    

    baseY = baseY - spacing;
    uicontrol('Style', 'text', 'Position', [950, baseY, 120, 20], 'String', 'Translation Y:');
    hTransY = uicontrol('Style', 'slider', 'Position', [950, baseY-20, 200, 20], ...
        'Min', transYRange(1), 'Max', transYRange(2), 'Value', transYInit, ...
        'SliderStep', calculateSliderStep(transYRange(1), transYRange(2), 1.0), ...
        'Callback', @updateVisualization);
    hTransYText = uicontrol('Style', 'text', 'Position', [1160, baseY-20, 50, 20], ...
        'String', num2str(transYInit, '%.1f'));
    

    baseY = baseY - spacing;
    hLinkScales = uicontrol('Style', 'checkbox', 'Position', [950, baseY-20, 200, 20], ...
        'String', 'Link X, Y, and Z Scales', 'Value', 1, 'Callback', @linkScalesCallback);
    

    baseY = baseY - spacing;
    uicontrol('Style', 'pushbutton', 'Position', [950, baseY-20, 200, 30], ...
        'String', 'Save Triangulation', 'Callback', @saveTriangulation);
    

    baseY = baseY - spacing;
    hStatus = uicontrol('Style', 'text', 'Position', [950, baseY-20, 250, 30], ...
        'String', 'Ready', 'HorizontalAlignment', 'left');
    

    hModel = [];
    
    % Initial visualization
    updateVisualization([]);
    
    % Link scales function
    function linkScalesCallback(~, ~)
        if get(hLinkScales, 'Value') == 1
            set(hScaleY, 'Value', get(hScaleX, 'Value'));
            set(hScaleYText, 'String', get(hScaleXText, 'String'));
            set(hScaleZ, 'Value', get(hScaleX, 'Value'));
            set(hScaleZText, 'String', get(hScaleXText, 'String'));
            updateVisualization([]);
        end
    end
    
    % Save triangulation function
    function saveTriangulation(~, ~)
        % Get current parameter values
        scaleX = get(hScaleX, 'Value');
        scaleY = get(hScaleY, 'Value');
        scaleZ = get(hScaleZ, 'Value');
        rotZ = get(hRotZ, 'Value');
        transX = get(hTransX, 'Value');
        transY = get(hTransY, 'Value');
        
        % Apply transformations
        scaledPoints = originalSwappedPoints;
        scaledPoints(:, 1) = originalSwappedPoints(:, 1) * scaleX;
        scaledPoints(:, 2) = originalSwappedPoints(:, 2) * scaleY;
        scaledPoints(:, 3) = originalSwappedPoints(:, 3) * scaleZ;
        
        % Rotate around Z only
        rotZRad = deg2rad(rotZ);
        Rz = [cos(rotZRad), -sin(rotZRad), 0;
              sin(rotZRad), cos(rotZRad), 0;
              0, 0, 1];
        
        % Apply rotation
        rotatedPoints = (Rz * scaledPoints')';
        
        % Apply floor flattening
        originalMinZ = min(scaledPoints(:, 3));
        originalPointsAtMinZ = abs(scaledPoints(:, 3) - originalMinZ) < 1e-6;
        minZAfterRotation = min(rotatedPoints(originalPointsAtMinZ, 3));
        rotatedPoints(:, 3) = rotatedPoints(:, 3) - minZAfterRotation;
        
        % Apply translation (X and Y only)
        transformedPoints = rotatedPoints;
        transformedPoints(:, 1) = transformedPoints(:, 1) + transX;
        transformedPoints(:, 2) = transformedPoints(:, 2) + transY;
        
        % Create the triangulation object
        transformedTriangulation = triangulation(filteredTriangles, transformedPoints);
        
        % Save to MAT file with timestamp
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        filename = ['transformed_triangulation_', timestamp, '.mat'];
        save(filename, 'transformedTriangulation');
        
        % Update status
        set(hStatus, 'String', ['Saved to ', filename]);
    end
    
    % Update visualization function
    function updateVisualization(src, ~)
        % Get current parameter values
        scaleX = get(hScaleX, 'Value');
        scaleY = get(hScaleY, 'Value');
        scaleZ = get(hScaleZ, 'Value');
        rotZ = get(hRotZ, 'Value');
        transX = get(hTransX, 'Value');
        transY = get(hTransY, 'Value');
        
        % Handle linked scales if toggled
        if ~isempty(src)
            if isequal(src, hScaleX) && get(hLinkScales, 'Value') == 1
                set(hScaleY, 'Value', scaleX);
                set(hScaleYText, 'String', num2str(scaleX, '%.2f'));
                set(hScaleZ, 'Value', scaleX);
                set(hScaleZText, 'String', num2str(scaleX, '%.2f'));
                scaleY = scaleX;
                scaleZ = scaleX;
            elseif isequal(src, hScaleY) && get(hLinkScales, 'Value') == 1
                set(hScaleX, 'Value', scaleY);
                set(hScaleXText, 'String', num2str(scaleY, '%.2f'));
                set(hScaleZ, 'Value', scaleY);
                set(hScaleZText, 'String', num2str(scaleY, '%.2f'));
                scaleX = scaleY;
                scaleZ = scaleY;
            elseif isequal(src, hScaleZ) && get(hLinkScales, 'Value') == 1
                set(hScaleX, 'Value', scaleZ);
                set(hScaleXText, 'String', num2str(scaleZ, '%.2f'));
                set(hScaleY, 'Value', scaleZ);
                set(hScaleYText, 'String', num2str(scaleZ, '%.2f'));
                scaleX = scaleZ;
                scaleY = scaleZ;
            end
        end
        
        % Update text displays
        set(hScaleXText, 'String', num2str(scaleX, '%.2f'));
        set(hScaleYText, 'String', num2str(scaleY, '%.2f'));
        set(hScaleZText, 'String', num2str(scaleZ, '%.2f'));
        set(hRotZText, 'String', num2str(rotZ, '%.2f'));
        set(hTransXText, 'String', num2str(transX, '%.1f'));
        set(hTransYText, 'String', num2str(transY, '%.1f'));
        
        % Apply transformations
        scaledPoints = originalSwappedPoints;
        scaledPoints(:, 1) = originalSwappedPoints(:, 1) * scaleX;
        scaledPoints(:, 2) = originalSwappedPoints(:, 2) * scaleY;
        scaledPoints(:, 3) = originalSwappedPoints(:, 3) * scaleZ;
        
        % Rotate around Z only
        rotZRad = deg2rad(rotZ);
        Rz = [cos(rotZRad), -sin(rotZRad), 0;
              sin(rotZRad), cos(rotZRad), 0;
              0, 0, 1];
        
        % Apply rotation
        rotatedPoints = (Rz * scaledPoints')';
        
        % Apply floor flattening
        originalMinZ = min(scaledPoints(:, 3));
        originalPointsAtMinZ = abs(scaledPoints(:, 3) - originalMinZ) < 1e-6;
        minZAfterRotation = min(rotatedPoints(originalPointsAtMinZ, 3));
        rotatedPoints(:, 3) = rotatedPoints(:, 3) - minZAfterRotation;
        
        % Apply translation (X and Y only)
        transformedPoints = rotatedPoints;
        transformedPoints(:, 1) = transformedPoints(:, 1) + transX;
        transformedPoints(:, 2) = transformedPoints(:, 2) + transY;
        
        % Update visualization
        axes(ax);
        
        % Remove previous model if it exists
        if ~isempty(hModel) && isvalid(hModel)
            delete(hModel);
        end
        
        % Display the new model
        combinedTriangulation = triangulation(filteredTriangles, transformedPoints);
        hModel = trisurf(combinedTriangulation, 'FaceColor', 'none', 'EdgeColor', 'y', 'LineWidth', 0.4);
        
        % Set axes properties
        axis equal;
        set(gca, 'XLim', [0, imgWidth], 'YLim', [0, imgHeight]);
        
        % Update status
        set(hStatus, 'String', 'Model updated');
        drawnow;
    end
end