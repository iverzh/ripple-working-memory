function b = bar_striped(x, data, barWidth, stripeSpacing, theta)

% Define stripe angle in degrees (e.g., 45 for 45? stripes)
m = tan(deg2rad(theta));  % slope of the lines

b = bar(x, data, barWidth);
hold on;


% Loop over each bar (assumes a single bar series)
for i = 1:length(b.XData)
    % Determine the x and y limits for this bar.
    x_center = b.XData(i);
    x0 = x_center - barWidth/2;
    x1 = x_center + barWidth/2;
    y0 = 0;
    y1 = data(i);
    
    % Compute candidate c values for each corner of the bar.
    % For a line: y = m*x + c, the intercept is c = y - m*x.
    c_candidates = [y0 - m*x0, y1 - m*x0, y0 - m*x1, y1 - m*x1];
    c_min = min(c_candidates);
    c_max = max(c_candidates);
    
    
    % Use a while loop to ensure the first stripe is plotted
    c_val = c_min;
    while c_val <= c_max
        pts = [];  % to store intersection points
        
        % Intersection with left vertical edge (x = x0)
        y_left = m*x0 + c_val;
        if y_left >= y0 && y_left <= y1
            pts = [pts; x0, y_left];
        end
        
        % Intersection with right vertical edge (x = x1)
        y_right = m*x1 + c_val;
        if y_right >= y0 && y_right <= y1
            pts = [pts; x1, y_right];
        end
        
        % Intersection with bottom horizontal edge (y = y0)
        if m ~= 0
            x_bottom = (y0 - c_val) / m;
            if x_bottom >= x0 && x_bottom <= x1
                pts = [pts; x_bottom, y0];
            end
        end
        
        % Intersection with top horizontal edge (y = y1)
        if m ~= 0
            x_top = (y1 - c_val) / m;
            if x_top >= x0 && x_top <= x1
                pts = [pts; x_top, y1];
            end
        end
        
        % Draw the stripe if exactly 2 intersection points are found
        if size(pts,1) == 2
            plot(pts(:,1), pts(:,2), 'k-');  % plot in black (adjust color as needed)
        end
        
        c_val = c_val + stripeSpacing;
    end
end

hold off;

