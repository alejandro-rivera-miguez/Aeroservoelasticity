%% =======================================================================
% FUNCTION: w_interpolate
% DESCRIPTION: Weighted interpolation function defined over the Puma blade 
%              grid. It ensures that the total integral of the interpolated 
%              quantity (mass, stiffness, etc.) is conserved between the 
%              original and the new discrete meshes.
%
% INPUTS:
%   in_grid  - Original spatial grid array
%   in_data  - Original data array corresponding to in_grid
%   out_grid - Target spatial grid array for interpolation
%
% OUTPUTS:
%   out_data - Interpolated data conserving the integral quantities
%   out_g2   - Expanded target grid used for piecewise constant assignment
%% =======================================================================
function [out_data, out_g2] = w_interpolate(in_grid, in_data, out_grid)
   
    % Size the arrays 
    % The output array out_g2 is the expanded grid for piecewise constant data
    out_g2 = zeros(2 * (length(out_grid) - 1), 1); 
    out_g2([1:2:end, end]) = out_grid;
    out_g2(2:2:end-1)      = out_grid(2:end-1);
    
    out_data = zeros(length(out_g2), 1); 
    
    %% First Element Processing
    % Nodes of the input grid that are larger than or equal to the final node of the first output element 
    idx_1 = find(in_grid >= out_grid(2));
    endn = in_grid(idx_1(1));
    
    % Nodes of the input grid that are smaller than or equal to the first node of the output element
    idx_2 = find(in_grid <= out_grid(1));
    inn = in_grid(idx_2(end));
    
    if (idx_2(end) + 1 == idx_1(1))
        % The output element is fully contained in the input element
        delta = (in_data(idx_1(1)) - in_data(idx_2(end))) / (endn - inn);
        out_data(1) = (2 * in_data(idx_2(end)) + delta * (out_grid(1) - inn) + delta * (out_grid(2) - inn)) / 2;
        out_data(2) = out_data(1);
    else
        % Fractional area integration
        dt1 = trapz(in_grid(idx_2(end):idx_1(1)-1), in_data(idx_2(end):idx_1(1)-1));
        scale2 = endn - in_grid(idx_1(1)-1);
        dt2 = trapz(in_grid(idx_1(1)-1:idx_1(1)), in_data(idx_1(1)-1:idx_1(1)));
        
        out_data(1) = (dt1 + dt2 * (out_grid(2) - in_grid(idx_1(1)-1)) / scale2) / (out_grid(2) - out_grid(1));
        out_data(2) = out_data(1);
    end
   
    %% Intermediate Elements Processing
    for i = 3:length(out_grid)-1
        idx_1 = find(in_grid >= out_grid(i));
        endn = in_grid(idx_1(1));
        
        idx_2 = find(in_grid <= out_grid(i-1));
        inn = in_grid(idx_2(end));
        
        if (idx_2(end) + 1 == idx_1(1))
            % The output element is fully contained in the input element
            delta = (in_data(idx_1(1)) - in_data(idx_2(end))) / (endn - inn);
            out_data((i-1)*2-1) = (2 * in_data(idx_2(end)) + delta * (out_grid(i-1) - inn) + delta * (out_grid(i) - inn)) / 2;
            out_data((i-1)*2) = out_data((i-1)*2-1);
        else    
            % Fractional area integration over multiple input sub-elements
            scale0 = in_grid(idx_2(end)+1) - inn; 
            dt0 = trapz(in_grid(idx_2(end):idx_2(end)+1), in_data(idx_2(end):idx_2(end)+1));
            
            dt1 = trapz(in_grid(idx_2(end)+1:idx_1(1)-1), in_data(idx_2(end)+1:idx_1(1)-1));
            
            scale2 = endn - in_grid(idx_1(1)-1); 
            dt2 = trapz(in_grid(idx_1(1)-1:idx_1(1)), in_data(idx_1(1)-1:idx_1(1)));
            
            out_data((i-1)*2-1) = (dt0 * (in_grid(idx_2(end)+1) - out_grid(i-1)) / scale0 + dt1 + dt2 * (out_grid(i) - in_grid(idx_1(1)-1)) / scale2) / (out_grid(i) - out_grid(i-1));
            out_data((i-1)*2) = out_data((i-1)*2-1);
        end
    end
   
    %% Last Element Processing
    idx_1 = find(in_grid >= out_grid(end));
    endn = in_grid(idx_1(1));
    
    idx_2 = find(in_grid <= out_grid(end-1));
    inn = in_grid(idx_2(end));
    
    if (idx_2(end) + 1 == idx_1(1))
        % The output element is fully contained in the input element
        delta = (in_data(idx_1(1)) - in_data(idx_2(end))) / (endn - inn);
        out_data(end-1) = (2 * in_data(idx_2(end)) + delta * (out_grid(end-1) - inn) + delta * (out_grid(end) - inn)) / 2;
        out_data(end) = out_data(end-1);
    else
        % Fractional area integration
        scale0 = in_grid(idx_2(end)+1) - inn;
        dt0 = trapz(in_grid(idx_2(end):idx_2(end)+1), in_data(idx_2(end):idx_2(end)+1));
        
        dt1 = trapz(in_grid(idx_2(end)+1:idx_1(1)), in_data(idx_2(end)+1:idx_1(1)));
        
        out_data(end-1) = (dt0 * (in_grid(idx_2(end)+1) - out_grid(end-1)) / scale0 + dt1) / (out_grid(end) - out_grid(end-1));
        out_data(end) = out_data(end-1);
    end
end