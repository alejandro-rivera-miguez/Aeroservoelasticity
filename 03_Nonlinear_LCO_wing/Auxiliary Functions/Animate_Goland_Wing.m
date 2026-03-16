%% =======================================================================
% FUNCTION: Animate_Goland_Wing
% DESCRIPTION: Generates a 3D video animation (.avi) of the Goland wing 
%              undergoing Limit Cycle Oscillations (LCO). It projects the 
%              FEM modal shapes into a 3D surface grid.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Alberto Rivero García
%% =======================================================================
function Animate_Goland_Wing(qr, qF, omega, video_filename)
    
    fprintf('--- STARTING VIDEO GENERATION ---\n');
    
    %% 1. VALIDATION AND PREPARATION
    % Ensure omega is real and positive
    omega = abs(real(omega));
    if omega < 1e-3
        omega = 1.0; 
        fprintf('WARNING: Omega is too low or null. Using omega=1 for visualization.\n');
    end
    
    % Goland geometric data
    L = 6.096; c = 1.829; xe = 0.33 * c;
    yL = 3.35; yU = 3.35 + 1.8; cf = 0.45; xf = c - cf;
    
    % Total deformation calculation
    Q_total = qr * qF; 
    real_scale = 50; % Scaling factor for visual exaggeration
    
    %% 2. MESH GENERATION
    x_vec = unique(sort([linspace(0, c, 25), xf])); 
    y_vec = unique(sort([linspace(0, L, 50), yL, yU]));
    [X_mesh, Y_mesh] = meshgrid(x_vec, y_vec);
    
    tol = 1e-5;
    is_flap_region = (Y_mesh >= yL-tol) & (Y_mesh <= yU+tol) & (X_mesh >= xf-tol);
    
    % FEM Interpolation
    Ne = 10;
    y_fem_bend = linspace(0, L, Ne + 1);
    y_fem_tors = linspace(0, L, 2*Ne + 1);
    w_fem = Q_total(1:2:21);
    phi_fem = Q_total(23:43);
    beta_val = Q_total(44);
    
    %% 3. GRAPHICS CONFIGURATION
    % 'Visible', 'on' ensures the frame is captured correctly by VideoWriter
    f = figure('Color','w', 'Units', 'pixels', 'Position', [100 100 1280 720], 'Visible', 'on');
    
    % --- Main Wing ---
    Z_wing_init = zeros(size(X_mesh)); 
    Z_wing_init(is_flap_region) = NaN; 
    
    hWing = surf(X_mesh, Y_mesh, Z_wing_init, ...
        'FaceColor', [0.15 0.15 0.15], ...  % Dark gray
        'EdgeColor', 'none', ...
        'FaceLighting', 'gouraud', ...
        'AmbientStrength', 0.6, ...         
        'SpecularStrength', 0.9);           % For 3D lighting reflections
    hold on;
    
    % --- Aileron (Flap) ---
    Z_flap_init = zeros(size(X_mesh)); Z_flap_init(~is_flap_region) = NaN;
    hFlap = surf(X_mesh, Y_mesh, Z_flap_init, 'FaceColor', [0.8 0.2 0.2], ...
        'EdgeColor', 'none', 'FaceLighting', 'gouraud', 'AmbientStrength', 0.5);
        
    % --- Hinge Line ---
    hHinge = plot3([xf xf], [yL yU], [0 0], 'b-', 'LineWidth', 2);
    
    axis equal; grid on; box on; view(130, 25);
    z_limit = 2.5; 
    zlim([-z_limit z_limit]); xlim([0 c]); ylim([0 L]);
    xlabel('Chord [m]'); ylabel('Wingspan [m]'); zlabel('Z [m]');
    title(sprintf('LCO Simulation @ %.2f rad/s', omega));
    light('Position', [-5 -5 10], 'Style', 'local');
    light('Position', [5 5 5], 'Style', 'local');
    
    %% 4. VIDEO GENERATION (TRY-CATCH BLOCK)
    % Motion JPEG AVI is widely compatible
    vidObj = VideoWriter(video_filename, 'Motion JPEG AVI');
    vidObj.FrameRate = 30; 
    vidObj.Quality = 95;
    
    try
        open(vidObj);
        
        T = 2*pi/omega; 
        num_frames = 150;
        t_vec = linspace(0, 4*T, num_frames); 
        
        fprintf('Generating %d frames...\n', num_frames);
        
        for i = 1:num_frames
            osc = exp(1i * omega * t_vec(i));
            
            % Instantaneous deformation
            w_inst = real(w_fem * osc) * real_scale;
            phi_inst = real(phi_fem * osc) * real_scale;
            beta_inst = real(beta_val * osc) * real_scale;
            
            W_grid = interp1(y_fem_bend, w_inst, Y_mesh, 'spline');
            Phi_grid = interp1(y_fem_tors, phi_inst, Y_mesh, 'spline');
            
            Z_base = W_grid - (X_mesh - xe) .* Phi_grid;
            
            % Update Wing graphics
            Z_w = Z_base; Z_w(is_flap_region) = NaN;
            set(hWing, 'ZData', Z_w);
            
            % Update Flap graphics
            Z_f = Z_base;
            Z_f(is_flap_region) = Z_f(is_flap_region) - (X_mesh(is_flap_region) - xf) * beta_inst;
            Z_f(~is_flap_region) = NaN;
            set(hFlap, 'ZData', Z_f);
            
            % Update Hinge
            y_h = linspace(yL, yU, 10);
            w_h = interp1(y_fem_bend, w_inst, y_h, 'spline');
            phi_h = interp1(y_fem_tors, phi_inst, y_h, 'spline');
            z_h = w_h - (xf - xe) .* phi_h;
            set(hHinge, 'YData', y_h, 'ZData', z_h);
            
            % Force rendering and capture frame
            drawnow limitrate; 
            writeVideo(vidObj, getframe(f));
            
            if mod(i, 30) == 0
                fprintf('  -> Frame %d / %d processed.\n', i, num_frames);
            end
        end
        
        fprintf('--- SUCCESS: Video saved as "%s" ---\n', video_filename);
        
    catch ME
        fprintf('\n!!! CRITICAL ERROR DURING VIDEO GENERATION !!!\n');
        fprintf('Message: %s\n', ME.message);
        % Close video object safely to avoid file corruption
        close(vidObj);
        close(f);
        rethrow(ME);
    end
    
    % Clean up
    close(vidObj);
    close(f);
end