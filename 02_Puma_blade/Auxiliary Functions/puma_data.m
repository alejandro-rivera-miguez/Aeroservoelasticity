%% =======================================================================
% FUNCTION: puma_data
% DESCRIPTION: Contains the structural, geometric, and inertial properties 
%              of the Puma helicopter main rotor blade.
%
% OUTPUTS:
%   B       - Struct containing all distributed properties
%   datapos - Index array for parameters
%% =======================================================================
function [B, datapos] = puma_data()

    blade_chord = 0.537; % [m]
    
    %% Geometric and Kinematic Constants
    B.blade_radius    = 7.490;      % [m]
    B.blade_cutout    = 2.080;      % [m]
    B.blade_chord     = 0.537;      % [m]
    B.lag_hinge_axis  = 0.269;      % [m]
    B.flap_hinge_axis = 0.289;      % [m]
    B.pitch_link      = 0.289;      % [m]
    B.pitch_bearing   = 0.432;      % [m]
    B.omega           = 4.5 * 2*pi; % Rotational Speed [rad/s]
    B.initial_r       = 0.28;       % [m]
    
    eps = 1e-7; % Tolerance to avoid duplicate exact coordinates
    
    %% Distributed Properties Arrays: [Station (m), Value]
    
    % 1. Blade Twist [deg]
    B.blade_twist = [
        0.280,  0.000
        0.760,  0.000
        0.760,  0.000
        1.757,  0.000
        1.757,  0.000
        2.007, -0.100
        2.007, -0.100
        2.257, -0.300
        2.257, -0.300
        2.507, -0.617
        2.507, -0.617
        2.757, -0.983
        2.757, -0.983
        3.007, -1.283
        3.007, -1.283
        3.257, -1.600
        3.257, -1.600
        3.507, -1.867
        3.507, -1.867
        6.257, -4.800
        6.257, -4.800
        7.490, -6.117
    ];
    for i = 2:length(B.blade_twist)
        if abs(B.blade_twist(i, 1) - B.blade_twist(i-1, 1)) <= eps
            B.blade_twist(i, 1) = B.blade_twist(i, 1) + eps;
        end
    end
    
    % 2. Blade Mass [kg/m]
    B.blade_mass = [
        0.280, 58.400
        0.610, 58.400
        0.610, 50.000
        0.730, 50.000
        0.730, 16.175
        0.754, 16.175
        0.754, 53.333
        0.760, 53.333
        0.760, 24.949
        0.800, 24.949
        0.800, 33.100
        0.840, 33.100
        0.840, 22.400
        1.040, 22.400
        1.040, 17.110
        1.110, 17.110
        1.110, 11.225
        1.260, 11.225
        1.260,  7.150
        1.770,  7.150
        1.770,  7.150
        1.887,  8.929
        1.887,  8.929
        7.070,  8.929
        7.070, 12.754
        7.390, 12.754
        7.390, 34.600
        7.402, 34.600
        7.402,  6.930
        7.490,  6.930
    ];
    for i = 2:length(B.blade_mass)
        if abs(B.blade_mass(i, 1) - B.blade_mass(i-1, 1)) <= eps
            B.blade_mass(i, 1) = B.blade_mass(i, 1) + eps;
        end
    end
    
    % 3. Blade CG Offset (Chordwise) [non-dimensional y/c]
    B.blade_cg_offset = [
        0.280,  0.000
        1.887,  0.000
        1.887, -0.010
        7.070, -0.010
        7.070,  0.030
        7.390,  0.030
        7.390,  0.147
        7.402,  0.147
        7.402,  0.000
        7.490,  0.000
    ];
    for i = 2:length(B.blade_cg_offset)
        if abs(B.blade_cg_offset(i, 1) - B.blade_cg_offset(i-1, 1)) <= eps
            B.blade_cg_offset(i, 1) = B.blade_cg_offset(i, 1) + eps;
        end
    end
    
    % Convert non-dimensional y/c to actual meters
    B.blade_cg_offset(:, 2) = B.blade_cg_offset(:, 2) * blade_chord;
    
    % 4. Blade Torsional Inertia [kg*m]
    B.blade_torsional_inertia = [
        0.280, 0.000
        0.604, 0.000
        0.604, 0.192
        0.610, 0.192
        0.610, 0.164
        0.730, 0.178
        0.730, 0.116
        0.754, 0.130
        0.754, 0.427
        0.760, 0.438
        0.760, 0.205
        0.800, 0.240
        0.800, 0.318
        0.836, 0.359
        0.836, 0.359
        0.837, 0.120
        0.837, 0.120
        0.840, 0.121
        0.840, 0.082
        1.017, 0.121
        1.017, 0.121
        1.040, 0.098
        1.040, 0.075
        1.085, 0.040
        1.085, 0.040
        1.110, 0.040
        1.110, 0.026
        1.260, 0.026
        1.260, 0.017
        1.770, 0.037
        1.770, 0.087
        1.887, 0.109
        1.887, 0.109
        7.070, 0.109
        7.070, 0.156
        7.390, 0.156
        7.390, 0.337
        7.402, 0.067
        7.402, 0.067
        7.434, 0.067
        7.434, 0.118
        7.490, 0.118
    ];
    for i = 2:length(B.blade_torsional_inertia)
        if abs(B.blade_torsional_inertia(i, 1) - B.blade_torsional_inertia(i-1, 1)) <= eps
            B.blade_torsional_inertia(i, 1) = B.blade_torsional_inertia(i, 1) + eps;
        end
    end
    
    % 5. Blade Extensional Stiffness EA [N]
    B.blade_extensional_stiffness = [
        0.280, 5.69e8
        0.760, 5.69e8
        0.760, 5.63e8
        0.820, 5.63e8
        0.820, 3.77e8
        0.873, 3.77e8
        0.873, 5.69e8
        0.965, 5.69e8
        0.965, 5.50e8
        1.040, 5.50e8
        1.040, 4.74e8
        1.111, 4.74e8
        1.111, 3.74e8
        1.260, 3.74e8
        1.260, 1.70e8
        1.697, 1.69e8
        1.697, 1.69e8
        1.757, 1.65e8
        1.757, 1.65e8
        1.880, 1.71e8
        1.880, 1.71e8
        2.128, 1.70e8
        2.128, 1.70e8
        2.278, 1.51e8
        2.278, 1.51e8
        5.250, 1.45e8
        5.250, 1.45e8
        5.800, 1.44e8
        5.800, 1.44e8
        6.110, 1.43e8
        6.110, 1.43e8
        7.350, 1.43e8
        7.350, 1.61e8
        7.382, 1.61e8
        7.382, 1.16e8
        7.490, 1.16e8
    ];
    for i = 2:length(B.blade_extensional_stiffness)
        if abs(B.blade_extensional_stiffness(i, 1) - B.blade_extensional_stiffness(i-1, 1)) <= eps
            B.blade_extensional_stiffness(i, 1) = B.blade_extensional_stiffness(i, 1) + eps;
        end
    end
    
    % 6. Blade Flap/Chord Stiffness EJ [N*m^2] (Assuming identical values for hover study)
    B.blade_flap_chord_stiffness = [
        0.280, 178.00e4, 178.00e4
        0.600, 178.00e4, 178.00e4
        0.600, 178.00e4, 178.00e4
        0.610, 137.00e4, 137.00e4
        0.610, 137.00e4, 137.00e4
        0.800, 137.00e4, 137.00e4
        0.800, 137.00e4, 137.00e4
        0.810,  41.20e4, 178.00e4
        0.810,  41.20e4, 178.00e4
        1.240,  41.00e4, 178.00e4
        1.240,  41.00e4, 178.00e4
        1.250,   8.10e4, 153.00e4
        1.250,   8.10e4, 153.00e4
        7.300,   8.10e4, 144.00e4
        7.300,   8.10e4, 144.00e4
        7.310,   8.20e4,  71.50e4
        7.310,   8.20e4,  71.50e4
        7.490,   8.20e4,  71.50e4
    ];
    for i = 2:length(B.blade_flap_chord_stiffness)
        if abs(B.blade_flap_chord_stiffness(i, 1) - B.blade_flap_chord_stiffness(i-1, 1)) <= eps
            B.blade_flap_chord_stiffness(i, 1) = B.blade_flap_chord_stiffness(i, 1) + eps;
        end
    end
    
    % 7. Blade Torsional Stiffness GJ [N*m^2]
    B.blade_torsional_stiffness = [
        0.280,  84.00e4
        0.725,  84.00e4
        0.725, 226.00e4
        0.828, 226.00e4
        0.828,  50.50e4
        1.017,  50.50e4
        1.017,   8.50e4
        7.241,   8.50e4
        7.241,   8.70e4
        7.490,   8.70e4
    ];
    for i = 2:length(B.blade_torsional_stiffness)
        if abs(B.blade_torsional_stiffness(i, 1) - B.blade_torsional_stiffness(i-1, 1)) <= eps
            B.blade_torsional_stiffness(i, 1) = B.blade_torsional_stiffness(i, 1) + eps;
        end
    end
    
    %% Metadata Definitions
    datapos = 1:16;
    
    B.name = {
        'Station',
        'Twist',
        'Mass',
        'CG Position Y',
        'CG Position Z',
        'EA',
        'EJ Flap',
        'EJ Lag',
        'Neutral axis Y',
        'Neutral axis Z',
        'GJ',
        'Elastic axis Y',
        'Elastic axis Z',
        'Torsional inertia J',
        'alpha',
        'beta'
    };

    B.Units = {
        'm', 'deg', 'kg/m', 'y/c adim.', 'z/c adim.', 'N', 'N m2', 'N m2', ...
        'y/c adim.', 'z/c adim.', 'N m2', 'y/c adim.', 'z/c adim.', 'kg m', 'deg', 'deg'
    };

end