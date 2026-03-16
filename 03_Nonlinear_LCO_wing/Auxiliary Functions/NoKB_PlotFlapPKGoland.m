%% =======================================================================
% FUNCTION: NoKB_PlotFlapPKGoland
% DESCRIPTION: Wrapper to plot flutter results sweeping up to Umax.
%
% AUTHORS: Francisco Javier Martin Lopez
%          Alejandro Rivera Miguez
%          Alberto Rivero García
% =======================================================================
function NoKB_PlotFlapPKGoland(Ne, Nd, Umax, Nu)
    % Directly calls the solver and enables the Display flag.
    [~, ~, ~] = NoKB_FlapPKGoland(Ne, Nd, Umax, Nu, 1);
end