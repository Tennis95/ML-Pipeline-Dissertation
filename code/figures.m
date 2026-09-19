% generate_extra_figures.m
% Generates: (1) flow field schematic, (2) spurious vs genuine detection example
% Run in MATLAB Online from /MATLAB Drive/dissertation-matlab/

cd('/MATLAB Drive/dissertation-matlab');

%% ── FIGURE 1: Flow Field Schematic ────────────────────────────────────────
% Draws a schematic of the mixing layer showing vortex roll-up
fig1 = figure('Position',[100 100 900 400], 'Color','white');
hold on; axis off;