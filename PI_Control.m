function [Vc_mc,Vg_mc_track,Vg_mc_track_0,V_load,V_mc_load_0,PI,DFB]=PI_Control(PI,Vrf_ideal,Vg_mc_track_0,...
    V_mc_load_0,V_load_cpu,TbAng_coef_mc,pattern,DFB)
% PI_CONTROL Realistic PI feedback and optional Direct RF feedback control
%
% Usage:
%   % Legacy PI-only mode (100% backward compatible):
%   [Vc_mc,Vg_mc_track,Vg_mc_track_0,V_load,V_mc_load_0,PI] = ...
%       PI_Control(PI,Vrf_ideal,Vg_mc_track_0,V_mc_load_0,V_load_cpu,TbAng_coef_mc,pattern);
%
%   % Combined PI + Direct RF Feedback mode:
%   [Vc_mc,Vg_mc_track,Vg_mc_track_0,V_load,V_mc_load_0,PI,DFB] = ...
%       PI_Control(PI,Vrf_ideal,Vg_mc_track_0,V_mc_load_0,V_load_cpu,TbAng_coef_mc,pattern,DFB);

if nargin < 8
    DFB = [];
end

[Vc_mc, Vg_mc_track, Vg_mc_track_0, V_load, V_mc_load_0, PI, DFB] = MC_Cavity_Control(...
    PI, DFB, Vrf_ideal, Vg_mc_track_0, V_mc_load_0, V_load_cpu, TbAng_coef_mc, pattern);

end
