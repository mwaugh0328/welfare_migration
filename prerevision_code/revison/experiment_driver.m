function [sim_expr, sim_cntr] = experiment_driver(assets, move, assets_temp,...
    move_temp, cons_eqiv, params, trans_shocks, monga_index,...
    states, pref_shocks, sim_data)%#codegen

Nsims = length(monga_index);

sim_expr = zeros(Nsims,11,11);
% sim_surv = zeros(Nsims,10,2);
sim_cntr = zeros(Nsims,9,11);

% This will implement the experiment in the monga season.

for xxx = 1:Nsims
    
    index = monga_index(xxx);
    
    state_at_expr = states(index,1:3);
    
    shock_states = states(index:index+10,4);
    
    p_shocks = pref_shocks(index:index+10,1);
       
    [panel_expr] = simmulate_experiment(assets, move, assets_temp, move_temp, cons_eqiv, params,...
    state_at_expr, trans_shocks, shock_states, p_shocks);
            
    sim_expr(xxx,:,1) = panel_expr(1,:);
    sim_expr(xxx,:,2) = panel_expr(2,:);
    sim_expr(xxx,:,3) = panel_expr(3,:);
    sim_expr(xxx,:,4) = panel_expr(4,:);
    sim_expr(xxx,:,5) = panel_expr(5,:);
    sim_expr(xxx,:,6) = panel_expr(6,:);
    sim_expr(xxx,:,7) = panel_expr(7,:);
    sim_expr(xxx,:,8) = panel_expr(8,:);
    sim_expr(xxx,:,9) = panel_expr(9,:);
    sim_expr(xxx,:,10) = panel_expr(10,:);
    sim_expr(xxx,:,11) = panel_expr(11,:);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sim_cntr(xxx,:,1) = sim_data(index,:); 
    sim_cntr(xxx,:,2) = sim_data(index+1,:); 
    sim_cntr(xxx,:,3) = sim_data(index+2,:); 
    sim_cntr(xxx,:,4) = sim_data(index+3,:); 
    sim_cntr(xxx,:,5) = sim_data(index+4,:); 
    sim_cntr(xxx,:,6) = sim_data(index+5,:); 
    sim_cntr(xxx,:,7) = sim_data(index+6,:); 
    sim_cntr(xxx,:,8) = sim_data(index+7,:); 
    sim_cntr(xxx,:,9) = sim_data(index+8,:); 
    sim_cntr(xxx,:,10) = sim_data(index+9,:);
    sim_cntr(xxx,:,11) = sim_data(index+10,:);

    % So as setup, this is just taking the choices when the experiment is
    % conducted. 
 
end

% if sum(sim_cntr(:,end-1)) ~= length(sim_expr(:,end-1))
%     disp('something is wrong')
% test = 1
% end
    
end

