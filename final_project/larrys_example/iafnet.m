function [output_events, unit_potentials] = iafnet(input_events)
% output_events(clock, output_unit) -- Returns boolean output events for each output unit.
% output_potentials(clock, network_unit) -- Returns network unit potentials.
% input_events(clock, input_unit) -- Boolean event for each input at each time step.
%
% Integrate and Fire Network
% Written by Larry Shupe
% University of Washington
% Dept of Physiology and Biophysics
% August 2017
%
% An example of input_events where each input has a 1% chance at any clock
% to add PSPs to connected network units:
% >> input_events = random('bino', 1, 0.01, time_steps, input_units);
%
% This function runs a simulation of an integrate and fire neural network.
% Each simulated unit can receive spikes from three sources:
% 1) input_events: These are external input spikes from other systems.
%    For example, these might be spikes discriminated from recorded cells.
% 2) bias inputs: These are randomly generated spikes that will keep
%    a unit's potential above zero but mostly below threshold.
% 3) other units: Any network unit may be connected to any other.
%
% When a unit's potential reaches threshold, a PSP is delivered to
% all connected units at a short delay time, and the unit's potential
% is cleared.  Network units that are defined as outputs may be used
% to trigger a stimulus output.
%
% When a PSP arrives at a unit, the corresponding connection_strength
% is added into both the slow and fast decay potentials for that unit.
% A unit's potential = slow_potential(unit) - fast_potential(unit).
% At each timestep, the slow and fast potential values decay towards
% zero by multiplication by the slow and fast decay constants.
% Note that final PSP height will be lower that connection_strength
% because of the expontial decays.  For example, with slow decay = 31/32,
% and fast decay = 7/8, the PSP height is about half of the corresponding
% connection strength.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation runs for the number of time steps given in the input array.

[time_steps, input_units] = size(input_events);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup parameters
% Network units have indexes 1 to network_units.
% Output units have indexes 1 to output_units.
% Input units have indexes network_units+1 to total_units.

output_units = 6;       % Number of network units whose output spikes drive stimulators.
network_units = 8 + output_units; % Number of units for which we calculate unit potentials.
psp_slow_decay = 31/32; % Slow decay constant for PSPs.
psp_fast_decay = 7/8;   % Slow decay constant for PSPs.
bias_rate = 2000;       % Randomly distributed spikes per second to each network unit
bias_strength = 500;    % Amplitude of bias spike exponential functions (actual PSP height will be smaller)
unit_threshold = 5000;  % Unit potential threshold.
connection_delay = 30;  % Timestep spike conduction time, [0..31] due to implementation as a 32-bit FIFO queue.
sample_rate = 10000;    % 10000 time_steps per second.
total_units = network_units + input_units;
connection_strength = 500 * ones(total_units, network_units); % connection_strength(from_unit, to_unit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_events = zeros(time_steps, output_units);    % Allocate output arrays
unit_potentials = zeros(time_steps, network_units);

slow_potential = zeros(network_units, 1);
fast_potential = zeros(network_units, 1);
unit_spike_queue = uint32(zeros(network_units, 1)); % Use uint64 if longer connection_delay is needed.
delay_bit = uint32(2^connection_delay); % The bit to check in the unit_spike_queue() values.

for clock = 1:time_steps
    bitshift(unit_spike_queue, 1);  % All spike queues advance
    fast_potential = psp_fast_decay .* fast_potential;  % All unit potentials decay
    slow_potential = psp_slow_decay .* slow_potential;  % All unit potentials decay
    
    % Check for external input and spike delivery to network units.
    
    for unit = 1:total_units
        if unit > network_units
            % External input spikes go directly to network units.
            spike = input_events(clock, unit - network_units);
        else
            % Bias intput to the network units
            if random('uniform', 0, sample_rate) < bias_rate
                slow_potential(unit) = slow_potential(unit) + bias_strength;
                fast_potential(unit) = fast_potential(unit) + bias_strength;
            end
            % Check if a spike on this unit impacts other network units now.
            spike = bitand(unit_spike_queue(unit), delay_bit);
        end
        
        if spike
            % Time to deliver psp to all downstream units
            for dest_unit = 1:network_units
                slow_potential(dest_unit) = slow_potential(dest_unit) + connection_strength(unit, dest_unit);
                fast_potential(dest_unit) = fast_potential(dest_unit) + connection_strength(unit, dest_unit);
            end
        end
    end % for unit
    
    % Check for units that fire when threshold is reached
    
    for unit = 1:network_units
        unit_potentials(clock, unit) = slow_potential(unit) - fast_potential(unit);
        if unit_potentials(clock, unit) >= unit_threshold
            % Destination unit fires.
            unit_spike_queue(unit) = bitor(unit_spike_queue(unit), 1); % Start destination unit spike
            fast_potential(unit) = 0;  % Reset unit potential
            slow_potential(unit) = 0;
            if unit <= output_units
                output_events(clock, unit) = 1;  % Output event occurred
                disp(['Clock ' num2str(clock)  ', Output ' num2str(unit)]);
            end
        end
    end % for unit
    
end % for clock

