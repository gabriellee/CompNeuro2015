%FinalProject
%Gabrielle Ewall, Kirsi Goldynia
%4-16-15

%determine the effect of varying inhibition on ocular dominance competition

%DEFINE PARAMETERS
t_end = 100;%ms
dt =.1;
tau_e = 20%ms
tau_weights = 20;
w_post_pre = zeros(1,t_end);%weight of the presynaptic cell onto the postsynaptic cell
w_post_pre(1) = .5;
stim_in = .5;%c_c: strength of input from retinal ganglion cells
t_vec = 0:dt:t_end;
r_ex =zeros(1,length(t_end))
r_mean = 5;%hz
A_plus = 0.003;%magnitude of LTP
A_minus = .003/.99;%magnitude of LTD
delta_w = 0;

%define temporal change in PSPs (epsilon)
%generate poisson spike distribution using the mean rate to determine when
%the next spike will occur
spike_rate_vec1 = poissrnd(r_mean, 1, t_end*dt);
for i=1:length(spike_rate_vec1)
    spike_vec1(i + (1/spike_rate_vec1(i))/dt) = 1;
    if length(spike_vec1) > t_end/dt
        return
    end
end

%finding rates
for i=1:t_end*dt
    temp_sum = 0;
    t = t_vec(i);
    %changing the value of PSP_var (E)
    if t>=0
        tau_e = 20%ms
        PSP_var(i) = (t/tau_e^2)*e^(-t/tau_e);
        
    else
        PSP_var(i) = 0;
    end
    %changing weights
    if t > 0:
        delta_w = A_plus * exp(-dt/tau_weights);
    elseif t < 0:
        delta_w = -A_minus * exp(dt/tau_weights);
    else:
        delta_w = 0;
    end
    
    for f = 0:sum(spike_vec1(i))
        temp_sum = temp_sum + PSP_var(i - (1/spike_rate_vec1(i))/dt) ...
            + r_ex_base;
    end
    
    r_ex(i) = temp_sum * stim_in;
end
    
%cover the second group of excitatory
%do basically the same for inhibitory
%integrate over #spikes and #synapses

%define the changing weights
%postsynaptic rates
%integrate over #spikes and #synapses




%run experiment
