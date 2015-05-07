%FinalProject
%Gabrielle Ewall, Kirsi Goldynia
%4-16-15

%generate poisson spike trains for the first layer of input from each
%eye
%calculate the firing rate at each point in time of both 'Excitatory Inputs'
    %this is done using the first layer poisson spike trains
%generate poisson trains for each group of layer 2 'Excitatory Inputs'
%using the rates
%use spike times to find conductance values vector
%plug conductances into v=IR for both groups of excitatory, andinhibitory
%inputs onto the postsynaptic cell
%calculate the total current onto postsynaptic cell
%plug into the LIF update rule

%determine the effect of varying inhibition on ocular dominance competition
clear
%DEFINE PARAMETERS
t_end = 5000;%ms
dt =.1;%ms
tau_e = 20%ms
tau_weights = 20;
w_post_pre = zeros(1,t_end);%weight of the presynaptic cell onto the postsynaptic cell
w_post_pre(1) = .5;
stim_in = .5;%c_c: strength of input from retinal ganglion cells
t_vec = 0:dt:t_end;
r_ex =zeros(1,length(t_end))
r_mean = 5;%hz
%rate_vec = r_mean * exp(-.5 * ((angle_vec - theta_max)./sigma_r) .^2);
A_plus = 0.003;%magnitude of LTP
A_minus = .003/.99;%magnitude of LTD
delta_w = 0;
delta_g_ex = .6;% made up
delta_g_inh = .5;% made up

%define temporal change in PSPs (epsilon)
%generate poisson spike distribution using the mean rate to determine when
%the next spike will occur
%use code from PS4

%GENERATE POISSON MODEL SPIKE TRAIN from LEFT EYE
%figure;
spike_train_vec_L = zeros(1,t_end/dt);
num_vec = rand(1,t_end/dt);
spike_train_vec_L((r_mean/1000)*dt > num_vec) = 1;
%plot(dt:dt:t_end, spike_train_vec_L)

%GENERATE POISSON MODEL SPIKE TRAIN from RIGHT EYE
%figure;
spike_train_vec_R = zeros(1,t_end/dt);
num_vec = rand(1,t_end/dt);
spike_train_vec_R((r_mean/1000)*dt > num_vec) = 1;
%plot(dt:dt:t_end, spike_train_vec_R)



% %finding rates
% for i=1:t_end*dt
%     temp_sum = 0;
%     t = t_vec(i);
%     %changing the value of PSP_var (E)
%     if t>=0
%         tau_e = 20%ms
%         PSP_var(i) = (t/tau_e^2)*e^(-t/tau_e);
%         
%     else
%         PSP_var(i) = 0;
%     end
%     %changing weights
%     if t > 0:
%         delta_w = A_plus * exp(-dt/tau_weights);%delta t is btw pre and psot
%     elseif t < 0:
%         delta_w = -A_minus * exp(dt/tau_weights);
%     else:
%         delta_w = 0;
%     end
%     
%     
%     for f = 0:sum(spike_vec1(i))
%         temp_sum = temp_sum + PSP_var(i - (1/spike_rate_vec1(i))/dt) ...
%             + r_ex_base;
%     end
%     
%     r_ex(i) = temp_sum * stim_in;
%     %calculate r_inh also in this loop
%     
%     %define change in conductances here
%     %define V of the postsynaptic cell here
% end
%V  tracks the post synaptic cell
    
%cover the second group of excitatory
%do basically the same for inhibitory
%integrate over #spikes and #synapses

%use the cha

%define the changing weights
%postsynaptic rates
%integrate over #spikes and #synapses




%run experiment
