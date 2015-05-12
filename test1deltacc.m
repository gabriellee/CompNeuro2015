%FinalProject
%Gabrielle Ewall, Kirsi Goldynia
%4-16-15

% SIMPLIFYING ASSUMPTION: only one left input, only one right input, only
% one cell in each group

%generate poisson spike trains for the first layer of input from each
%eye
%calculate the firing rate at each point in time of both 'Excitatory Inputs'
    %this is done using the first layer poisson spike trains
%generate poisson trains for each group of layer 2 'Excitatory Inputs'
%   using the rates
%use spike times to find conductance values vector
%plug conductances into v=IR for both groups of excitatory, andinhibitory
%inputs onto the postsynaptic cell
%calculate the total current onto postsynaptic cell
%plug into the LIF update rule

%do monocular deprivation

%determine the effect of varying inhibition on ocular dominance competition
clear
%DEFINE PARAMETERS
%t_dep*5 = 5000;%ms
t_dep = 500*1000;
dt =.1;%ms
tau_e = 20;%ms
tau_weights = 20;
w_L = zeros(1,t_dep*5/dt);%weight of the presynaptic cell onto the postsynaptic cell
w_L(1) = .5;
w_R = zeros(1,t_dep*5/dt);
w_R(1) = .5;
%PSP_var_0 = (0/tau_e^2)*exp(0/tau_e);
tau_m = 20;
delta_cc = .5;
cf_test_vec = [0 .2 .5 .8 1];
delta_cc_test = [.2 .4 .5];

%USER MODIFIED INPUTS

cc = .5;%c_c: strength of input from retinal ganglion cells
cf = .5;%level of inhibitory input to postsynaptic cell

t_vec = 0:dt:t_dep*5;
r_ex_L =zeros(1,length(t_dep*5));
r_ex_R =zeros(1,length(t_dep*5));
r_mean = 5;%hz
%rate_vec = r_mean * exp(-.5 * ((angle_vec - theta_max)./sigma_r) .^2);
A_plus = 0.003;%magnitude of LTP
A_minus = .003/.99;%magnitude of LTD
delta_w = 0;
delta_g_ex = .6;% made up
delta_g_inh = .5;% made up
r_base = 10*(1-cf); %base inhibitory rate 
num_ex = 1; %number of excitatory synapses onto inhibitory/postsynaptic cell
num_inh = 1; %number of inhibitory synapses onto postsynaptic cell
gbar_ex = .015;
tau_ex = 5;
gbar_inh = .005;
tau_inh = 10; %ms
E_ex = 0;%mV
E_inh = -80;%mV
g_leak = 1;
E_leak = -74;%mV
t_dep = 5*1000;
tau_e = 20;%ms
spike_train_L_lay2 = zeros(1,t_dep*5/dt);
spike_train_R_lay2 = zeros(1,t_dep*5/dt);


%changing the value of PSP_var (E)
PSP_var = zeros(length(t_dep*5/dt));
PSP_var(t_vec >= 0) = (t_vec(t_vec >= 0)/tau_e^2).*exp(-t_vec(t_vec >= 0)/tau_e);
PSP_var = [(0/tau_e^2)*exp(0/tau_e) PSP_var];


%GENERATE POISSON MODEL SPIKE TRAIN from LEFT EYE
%figure;
num_vec = rand(1,t_dep*5/dt);
spike_train_vec_L((r_mean/1000)*dt > num_vec) = 1;
%plot(dt:dt:t_dep*5, spike_train_vec_L)

%GENERATE POISSON MODEL SPIKE TRAIN from RIGHT EYE
%figure;
num_vec = rand(1,t_dep*5/dt);
spike_train_vec_R((r_mean/1000)*dt > num_vec) = 1;
%plot(dt:dt:t_dep*5, spike_train_vec_R)
delta_cc = delta_cc_test(1)
tic
cc_vec_L = zeros(1,t_dep*5/dt);
cc_vec_R = zeros(1,t_dep*5/dt);

cc_vec_L(1:t_dep/dt) = cc;
cc_vec_L(t_dep/dt:t_dep*2/dt) = cc - delta_cc;
cc_vec_L(t_dep*2/dt:t_dep*3/dt) = cc;
cc_vec_L(t_dep*3/dt:t_dep*4/dt) = cc;
cc_vec_L(t_dep*4/dt:t_dep*5/dt) = cc;

cc_vec_R(1:t_dep/dt) = cc;
cc_vec_R(t_dep/dt:t_dep*2/dt) = cc;
cc_vec_R(t_dep*2/dt:t_dep*3/dt) = cc;
cc_vec_R(t_dep*3/dt:t_dep*4/dt) = cc - delta_cc;
cc_vec_R(t_dep*4/dt:t_dep*5/dt) = cc;



%finding rates

for i=1:(t_dep*5/dt - 1)
    temp_sum = 0;
    t = t_vec(i);


    sum_ex_L = 0;
    sum_ex_R = 0;
    r_ex_L(i) = sum(PSP_var(i - spike_train_vec_L(1:i).*[1:i]+1) + r_base) * cc_vec_L(i);
    %spike_train_vec_L(1:i)*[1:i] results in a vector with the indices
    %of spikes
    %add one to the index so that you don't have to index at zero when
    %the spike is at time i*dt


    r_ex_R(i) = sum(PSP_var(i - spike_train_vec_R(1:i).*[1:i]+1) + r_base) * cc_vec_R(i);

end
disp('139')
toc