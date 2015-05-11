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
PSP_var_0 = (0/tau_e^2)*exp(0/tau_e);
tau_m = 20;
delta_cc = .5;
cf_test_vec = [0 .2 .5 .8 1];
delta_cc_test = [.2 .4 .5];

%USER MODIFIED INPUTS

cc = .5;%c_c: strength of input from retinal ganglion cells
cf = .5;%level of inhibitory input to postsynaptic cell

t_vec = 0:dt:t_dep*5;
r_ex_L =zeros(1,length(t_dep*5));
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

%define temporal change in PSPs (epsilon)
%generate poisson spike distribution using the mean rate to determine when
%the next spike will occur
%use code from PS4
for delta_cc = delta_cc_test
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

    
    %finding rates

    for i=1:(t_dep*5/dt - 1)
        temp_sum = 0;
        t = t_vec(i);

       
        sum_ex_L = 0;
        sum_ex_R = 0;
        for f = 1:i
            if spike_train_vec_L(f) == 1
                if f == i
                    sum_ex_L = PSP_var_0 + r_base + sum_ex_L;
                else
                    sum_ex_L = PSP_var(i - f) + r_base + sum_ex_L;
                end
            end
            if spike_train_vec_R(f) == 1
                if f == i
                    sum_ex_R = PSP_var_0 + r_base + sum_ex_R;
                else
                    sum_ex_R = PSP_var(i - f) + r_base + sum_ex_R;
                end
            end
            disp('139')
            toc
        end
        disp('149')
        toc

        r_ex_L(i) = sum_ex_L * cc_vec_L(i);
        r_ex_R(i) = sum_ex_R * cc_vec_R(i);
    end

    %GENERATE POISSON MODEL SPIKE TRAIN from LEFT EX INPUTS
    %figure;
    num_vec = rand(1,t_dep*5/dt);
    spike_train_L_lay2((r_ex_L/1000)*dt > num_vec) = 1;
    %plot(dt:dt:t_dep*5, spike_train_vec_L)

    %GENERATE POISSON MODEL SPIKE TRAIN from RIGHT EX INPUTS
    %figure;
    num_vec = rand(1,t_dep*5/dt);
    spike_train_R_lay2((r_ex_R/1000)*dt > num_vec) = 1;
    %plot(dt:dt:t_dep*5, spike_train_vec_L)
        
    for i = 1:(t_dep*5/dt -1)

        %calculate RATE OF INHIBITORY GROUP
        sum_inh_L = 0;
        sum_inh_R = 0;
        for syn = 1:num_ex
            for f = 1:i
                if spike_train_L_lay2(f) ==1
                    if f == i
                        sum_inh_L = PSP_var_0 + r_base + sum_inh_L;
                    else
                        sum_inh_L = PSP_var(i - f) + r_base + sum_inh_L;
                    end
                end

                if spike_train_R_lay2(f) == 1
                    if spike_train_L_lay2(f) == 1
                        if f == i
                            sum_inh_R = PSP_var_0 + r_base + sum_inh_R;
                        else
                            sum_inh_R = PSP_var(i - f) + r_base + sum_inh_R;
                        end
                    end
                end
            end
        end

        sum_inh = sum_inh_L +sum_inh_R;
        r_inh(i) = sum_inh *cf/(num_ex*2); %2 because 2 groups of excitatory inputs
        %it makes sense that there are more spikes later on, because as t is
        %farther from arrival time of f, epsilon gets smaller, so at the
        %beginning of the spike train, you have rate =  a bunch of higher
        %values/#syn.  Later, you have rate = even more higher values (b/c lay2 ex cells also increase rate for the same reasons) and many
        %low values (from long past spikes)/#syn

        %GENERATE POISSON MODEL SPIKE TRAIN from INH INPUTS
        %figure;
        spike_train_inh = zeros(1,t_dep*5/dt);
        num_vec = rand(1,t_dep*5/dt);
        spike_train_inh((r_inh/1000)*dt > num_vec(1:i)) = 1;
        plot(dt:dt:t_dep*5, spike_train_inh)
    end



    %CALCULATE POSTSYNAPTIC VOLTAGE
    %assume that a change in weight only affects the next time step
    for i = 1:t_dep*5/dt - 1
        V = zeros(1,t_dep*5/dt);
        V(1) = -60;%mV, made up
        for synE = 1:num_ex

            if V(i) >= -54 %mV
                V(i+1) = -60; %mV

                %LEFT WEIGHT
                delta_t = dt* (find(fliplr(spike_train_L_lay2),1) - 1);
                if delta_t == 0
                    delta_w = 0;
                elseif delta_t > 0
                    delta_w = A_plus * exp(-delta_t/tau_weights);
                elseif delta_t < 0
                    delta_w = -A_minus*exp(delta_t/tau_weights);
                w_L(i+1) = delta_w + w_L(i);
                delta_w
                end


                %RIGHT WEIGHT
                delta_t = dt* (find(fliplr(spike_train_R_lay2),1) - 1);
                if delta_t == 0
                    delta_w = 0;
                elseif delta_t > 0
                    delta_w = A_plus * exp(-delta_t/tau_weights);
                elseif delta_t < 0
                    delta_w = -A_minus*exp(delta_t/tau_weights);
                w_R(i+1) = delta_w + w_R(i);
                end
            end

            for synI = 1:num_inh
                %calculate gs (for each ex and inh synapse)
                g_ex_L = gbar_ex * w_L(i +1) * exp(-t/tau_ex);
                g_ex_R = gbar_ex * w_R(i+1) * exp(-t/tau_ex);

                g_inh = gbar_inh * (exp(1)/tau_inh)*t*exp(-t/tau_inh);


                %calculate Is (for each ex and inh synapse)
                I_ex_L = -g_ex_L *(V(i+1) -E_ex);
                I_ex_R = -g_ex_R * (V(i+1) -E_ex);
                I_inh = -g_inh * (V(i+1) -E_inh);
            end

        end
        %calculate voltage!!!
        I_total = I_inh+I_ex_L +I_ex_R;
        V_inf = -I_total/(g_leak*(-E_leak + 1));
        tau_eff = tau_m/(g_leak*(-E_leak + 1));
        V(i+1) = V_inf + (V(1) - V_inf)*exp(-t/tau_eff);

        %strategy: use first weight to calculate initial voltage and see if
        %spike.  Kepp doing this until there is a spike, then use the weight
        %update rule



        %define change in conductances here
        %define V of the postsynaptic cell here
    end
    figure;
    hold on
    plot(0:dt:t_dep*5, w_L)
    plot(0:dt:t_dep*5, w_R)
    legend('Left','Right')
    title(strcat('delta cc = ',num2str(delta_cc)));
    figure;
    hold on;
    plot(0:dt:t_dep*5, w_L,'o')
    plot(0:dt:t_dep*5, w_R,'o')
    legend('Left','Right')
end


delta_cc = .5;
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
for cf = cf_test_vec
    %GENERATE POISSON MODEL SPIKE TRAIN from LEFT EYE
    %figure;
    spike_train_vec_L = zeros(1,t_dep*5/dt);
    num_vec = rand(1,t_dep*5/dt);
    spike_train_vec_L((r_mean/1000)*dt > num_vec) = 1;
    %plot(dt:dt:t_dep*5, spike_train_vec_L)

    %GENERATE POISSON MODEL SPIKE TRAIN from RIGHT EYE
    %figure;
    spike_train_vec_R = zeros(1,t_dep*5/dt);
    num_vec = rand(1,t_dep*5/dt);
    spike_train_vec_R((r_mean/1000)*dt > num_vec) = 1;
    %plot(dt:dt:t_dep*5, spike_train_vec_R)


    %finding rates
    for i=1:(t_dep*5/dt - 1)
        temp_sum = 0;
        t = t_vec(i);
        %changing weights
        if t > 0
            delta_w = A_plus * exp(-dt/tau_weights);%delta t is btw pre and psot
        elseif t < 0
            delta_w = -A_minus * exp(dt/tau_weights);
        else
            delta_w = 0;
        end

        sum_ex_L = 0;
        sum_ex_R = 0;
        for f = 1:i
            if spike_train_vec_L(f) == 1
                if f == i
                    sum_ex_L = PSP_var_0 + r_base + sum_ex_L;
                else
                    sum_ex_L = PSP_var(i - f) + r_base + sum_ex_L;
                end
            end
            if spike_train_vec_R(f) == 1
                if f == i
                    sum_ex_R = PSP_var_0 + r_base + sum_ex_R;
                else
                    sum_ex_R = PSP_var(i - f) + r_base + sum_ex_R;
                end
            end
        end

        r_ex_L(i) = sum_ex_L * cc_vec_L;
        r_ex_R(i) = sum_ex_R * cc_vec_R;

        %GENERATE POISSON MODEL SPIKE TRAIN from LEFT EX INPUTS
        %figure;
        spike_train_L_lay2 = zeros(1,t_dep*5/dt);
        num_vec = rand(1,t_dep*5/dt);
        spike_train_L_lay2((r_ex_L/1000)*dt > num_vec(1:i)) = 1;
        %plot(dt:dt:t_dep*5, spike_train_vec_L)

        %GENERATE POISSON MODEL SPIKE TRAIN from RIGHT EX INPUTS
        %figure;
        spike_train_R_lay2 = zeros(1,t_dep*5/dt);
        num_vec = rand(1,t_dep*5/dt);
        spike_train_R_lay2((r_ex_R/1000)*dt > num_vec(1:i)) = 1;
        %plot(dt:dt:t_dep*5, spike_train_vec_L)

        %calculate RATE OF INHIBITORY GROUP
        sum_inh_L = 0;
        sum_inh_R = 0;
        for syn = 1:num_ex
            for f = 1:i
                if spike_train_L_lay2(f) ==1
                    if f == i
                        sum_inh_L = PSP_var_0 + r_base + sum_inh_L;
                    else
                        sum_inh_L = PSP_var(i - f) + r_base + sum_inh_L;
                    end
                end

                if spike_train_R_lay2(f) == 1
                    if spike_train_L_lay2(f) == 1
                        if f == i
                            sum_inh_R = PSP_var_0 + r_base + sum_inh_R;
                        else
                            sum_inh_R = PSP_var(i - f) + r_base + sum_inh_R;
                        end
                    end
                end
            end
        end

        sum_inh = sum_inh_L +sum_inh_R;
        r_inh(i) = sum_inh *cf/(num_ex*2); %2 because 2 groups of excitatory inputs
        %it makes sense that there are more spikes later on, because as t is
        %farther from arrival time of f, epsilon gets smaller, so at the
        %beginning of the spike train, you have rate =  a bunch of higher
        %values/#syn.  Later, you have rate = even more higher values (b/c lay2 ex cells also increase rate for the same reasons) and many
        %low values (from long past spikes)/#syn

        %GENERATE POISSON MODEL SPIKE TRAIN from INH INPUTS
        %figure;
        spike_train_inh = zeros(1,t_dep*5/dt);
        num_vec = rand(1,t_dep*5/dt);
        spike_train_inh((r_inh/1000)*dt > num_vec(1:i)) = 1;
        plot(dt:dt:t_dep*5, spike_train_inh)
    end



    %CALCULATE POSTSYNAPTIC VOLTAGE
    %assume that a change in weight only affects the next time step
    for i = 1:t_dep*5/dt - 1
        V = zeros(1,t_dep*5/dt);
        V(1) = -60;%mV, made up
        for synE = 1:num_ex

            if V(i) >= -54 %mV
                V(i+1) = -60; %mV

                %LEFT WEIGHT
                delta_t = dt* (find(fliplr(spike_train_L_lay2),1) - 1);
                if delta_t == 0
                    delta_w = 0;
                elseif delta_t > 0
                    delta_w = A_plus * exp(-delta_t/tau_weights);
                elseif delta_t < 0
                    delta_w = -A_minus*exp(delta_t/tau_weights);
                w_L(i+1) = delta_w + w_L(i);
                end


                %RIGHT WEIGHT
                delta_t = dt* (find(fliplr(spike_train_R_lay2),1) - 1);
                if delta_t == 0
                    delta_w = 0;
                elseif delta_t > 0
                    delta_w = A_plus * exp(-delta_t/tau_weights);
                elseif delta_t < 0
                    delta_w = -A_minus*exp(delta_t/tau_weights);
                w_R(i+1) = delta_w + w_R(i);
                end
            end

            for synI = 1:num_inh
                %calculate gs (for each ex and inh synapse)
                g_ex_L = gbar_ex * w_L(i +1) * exp(-t/tau_ex);
                g_ex_R = gbar_ex * w_R(i+1) * exp(-t/tau_ex);

                g_inh = gbar_inh * (exp(1)/tau_inh)*t*exp(-t/tau_inh);


                %calculate Is (for each ex and inh synapse)
                I_ex_L = -g_ex_L *(V(i+1) -E_ex);
                I_ex_R = -g_ex_R * (V(i+1) -E_ex);
                I_inh = -g_inh * (V(i+1) -E_inh);
            end

        end
        %calculate voltage!!!
        I_total = I_inh+I_ex_L +I_ex_R;
        V_inf = -I_total/(g_leak*(-E_leak + 1));
        tau_eff = tau_m/(g_leak*(-E_leak + 1));
        V(i+1) = V_inf + (V(1) - V_inf)*exp(-t/tau_eff);

        %strategy: use first weight to calculate initial voltage and see if
        %spike.  Kepp doing this until there is a spike, then use the weight
        %update rule



        %define change in conductances here
        %define V of the postsynaptic cell here
    end
    figure;
    hold on
    plot(0:dt:t_dep*5, w_L)
    plot(0:dt:t_dep*5, w_R)
    legend('Left','Right')
    title(strcat('cf = ',num2str(cf)));
    figure;
    hold on;
    plot(0:dt:t_dep*5, w_L,'o')
    plot(0:dt:t_dep*5, w_R,'o')
    legend('Left','Right')
end

%V  tracks the post synaptic cell
    
%cover the second group of excitatory
%do basically the same for inhibitory
%integrate over #spikes and #synapses

%use the cha

%define the changing weights
%postsynaptic rates
%integrate over #spikes and #synapses




%run experiment
