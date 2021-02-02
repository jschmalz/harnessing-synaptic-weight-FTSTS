function [time v_E v_I S_EI S_IE X_EI X_IE Apost Apre W_IE spike_E spike_I ref_E ref_I synchrony,spike_E_time,spike_I_time,ISI_EI] = ode_neuron_model(plast_on,ON1,vE0,vI0,S_EI0,S_IE0,X_EI0,X_IE0,Apost0,Apre0,W_IE0,W_EI0,mew_e,sigma_e,ue,ui,mew_i,sigma_i,J_E,J_I,C_E,C_I,tau_LTP,tau_LTD,step,sample_duration,N_E,N_I,S_key_EI,S_key_IE,leftover_S_EI,leftover_S_IE,ref_E,ref_I,Energy0,tau_E_m,tau_I_m, percent_V_stim,spike_E_time0,spike_I_time0,ISI_EI0,comp_time)
% start timer
    
% -- run parameters --
    duration_step = sample_duration/step;
 
% -- neuron parameters --
    tau_d = 1;
    tau_r = 1;
    vrest = 0; %mV
    whitenoise_E = randn(duration_step+1,N_E); % gaussian white noise
    whitenoise_I = randn(duration_step+1,N_I); % gaussian white noise
    vreset = 14;
    refractory = 2; % ms     
  
     
% -- Synaptic Parameters --    
    WEI = W_EI0;
    syn_delay = 5; %ms
    num_synapses_IE = max(max(S_key_IE));
    num_synapses_EI = max(max(S_key_EI));
     
% -- Plasticity Parameters --
    dApre_0 = 0.005*1;
    dApost_0 = dApre_0*1;
    Apre_i = 0;
    Apost_i = 0;
    a_LTD = -1*plast_on*1.05;
    a_LTP = 1*plast_on*1;
    eta = 1;
    cP = 0.038*plast_on;
    cD = 0.02*plast_on;
    tauP = 10;
    tauD = 25;
       
% -- intialize vectors --
    % ode vectors
    dv_Edt = zeros(duration_step+1,N_E);
    dv_Idt = zeros(duration_step+1,N_I);
    dS_EIdt = zeros(duration_step+1,N_E);
    dS_IEdt = zeros(duration_step+1,N_I);
    dX_EIdt = zeros(duration_step+1,N_E);
    dX_IEdt = zeros(duration_step+1,N_I);
    dApostdt = zeros(duration_step+1,num_synapses_IE);
    dApredt = zeros(duration_step+1,num_synapses_IE);
     
    % state vectors    
    v_E = zeros(duration_step+1,N_E);
    v_I = zeros(duration_step+1,N_I);     
    S_EI = zeros(duration_step+1,N_E);
    S_IE = zeros(duration_step+1,N_I);    
    X_EI = zeros(duration_step+1,N_E);
    X_IE = zeros(duration_step+1,N_I);    
    Apost = zeros(duration_step+1,num_synapses_IE);
    Apre = zeros(duration_step+1,num_synapses_IE);    
    W_IE = zeros(duration_step+1,num_synapses_IE);    
    time = zeros(duration_step+1,1);
    spike_E = zeros(duration_step+1,N_E);
    spike_I = zeros(duration_step+1,N_I);
    spike_E_time = zeros(duration_step+1,N_E);
    spike_I_time = zeros(duration_step+1,N_I);
    ISI_EI = zeros(duration_step+1,num_synapses_IE);
    
     
    v_E(1,:) = vE0;
    v_I(1,:) = vI0;
     
    S_EI(1,:) = S_EI0;
    S_IE(1,:) = S_IE0;
     
    X_EI(1,:) = X_EI0;
    X_IE(1,:) = X_IE0;
     
    Apost(1,:) = Apost0;
    Apre(1,:) = Apre0;
     
    W_IE(1,:) = W_IE0;
    
    spike_E_time(1,:) = spike_E_time0;
    spike_I_time(1,:) = spike_I_time0;
    
    ISI_EI(1,:) = ISI_EI0;
    
for i = 1:duration_step
% -- time update --
    time(i+1,1) = time(i,1) + step;

% -- excitatory --
%     for k = 1:N_E
        % dxdt update
        if i > syn_delay/step
            dv_Edt(i,:) = (vrest - v_E(i,:) + J_E/C_E*S_EI(i-syn_delay/step,:) + mew_e + ON1.*ue(1,i) + sigma_e.*(tau_E_m(1,:).^0.5).*whitenoise_E(i,:) )./tau_E_m(1,:);
        else
            dv_Edt(i,:) = (vrest - v_E(i,:) + J_E/C_E*leftover_S_EI(i,:) + mew_e + ON1.*ue(1,i) + sigma_e.*(tau_E_m(1,:).^0.5).*whitenoise_E(i,:) )./tau_E_m(1,:);
        end
        dS_EIdt(i,:) = (-S_EI(i,:) + X_EI(i,:))/tau_d;
        dX_EIdt(i,:) = -X_EI(i,:)/tau_r;

        
        
        % x update
        v_E(i+1,:) = v_E(i,:) + step*dv_Edt(i,:);
        S_EI(i+1,:) = S_EI(i,:) + step*dS_EIdt(i,:);
        X_EI(i+1,:) = X_EI(i,:) + step*dX_EIdt(i,:);
         
        % update refractory
        ref_E(1,:) = ref_E(1,:) - step;
%     end
         
     
% -- inhibitory --
   
%         dxdt update
        if i > syn_delay/step
            dv_Idt(i,:) = (vrest - v_I(i,:) + J_I/C_I*S_IE(i-syn_delay/step,:) + mew_i + ON1.*ui(1,i) + sigma_i.*(tau_I_m(1,:).^0.5).*whitenoise_I(i,:) )./tau_I_m(1,:);
        else
            dv_Idt(i,:) = (vrest - v_I(i,:) + J_I/C_I*leftover_S_IE(i,:) + mew_i + ON1.*ui(1,i) + sigma_i.*(tau_I_m(1,:).^0.5).*whitenoise_I(i,:) )./tau_I_m(1,:);
        end
        dS_IEdt(i,:) = (-S_IE(i,:) + X_IE(i,:))/tau_d;
        dX_IEdt(i,:) = -X_IE(i,:)/tau_r;   
        
        
        % x update
        v_I(i+1,:) = v_I(i,:) + step*dv_Idt(i,:);
        S_IE(i+1,:) = S_IE(i,:) + step*dS_IEdt(i,:);
        X_IE(i+1,:) = X_IE(i,:) + step*dX_IEdt(i,:);        
         
        % update refractory
        ref_I(1,:) = ref_I(1,:) - step;
%     end
     
% --- update synaptic variables ---
 
     
    % dxdt update
        dApostdt(i,:) = -Apost(i,:)/tau_LTP;
        dApredt(i,:) = -Apre(i,:)/tau_LTD;
     
    % x update
        Apost(i+1,:) = Apost(i,:) + step*dApostdt(i,:);
        Apre(i+1,:) = Apre(i,:) + step*dApredt(i,:);
        W_IE(i+1,:) = W_IE(i,:);
        spike_E_time(i+1,:) = spike_E_time(i,:);
        spike_I_time(i+1,:) = spike_I_time(i,:);
        ISI_EI(i+1,:) = ISI_EI(i,:);
%     
% -- check for action potentials --
 
 
    % -- excitatory spike --
    for k = 1:N_E
        if v_E(i+1,k) >= 20 && v_E(i,k) < 20
             
            % reset & refractory
            v_E(i+1,k) = vreset;
            ref_E(1,k) = refractory;
             
            % spike monitor
            spike_E(i+1,k) = k;
            spike_E_time(i+1,k) = time(i+1,1) + comp_time;
             
            % check synaptic connection E to I
            for j = 1:N_I
                if S_key_IE(k,j) ~= 0
                    index = S_key_IE(k,j);
                    % synaptic input from E to I : _(IE)
                    X_IE(i+1,j) = X_IE(i,j) + W_IE(i,index);
                     
                    % plasticity update -
                    ISI_EI(i+1,index) = spike_E_time(i+1,k) - spike_I_time(i+1,j);                    
                    deltaS_EI = cP*exp(-ISI_EI(i+1,index)/tauP) - cD*exp(-ISI_EI(i+1,index)/tauD);
                    
                    W_IE(i+1,index) = W_IE(i,index) + deltaS_EI;
                    
                    % max synaptic weight check                       
                    if (J_I*W_IE(i+1,index)) < 10
                        W_IE(i+1,index) = 10/J_I;
                    end
                    % max synaptic weight check
                    if (J_I*W_IE(i+1,index)) > 290
                        W_IE(i+1,index) = 290/J_I;                   
                    end
                end
            end
             
        elseif ref_E(1,k) >= 0
             
            % check if in refractory period
            v_E(i+1,k) = vreset;
        elseif v_E(i+1,k) < vrest
            v_E(i+1,k) = vrest;
             
        end
    end
     
    % -- inhibitory neuron spike --
    for k = 1:N_I
        if v_I(i+1,k) >= 20 && v_I(i,k) < 20
             
            % reset & refractory
            v_I(i+1,k) = vreset;
            ref_I(1,k) = refractory;
             
            % spike monitor
            spike_I(i+1,k) = k+N_E;
            spike_I_time(i+1,k) = time(i+1,1) + comp_time;
             
            % check synaptic connection I to E
            for j = 1:N_E
                % synaptic input from I to E : _(EI)
                if S_key_EI(k,j) ~=0  
                    index = S_key_EI(k,j);
                    X_EI(i+1,j) = X_EI(i,j) - WEI;
                end    
                
                % add start
                % check synaptic connection E to I
            
                if S_key_IE(j,k) ~= 0
                    index = S_key_IE(j,k);
                   
                    % plasticity update -
                    ISI_EI(i+1,index) = spike_I_time(i+1,k) - spike_E_time(i+1,j);                    
                    deltaS_EI = cP*exp(-abs(ISI_EI(i+1,index))/tauP) - cD*exp(-abs(ISI_EI(i+1,index))/tauD);
                    
                    W_IE(i+1,index) = W_IE(i,index) + deltaS_EI;
                    
                    % max synaptic weight check                       
                    if (J_I*W_IE(i+1,index)) < 10
                        W_IE(i+1,index) = 10/J_I;
                    end
                    % max synaptic weight check
                    if (J_I*W_IE(i+1,index)) > 290
                        W_IE(i+1,index) = 290/J_I;                   
                    end
                end
                % add end
                
                
            end
             
        elseif ref_I(1,k) >= 0
             
            % check if in refractory period
            v_I(i+1,k) = vreset;
             
        elseif v_I(i+1,k) < vrest
            v_I(i+1,k) = vrest;
        end
    end
     
end
% --- Synchrony calculation ---

    N = N_E;
    Vcomb = zeros(sample_duration/step+1,N);
    Vcomb(:,1:N_E) = v_E;
    V1 = mean(Vcomb,2); 
    
    % variance of average volatage over whole run
    sigma_squ_v = mean(V1.^2) - (mean(V1))^2;   % sigma_v^2 = <(V(t)^2>t -[<V(t)>]^2
    
    % variance of volatage at each time step
    sigma_vi = zeros(1,N);
    sum_sig = 0;
    
    for i = 1:N
       sigma_vi(i) = mean(Vcomb(:,i).^2) - (mean(Vcomb(:,i)))^2; 
       sum_sig = sum_sig + sigma_vi(i);
    end
   
    % calculate synchrony
    syn_squ = sigma_squ_v/(sum_sig/N);
    synchrony = syn_squ^0.5;
 
end