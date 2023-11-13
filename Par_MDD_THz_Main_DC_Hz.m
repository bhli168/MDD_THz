function [] = Par_MDD_THz_Main_DC_Hz(Number_AP,Number_User,Regular,U_max,C_max,Bandwidth)
rng('shuffle')

Len = 100;
Wid = 100;
Roof_Height = 10;
CPU_height = Roof_Height;
User_height = 1;

N_cpu = 64; % 8*8 UPA
N_cpu_w = 8;
N_cpu_h = 8;
N_RF_cpu = 32;
N_ap = 16; % 4*4 UPA
N_ap_w = 4;
N_ap_h = 4;
N_RF_ap = 8;
Power_AP_dBm = 35;
Power_AP_W = (10^(0.1 * Power_AP_dBm) * 10^-3);
Power_CPU_dBm = 45;
Power_CPU_W = (10^(0.1 * Power_CPU_dBm) * 10^-3);

transmit_gain_dBi = 20 ; %dBi
transmit_gain = 10^(.1*transmit_gain_dBi) ;
receive_gain_dBi = 20 ; %dBi
receive_gain = 10^(.1*receive_gain_dBi) ;


Central_frequency = 200e9; % 200GHz
NumberOfSub = 32;
%Bandwidth = 1e9; % 1GHz
No_dB_THz = -174+10*log10(Bandwidth); % Noise power in dBm, -174 dBm/Hz
No_THz = 10^(.1*No_dB_THz) * 10^-3;
% No_dB = -75; %dBm
% No = 10^(.1*No_dB) * 10^-3;
SI = -10; %dB
SI_W = 10^(SI/10) * No_THz;
NumberOfDelay = 6;
NumberOfNLOS = 3;
SubTHz_Bandwidth = Bandwidth / NumberOfSub;
K_abs = 0.0033;
Fresnel = 0.15;
Sample_t = 7.8*1e-12;
roll_off = 1;
CP = 16;
Sub6G_Bandwidth = 1e8; % 100MHz
No_dB_6G = -174+10*log10(Sub6G_Bandwidth); % Noise power in dBm, -174 dBm/Hz
No_6G = 10^(.1*No_dB_6G) * 10^-3;
Sub_Fre = [];
for ii = 1:NumberOfSub/2
    Sub_Fre = [Sub_Fre Central_frequency+(SubTHz_Bandwidth/2)*ii];
    Sub_Fre = [Central_frequency-(SubTHz_Bandwidth/2)*ii Sub_Fre];
end

User_location = PLsetup(Number_User,Len,Wid,'user');
[AP_location,AP_height] = PLsetup(Number_AP,Len,Wid,'AP');
CPU_location = Len/2 + 1i * Wid/2;
F = zeros(NumberOfSub,NumberOfSub);
for m = 1:NumberOfSub
    for n = 1:NumberOfSub
        F(m,n) = exp(-1i * 2 * pi * (m-1) * (n-1) / NumberOfSub); %FFT
    end
end

parfor (par_iter = 1:40)


    %% Channel generation

    %% CPU-to-AP
    H_CC_t = cell(1,Number_AP);
    H_CC_f = cell(1,Number_AP);
    CC_complex_gain = zeros(Number_AP,1);
    for ll = 1:Number_AP %K_f=0.0033
        loss = zeros(1,NumberOfSub);
        distance = sqrt(abs(CPU_location-AP_location(ll))^2+(CPU_height-AP_height(ll))^2);
        AAoD_LoS = 2*pi*rand(1)-pi;
        EAoD_LoS = pi*rand(1)-pi/2;
        LoS_array_re = repmat([0:N_cpu_w-1],[N_cpu_h,1])*sin(AAoD_LoS)...
            + [0:N_cpu_h-1]'*sin(EAoD_LoS);
        LoS_array_re = exp(1i*LoS_array_re(:));
        LoS_delay = rand * (1*Sample_t);

        AAoD_NLoS = 2*pi*rand(1,NumberOfNLOS)-pi;
        EAoD_NLoS = pi*rand(1,NumberOfNLOS)-pi/2;
        NLoS_delay = rand(1,NumberOfNLOS) * (CP*Sample_t);
        NLoS = 0;
        for t = 0:NumberOfDelay-1

            LoS = pulse_shaping(t*Sample_t-LoS_delay,roll_off,Sample_t) * LoS_array_re;
            for n = 1:NumberOfNLOS
                NLoS_array_re = repmat([0:N_cpu_w-1],[N_cpu_h,1])*sin(AAoD_NLoS(n))...
                    + [0:N_cpu_h-1]'*sin(EAoD_NLoS(n));
                NLoS_array_re = exp(1i*NLoS_array_re(:)) * Fresnel * sqrt(1/NumberOfNLOS);
                NLoS_prc = pulse_shaping(t*Sample_t-NLoS_delay(n),roll_off,Sample_t); %% pulse shaping
                NLoS = NLoS + NLoS_array_re * NLoS_prc;
            end
            H_CC_t{1,ll}(:,t+1) = sqrt(transmit_gain*receive_gain) * (LoS + NLoS);

        end
        for ff = 1:NumberOfSub
            loss(ff) = ( 3 * 1e8 / ( 4 * pi * Sub_Fre(ff) * distance))* exp(-K_abs*distance/2);
        end
        H_CC_f{1,ll} = (diag(loss) * (F(:,1:NumberOfDelay) * H_CC_t{1,ll}.')).';
        CC_complex_gain(ll) =  ( 3 * 1e8 / ( 4 * pi * Central_frequency * distance))* exp(-K_abs*distance/2);
    end

    %% AP-to-AP
    H_CA_t = cell(Number_AP,Number_AP);
    H_CA_f = cell(Number_AP,Number_AP);
    CA_complex_gain = zeros(Number_AP,Number_AP);
    for ll = 1:Number_AP
        for lll = 1:Number_AP
            if ll==lll
                H_CA_t{ll,lll} = 0;
                H_CA_f{ll,lll} = 0;
                continue;
            end
            loss = zeros(1,NumberOfSub);
            distance = sqrt(abs(AP_location(lll)-AP_location(ll))^2+(AP_height(lll)-AP_height(ll))^2);
            AAoD_LoS = 2*pi*rand(1)-pi;
            EAoD_LoS = pi*rand(1)-pi/2;
            LoS_array_re = repmat([0:N_ap_w-1],[N_ap_h,1])*sin(AAoD_LoS)...
                + [0:N_ap_h-1]'*sin(EAoD_LoS);
            LoS_array_re = exp(1i*LoS_array_re(:));
            LoS_delay = rand * (1*Sample_t);

            AAoD_NLoS = 2*pi*rand(1,NumberOfNLOS)-pi;
            EAoD_NLoS = pi*rand(1,NumberOfNLOS)-pi/2;
            NLoS_delay = rand(1,NumberOfNLOS) * (CP*Sample_t);
            NLoS = 0;
            for t = 0:NumberOfDelay-1

                LoS = pulse_shaping(t*Sample_t-LoS_delay,roll_off,Sample_t) * LoS_array_re;
                for n = 1:NumberOfNLOS
                    NLoS_array_re = repmat([0:N_ap_w-1],[N_ap_h,1])*sin(AAoD_NLoS(n))...
                        + [0:N_ap_h-1]'*sin(EAoD_NLoS(n));
                    NLoS_array_re = exp(1i*NLoS_array_re(:)) * Fresnel * sqrt(1/NumberOfNLOS);
                    NLoS_prc = pulse_shaping(t*Sample_t-NLoS_delay(n),roll_off,Sample_t);
                    NLoS = NLoS + NLoS_array_re * NLoS_prc;
                end
                H_CA_t{ll,lll}(:,t+1) = sqrt(transmit_gain*receive_gain) * (LoS + NLoS);

            end
            for ff = 1:NumberOfSub
                loss(ff) = ( 3 * 1e8 / ( 4 * pi * Sub_Fre(ff) * distance))* exp(-K_abs*distance/2);
            end
            H_CA_f{ll,lll} = (diag(loss) * (F(:,1:NumberOfDelay) * H_CA_t{ll,lll}.')).';
            CA_complex_gain(ll,lll) = ( 3 * 1e8 / ( 4 * pi * Central_frequency * distance))* exp(-K_abs*distance/2);
        end
    end

    %% AP-to-User
    sha_F = 4;
    AU_complex_gain = zeros(Number_AP,Number_User);
    Beta_AP_MS = zeros(Number_AP,Number_User);
    AP_User_PL = (sqrt(abs(repmat(AP_location.',Number_User,1) - repmat(User_location,1,Number_AP)).^2 +...
        (repmat(AP_height.',Number_User,1)-User_height.*ones(Number_User,Number_AP)).^2)).';
    for ll = 1:Number_AP
        for d = 1:Number_User
            shadow = sha_F * randn();
            Beta_AP_MS(ll,d) = 10 ^ ((-30.5 - 36.7 * log10(AP_User_PL(ll,d))+shadow) / 10);
        end
    end
    H_DL = cell(Number_AP,Number_User);
    for m = 1:Number_AP
        for n = 1:Number_User
            H_DL{m,n} = sqrt(Beta_AP_MS(m,n)/2) .* (randn(N_ap,1) + 1i * randn(N_ap,1));
            AU_complex_gain(m,n) = H_DL{m,n}' * H_DL{m,n};
            H_DL{m,n} = H_DL{m,n} * sqrt(transmit_gain*receive_gain);
        end
    end
    Beta_AP_AP = zeros(Number_AP,Number_AP);
    AP_AP_PL = (sqrt(abs(repmat(AP_location.',Number_AP,1) - repmat(AP_location,1,Number_AP)).^2 + (repmat(AP_height.',Number_AP,1)-repmat(AP_height,1,Number_AP)).^2)).';
    AP_AP_PL_RE = 1000.*eye(Number_AP) + AP_AP_PL;


    DL_SE = zeros(1,6);
    Fronthaul_CA_MDD = zeros(1,6);
    Fronthaul_CC_MDD = zeros(1,6);
    Final_MDD = zeros(1,6);
    CC_CA_abs = zeros(6,5);
    CA_temp = zeros(6,5);
    CC_temp = zeros(6,5);

    ite = [4:4:24];
    for itee = 1:length(ite)
        Number_cluster = ite(itee);
        %% AP Clustering
        AU_cluster_gain = zeros(Number_cluster,Number_User);

        %% K-means
        Cluster_id = cell(1,Number_cluster);
        CAP_id = zeros(1,Number_cluster);
        AP_coordinate = [real(AP_location) imag(AP_location) AP_height];
        [AP_idx,Cen] = kmedoids(AP_coordinate,Number_cluster);
        for cc = 1:Number_cluster
            idx_temp = find(AP_idx==cc);
            CA_gain = sum(CA_complex_gain(idx_temp,:),2);
            CC_gain = CC_complex_gain(idx_temp);
            [~,idx_max] = max(CA_gain+CC_gain);
            CAP_id_temp = idx_temp(idx_max);
            CAP_id(cc) = CAP_id_temp;
            Cluster_id{1,cc} = [CAP_id_temp setdiff(idx_temp,CAP_id_temp).']; %% the first item is CAP
            if length(idx_temp)==1
                AU_cluster_gain(cc,:) = AU_complex_gain(CAP_id_temp,:);
            else
                AP_id_temp = setdiff(idx_temp,CAP_id_temp);
                AU_cluster_gain(cc,:) = sum(AU_complex_gain(AP_id_temp,:),1);
            end
        end
        CAP_coordinate = AP_coordinate(CAP_id,:);
        %     figure;
        %     plot(Len/2,Wid/2,'k^','MarkerSize',14)
        %     hold on
        %     plot(AP_coordinate(AP_idx==1,1),AP_coordinate(AP_idx==1,2),'r.','MarkerSize',12)
        %     hold on
        %     plot(AP_coordinate(AP_idx==2,1),AP_coordinate(AP_idx==2,2),'b.','MarkerSize',12)
        %     hold on
        %     plot(AP_coordinate(AP_idx==3,1),AP_coordinate(AP_idx==3,2),'g.','MarkerSize',12)
        %     hold on
        %     plot(AP_coordinate(AP_idx==4,1),AP_coordinate(AP_idx==4,2),'k.','MarkerSize',12)
        %     hold on
        %     plot(AP_coordinate(AP_idx==5,1),AP_coordinate(AP_idx==5,2),'r*','MarkerSize',12)
        %     hold on
        %     plot(AP_coordinate(AP_idx==6,1),AP_coordinate(AP_idx==6,2),'b*','MarkerSize',12)
        %     hold on
        %     plot(AP_coordinate(AP_idx==7,1),AP_coordinate(AP_idx==7,2),'g*','MarkerSize',12)
        %     hold on
        %     plot(AP_coordinate(AP_idx==8,1),AP_coordinate(AP_idx==8,2),'k*','MarkerSize',12)
        %     plot(CAP_coordinate(:,1),CAP_coordinate(:,2),'ro',...
        %          'MarkerSize',13,'LineWidth',3)

        %% User selection
        G_u = cell(1,Number_User);
        U_l = cell(1,Number_cluster);
        AU_cluster_gain_temp = AU_cluster_gain;
        User_gain_iterative = zeros(1,Number_User);
        for u = 1:Number_User %% Initialization
            [value_temp,C_id] = max(AU_cluster_gain_temp(:,u));
            G_u{1,u} = [G_u{1,u} C_id];
            User_gain_iterative(1,u) = value_temp;
            AU_cluster_gain_temp(C_id,u) = 0;
            U_l{1,C_id} = [U_l{1,C_id} u];
            if length(U_l{1,C_id})==U_max
                AU_cluster_gain_temp(C_id,:) = zeros(1,Number_User);
            end
        end
        while 1
            if isempty(find(AU_cluster_gain_temp>0))
                break
            end
            if min(User_gain_iterative)==1
                break
            end
            [~,U_min_id] = min(User_gain_iterative);
            if isempty(find(AU_cluster_gain_temp(:,U_min_id)>0)) || length(G_u{1,U_min_id})>=C_max
                User_gain_iterative(1,U_min_id) = 1;
            else
                [value_temp,C_id] = max(AU_cluster_gain_temp(:,U_min_id));
                G_u{1,U_min_id} = [G_u{1,U_min_id} C_id];
                AU_cluster_gain_temp(C_id,U_min_id) = 0;
                User_gain_iterative(1,U_min_id) = User_gain_iterative(1,U_min_id) + value_temp;
                U_l{1,C_id} = [U_l{1,C_id} U_min_id];
                if length(U_l{1,C_id})==U_max
                    AU_cluster_gain_temp(C_id,:) = zeros(1,Number_User);
                end
            end


        end



        CAP_AP_id = cell(1,Number_cluster);
        for cc = 1:Number_cluster %% acquire the channel gain of sub
            AP_id = Cluster_id{1,cc}(2:end);
            CAP_AP_id{1,cc} = AP_id;
        end
        %% Sub 1
        %%% RZF
        AP_Null_count = 0;
        Cluster_Null_id = [];
        H_cluster_user = cell(1,Number_cluster);
        V_RZF_cluster = cell(1,Number_cluster);
        V_MMSE_cluster = cell(1,Number_cluster);
        for cc = 1:Number_cluster
            serve_user = U_l{1,cc};
            if isempty(serve_user)
                AP_Null_count = AP_Null_count + 1;
                Cluster_Null_id = [Cluster_Null_id cc];
                continue
            end
            CAP = CAP_id(cc);
            AP_id = CAP_AP_id{1,cc};
            if isempty(AP_id)
                AP_Null_count = AP_Null_count + 1;
                Cluster_Null_id = [Cluster_Null_id cc];
                H_cluster_user{1,cc} = cell2mat(H_DL(CAP,serve_user));
            else
                H_cluster_user{1,cc} = cell2mat(H_DL(AP_id,serve_user));
            end
            V_RZF_cluster{1,cc} = H_cluster_user{1,cc}*pinv(H_cluster_user{1,cc}'*H_cluster_user{1,cc}+Regular*eye(length(serve_user)));
            V_RZF_cluster{1,cc} = V_RZF_cluster{1,cc} ./ vecnorm(V_RZF_cluster{1,cc});

        end

        DL_SE(itee) = Function_MDD_THz_sub1(H_DL,H_cluster_user,V_RZF_cluster,Cluster_id,G_u,U_l,Number_AP,N_ap,Power_AP_W,No_6G,Sub6G_Bandwidth);

        Max_sub = NumberOfSub*(1/2);

        CC_CA_Num_ite = 5;


        for CC_CA_ite = 1:CC_CA_Num_ite
            %% Sub 2
            %% Initialization M_CC & M_CA

            M_FH = 1:NumberOfSub;
            M_CA_MDD = [];

            CA_f_gain = zeros(NumberOfSub,Number_cluster);
            for cc = 1:Number_cluster %% acquire the channel gain of sub
                CAP = CAP_id(cc);
                AP_id = Cluster_id{1,cc}(2:end);
                temp = cell2mat(H_CA_f(CAP,AP_id));
                temp = diag(temp' * temp);
                CA_f_gain(:,cc) = sum(reshape(temp,NumberOfSub,[]),2);
            end
            CA_f_gain_temp = CA_f_gain;
            while length(M_CA_MDD)<Max_sub
                if Number_cluster==Number_AP
                    break
                end
                for cc = 1:Number_cluster
                    [value_temp,Sub_id] = max(CA_f_gain_temp(:,cc));
                    if value_temp==0 %%cluster may have only one AP
                        continue
                    else
                        CA_f_gain_temp(Sub_id,cc) = 0;
                    end
                    M_CA_MDD = union(M_CA_MDD,Sub_id);
                    if length(M_CA_MDD)==Max_sub
                        break
                    end
                end
            end
            M_CC_MDD = setdiff(M_FH,M_CA_MDD);
            %%% CPU-to-CAP
            H_CPU_CAP = cell(1,NumberOfSub);
            F_RZF = cell(1,NumberOfSub);
            CC_channel_gain = cell(1,NumberOfSub);
            for m = 1:length(M_FH)
                for ll = 1:Number_cluster
                    H_CPU_CAP{m} = [H_CPU_CAP{m} H_CC_f{CAP_id(ll)}(:,m)];
                end
                F_RZF{m} = H_CPU_CAP{m}*pinv(H_CPU_CAP{m}'*H_CPU_CAP{m} + Regular*eye(Number_cluster));
                F_RZF{m} = F_RZF{m} ./ vecnorm(F_RZF{m});
                CC_channel_gain{m} = H_CPU_CAP{m}' * F_RZF{m};
            end
            CC_channel_gain_MDD = cell(1,NumberOfSub);
            CC_channel_gain_MDD(1,M_CC_MDD) = CC_channel_gain(1,M_CC_MDD);
            CC_temp(itee,CC_CA_ite) = Function_MDD_THz_sub2_CC(M_CC_MDD,CAP_id,CC_channel_gain_MDD,G_u,U_l,Power_CPU_W,No_THz,SI_W,SubTHz_Bandwidth);

            %%% CAP-to-AP
            Real_cluster_id = setdiff(1:Number_cluster,Cluster_Null_id);
            if isempty(Real_cluster_id)
                Fronthaul_CA_MDD(itee) = 0;
            else
                cluster_count = Number_cluster - AP_Null_count;
                H_CAP_AP = cell(Number_cluster,NumberOfSub);
                W_RZF = cell(Number_cluster,NumberOfSub);
                CA_channel_gain = cell(Number_cluster,NumberOfSub);
                CA_channel_gain_MDD = cell(Number_cluster,NumberOfSub);
                for ll = 1:cluster_count
                    clusterId = Real_cluster_id(ll);
                    CAPId = CAP_id(clusterId);
                    AP_pool = CAP_AP_id{1,clusterId};
                    for m = 1:length(M_FH)
                        for lll = 1:length(AP_pool)
                            H_CAP_AP{clusterId,m} = [H_CAP_AP{clusterId,m} H_CA_f{CAPId,AP_pool(lll)}(:,m)];
                        end
                        W_RZF{clusterId,m} = H_CAP_AP{clusterId,m}*pinv(H_CAP_AP{clusterId,m}'*H_CAP_AP{clusterId,m} + Regular*eye(length(AP_pool)));
                        W_RZF{clusterId,m} = W_RZF{clusterId,m} ./ vecnorm(W_RZF{clusterId,m});
                        CA_channel_gain{clusterId,m} = H_CAP_AP{clusterId,m}' * W_RZF{clusterId,m};
                    end
                    CA_channel_gain_MDD(clusterId,M_CA_MDD) = CA_channel_gain(clusterId,M_CA_MDD);
                end

                %%%% Inter-cluster channel gain
                Inter_cluster_ChannelGain = cell(Number_AP,NumberOfSub);
                Inter_cluster_ChannelGain_MDD = cell(Number_AP,NumberOfSub);
                for c = 1:cluster_count
                    clusterId = Real_cluster_id(c);
                    AP_pool = CAP_AP_id{1,clusterId};
                    inter_clusterId = setdiff(Real_cluster_id,clusterId);
                    for ll = 1:length(AP_pool)
                        APId = AP_pool(ll);
                        for m = 1:length(M_FH)
                            for lll = 1:length(inter_clusterId)
                                Id = CAP_id(inter_clusterId(lll));
                                channel = H_CA_f{Id,APId}(:,m);
                                gain = channel' * W_RZF{inter_clusterId(lll),m};
                                Inter_cluster_ChannelGain{APId,m} = ...
                                    [Inter_cluster_ChannelGain{APId,m} gain];
                            end
                        end
                        Inter_cluster_ChannelGain_MDD(APId,M_CA_MDD) = Inter_cluster_ChannelGain(APId,M_CA_MDD);
                    end
                end


                CA_temp(itee,CC_CA_ite) = Function_MDD_THz_sub2_CA(M_CA_MDD,CAP_id,CAP_AP_id,Real_cluster_id,Inter_cluster_ChannelGain_MDD,...
                    CA_channel_gain_MDD,G_u,U_l,Number_AP,Power_AP_W,No_THz,SubTHz_Bandwidth);
            end

            if CC_temp(itee,CC_CA_ite) > CA_temp(itee,CC_CA_ite)
                Max_sub = Max_sub + 8*0.5^(CC_CA_ite-1);
            else
                Max_sub = Max_sub - 8*0.5^(CC_CA_ite-1);
            end
            CC_CA_abs(itee,CC_CA_ite) = min(CC_temp(itee,CC_CA_ite),CA_temp(itee,CC_CA_ite));
            if abs(CC_temp(itee,CC_CA_ite)*1e-9 - CA_temp(itee,CC_CA_ite)*1e-9) <= 0.1
                break
            end
        end
        [~,CC_CA_index] = max(CC_CA_abs(itee,:));
        Fronthaul_CC_MDD(itee) = CC_temp(itee,CC_CA_index);
        Fronthaul_CA_MDD(itee) = CA_temp(itee,CC_CA_index);
        Final_MDD(itee) = min([Fronthaul_CC_MDD(itee),Fronthaul_CA_MDD(itee), DL_SE(itee)]);

    end

    Par_MDD(par_iter,:) = max(Final_MDD);
    Par_MDD_CC(par_iter,:) = Fronthaul_CC_MDD;
    Par_MDD_CA(par_iter,:) = Fronthaul_CA_MDD;
    Par_DL_SE(par_iter,:) = DL_SE;

end
if Regular == 0
    temp_name = 1;
else
    temp_name = log10(Regular);
end
Bandwidth_name = Bandwidth*1e-9;

filename = ['RE20_HzMDDVs_',num2str(Number_AP),num2str(Number_User),'_',num2str(Bandwidth_name),'_','UC',num2str(U_max),num2str(C_max),'.txt'];
filename = sprintf(filename);
fid1 = fopen(filename,'at+');
fprintf(fid1,'1 ');
fprintf(fid1,'%6.6f ',Par_MDD_CC);
fprintf(fid1,'\n');
fprintf(fid1,'2 ');
fprintf(fid1,'%6.6f ',Par_MDD_CA);
fprintf(fid1,'\n');
fprintf(fid1,'5 ');
fprintf(fid1,'%6.6f ',Par_DL_SE);
fprintf(fid1,'\n');
fprintf(fid1,'6 ');
fprintf(fid1,'%6.6f ',Par_MDD);
fprintf(fid1,'\n');
fclose(fid1);
end