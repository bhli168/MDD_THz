function SE = Function_MDD_THz_sub2_CA(M_CA,CAP_id,CAP_AP_id,Real_cluster_id,Inter_cluster_ChannelGain,...
    CA_channel_gain,G_u,U_l,N_AP,Power_AP_W,No,SubTHz_Bandwidth)


N_CAP = length(find(CAP_id>0));
N_User = length(G_u);
N_sub = length(M_CA);
Psi = 1e2;
%thre = Power_CPU_W/(N_sub*N_User*N_CAP)/10;

%% Initialize p
for c = 1:length(CAP_AP_id)
    maxx(c) = length(CAP_AP_id{1,c});
end
p_ini = (Power_AP_W/(N_sub*max(maxx))) .* ones(N_AP,N_User,N_sub);
% Initialize \gamma
gamma_ini = zeros(N_AP,N_User,N_sub);
for c = 1:length(Real_cluster_id)
    clusterId = Real_cluster_id(c);
    AP_id = CAP_AP_id{1,clusterId};
    ser_id = U_l{clusterId};
    subpool = buffer(1:N_sub,length(ser_id));
    for ll = 1:length(AP_id)
        APId = AP_id(ll);
        for u = 1:length(ser_id)
            serId = ser_id(u);
            temp = subpool(u,find(subpool(u,:)>0));
            gamma_ini(APId,serId,temp) = 1;
        end
    end
end
p_ini = p_ini .* gamma_ini;
p_lo_ini = p_ini .* gamma_ini;
%% Compute z
z_ini = zeros(N_AP,N_User,N_sub);
for c = 1:length(Real_cluster_id)
    clusterId = Real_cluster_id(c);
    AP_id = CAP_AP_id{1,clusterId};
    ser_id = U_l{clusterId};
    inter_clu = setdiff(Real_cluster_id,clusterId);
    for m = 1:N_sub
        mId = M_CA(m);
        gain_temp = CA_channel_gain{clusterId,mId};
        for ll = 1:length(AP_id)
            APId = AP_id(ll);
            gain_A = abs(gain_temp(ll,ll)).^2;
            gain_B = abs(Inter_cluster_ChannelGain{APId,mId}).^2;
            for u = 1:length(ser_id)
                serId = ser_id(u);
                if p_ini(APId,serId,m)==0
                    continue
                else

                    inter_Id = cell2mat(CAP_AP_id(1,inter_clu));
                    if isempty(inter_clu)
                        z_ini(APId,serId,m) = sqrt(p_ini(APId,serId,m)*gain_A*1e12)...
                            / (1e6*No);
                    else

                        z_ini(APId,serId,m) = sqrt(p_ini(APId,serId,m)*gain_A*1e12)...
                            / (1e6*(gain_B*sum(p_ini(inter_Id,:,m),2)+No));
                    end

                end
            end
        end
    end
end


%% find upper bound
cvx_begin quiet
variable p(1,N_sub) nonnegative
expression R(1,N_sub)
for m = 1:N_sub
    R(m) = log(1 + p(1,m)*abs(CA_channel_gain{Real_cluster_id(1),M_CA(m)}(1,1))^2/No)/log(2);
end

maximize sum(R)

subject to
sum(p)<=Power_AP_W;
cvx_end

upper = sum(R);
lower = 0;
delta = 1;
iter_count = 0;

p_final = p_ini;

while 1
    while upper - lower > delta
        mid = (upper+lower) / 2;
        %%cvx_solver SeDuMi
        cvx_begin quiet
        %cvx_precision low
        variable p(N_AP,N_User,N_sub) nonnegative
        expression A(N_AP,N_User,N_sub)
        expression B(N_AP,N_User,N_sub)
        expression C(N_AP,N_User,N_sub)
        expression D(N_AP,N_User,N_sub)
        expression Lo(N_AP,N_sub)
        for c = 1:length(Real_cluster_id)
            clusterId = Real_cluster_id(c);
            AP_id = CAP_AP_id{1,clusterId};
            ser_id = U_l{clusterId};
            inter_clu = setdiff(Real_cluster_id,clusterId);
            for m = 1:N_sub
                mId = M_CA(m);
                gain_temp = CA_channel_gain{clusterId,mId};
                for ll = 1:length(AP_id)
                    APId = AP_id(ll);
                    gain_A = abs(gain_temp(ll,ll)).^2;
                    gain_B = abs(Inter_cluster_ChannelGain{APId,mId}).^2;
                    for u = 1:length(ser_id)
                        serId = ser_id(u);
                        inter_Id = cell2mat(CAP_AP_id(1,inter_clu));
                        A(APId,serId,m) = p(APId,serId,m)*gain_A;
                        if isempty(inter_clu)
                            B(APId,serId,m) = No;
                        else
                            B(APId,serId,m) = gain_B*sum(p(inter_Id,:,m),2)+No;
                        end
                        C(APId,serId,m) = 1*1e6 + 2*z_ini(APId,serId,m)*sqrt(A(APId,serId,m)*1e12)-...
                            z_ini(APId,serId,m)^2*B(APId,serId,m)*1e6;
                        D(APId,serId,m) = log(C(APId,serId,m)*1e-6) / log(2);
                        %                     Lo(APId,m) = Lo(APId,m) + (1-exp(-Psi*p_lo_ini(APId,serId,m)))...
                        %                         + Psi*exp(-Psi*p_lo_ini(APId,serId,m))*(p(APId,serId,m)-p_lo_ini(APId,serId,m));

                    end
                end
            end
        end

        minimize sum(p(:))
        subject to

        for c = 1:length(Real_cluster_id)
            clusterId = Real_cluster_id(c);
            AP_id = CAP_AP_id{1,clusterId};
            ser_id = U_l{clusterId};
            for ll = 1:length(AP_id)
                APId = AP_id(ll);
                for u = 1:length(ser_id)
                    serId = ser_id(u);
                    sum(D(APId,serId,:)) >= mid;
                end
            end
        end
        for c = 1:length(Real_cluster_id)
            clusterId = Real_cluster_id(c);
            AP_id = CAP_AP_id{1,clusterId};

            sum(sum(sum(p(AP_id,:,:)))) <= Power_AP_W;

        end
        %     for c = 1:length(Real_cluster_id)
        %         clusterId = Real_cluster_id(c);
        %         AP_id = CAP_AP_id{1,clusterId};
        %         for m = 1:N_sub
        %             for ll = 1:length(AP_id)
        %                 APId = AP_id(ll);
        %                 Lo(APId,m) <=1;
        %             end
        %         end
        %     end

        cvx_end
        if cvx_status(1)=='S'
            lower = mid;
            p_final = p;
            p_lo_ini = p;
            iter_count = iter_count + 1;
            if iter_count>10
                break
            end
            for c = 1:length(Real_cluster_id)
                clusterId = Real_cluster_id(c);
                AP_id = CAP_AP_id{1,clusterId};
                ser_id = U_l{clusterId};
                inter_clu = setdiff(Real_cluster_id,clusterId);
                for m = 1:N_sub
                    mId = M_CA(m);
                    gain_temp = CA_channel_gain{clusterId,mId};
                    for ll = 1:length(AP_id)
                        APId = AP_id(ll);
                        gain_A = abs(gain_temp(ll,ll)).^2;
                        gain_B = abs(Inter_cluster_ChannelGain{APId,mId}).^2;
                        for u = 1:length(ser_id)
                            serId = ser_id(u);
                            if p(APId,serId,m)<1e-5
                                continue
                            else

                                inter_Id = cell2mat(CAP_AP_id(1,inter_clu));
                                if isempty(inter_clu)
                                    z_ini(APId,serId,m) = sqrt(p(APId,serId,m)*gain_A*1e12)...
                                        / (1e6*No);
                                else
                                    z_ini(APId,serId,m) = sqrt(p(APId,serId,m)*gain_A*1e12)...
                                        / (1e6*(gain_B*sum(p(inter_Id,:,m),2)+No));
                                end

                            end
                        end
                    end
                end
            end
            CAP_count = length(Real_cluster_id);
            if upper-lower<=delta && CAP_count*Power_AP_W-sum(p(:))>1
                upper = upper+5;
            end

        else
            upper = mid;
        end

    end
    if isequal(p_ini,p_final) && delta>0.1
        delta = delta*0.5;
        upper = sum(R)*delta;
        lower = 0;
    else
        break
    end
end



SINR = zeros(N_AP,N_User,N_sub);
for c = 1:length(Real_cluster_id)
    clusterId = Real_cluster_id(c);
    AP_id = CAP_AP_id{1,clusterId};
    ser_id = U_l{clusterId};
    inter_clu = setdiff(Real_cluster_id,clusterId);
    for m = 1:N_sub
        mId = M_CA(m);
        gain_temp = CA_channel_gain{clusterId,mId};
        for ll = 1:length(AP_id)
            APId = AP_id(ll);
            gain_A = abs(gain_temp(ll,ll)).^2;
            gain_B = abs(Inter_cluster_ChannelGain{APId,mId}).^2;
            for u = 1:length(ser_id)
                serId = ser_id(u);
                inter_Id = cell2mat(CAP_AP_id(1,inter_clu));
                if isempty(inter_clu)
                    SINR(APId,serId,m) = log2(1 + p_final(APId,serId,m)*gain_A...
                        / No);
                else
                    SINR(APId,serId,m) = log2(1 + p_final(APId,serId,m)*gain_A...
                        / (gain_B*sum(p_final(inter_Id,:,m),2)+No));
                end

            end
        end
    end
end
temp = sum(SINR,3);
temp = temp(find(temp>0));
SE = SubTHz_Bandwidth * min(temp);

end
