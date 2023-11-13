function SE = Function_MDD_THz_sub2_CC(M_CC,CAP_id,CC_channel_gain,G_u,U_l,Power_CPU_W,No,SI_W,SubTHz_Bandwidth)

N_CAP = length(find(CAP_id>0));
N_User = length(G_u);
N_sub = length(M_CC);
Psi = 1e2;

AP_Null_count = 0;
Cluster_Null_id = [];
for cc = 1:N_CAP
    serve_user = U_l{1,cc};
    if isempty(serve_user)
        AP_Null_count = AP_Null_count + 1;
        Cluster_Null_id = [Cluster_Null_id cc];
    end
end
Real_cluster_id = setdiff(1:N_CAP,Cluster_Null_id);
N_CAP_real = length(Real_cluster_id);

%thre = Power_CPU_W/(N_sub*N_User*N_CAP)/10;

%% Initialize p
p_ini = (Power_CPU_W/(N_sub*N_CAP_real)) .* ones(N_CAP,N_User,N_sub);
% Initialize \gamma
gamma_ini = zeros(N_CAP,N_User,N_sub);
% for ll = 1:N_CAP
%     for m = 1:N_sub
%         ran = randi(length(U_l{ll}));
%         gamma_ini(ll,U_l{ll}(ran),m) = 1;
%     end
% end
for ll = Real_cluster_id
    ser = length(U_l{ll});
    subpool = buffer(1:N_sub,ser);
    for u = 1:ser
        temp = subpool(u,find(subpool(u,:)>0));
        gamma_ini(ll,U_l{ll}(u),temp) = 1;
    end
end
p_ini = p_ini .* gamma_ini;
p_lo_ini = p_ini .* gamma_ini;
%% Compute z
z_ini = zeros(N_CAP,N_User,N_sub);
for m = 1:N_sub
    gain_temp = CC_channel_gain{M_CC(m)};
    for ll = Real_cluster_id
        gain_temp_1 = gain_temp(ll,:);
        gain_A = abs(gain_temp_1(ll)).^2;
        inter_id = setdiff(Real_cluster_id,ll);
        gain_B = abs(gain_temp_1(inter_id)).^2;
        inter_power = sum(p_ini(:,:,m),2);
        inter_power = inter_power(inter_id);
        inter = gain_B * inter_power;
        for u = 1:length(U_l{ll})
            if p_ini(ll,U_l{ll}(u),m)==0
                continue
            else
                z_ini(ll,U_l{ll}(u),m) = sqrt(p_ini(ll,U_l{ll}(u),m) * gain_A)...
                    / (inter+No+SI_W);
            end
        end
    end
end

%% find upper bound
cvx_begin quiet
variable p(1,N_sub) nonnegative
expression R(1,N_sub)
for m = 1:N_sub
    R(m) = log(1 + p(1,m)*abs(CC_channel_gain{M_CC(m)}(1,1))^2/(No+SI_W))/log(2);
end

maximize sum(R)

subject to
sum(p)<=Power_CPU_W;
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
        variable p(N_CAP,N_User,N_sub) nonnegative
        expression A(N_CAP,N_User,N_sub)
        expression B(N_CAP,N_User,N_sub)
        expression C(N_CAP,N_User,N_sub)
        expression D(N_CAP,N_User,N_sub)
        expression Lo(N_CAP,N_sub)
        for m = 1:N_sub
            gain_temp = CC_channel_gain{M_CC(m)};
            for ll = Real_cluster_id
                gain_temp_1 = gain_temp(ll,:);
                gain_A = abs(gain_temp_1(ll)).^2;
                inter_id = setdiff(Real_cluster_id,ll);
                gain_B = abs(gain_temp_1(inter_id)).^2;
                inter_power = sum(p(:,:,m),2);
                inter_power = inter_power(inter_id);
                for u = 1:length(U_l{ll})
                    B(ll,U_l{ll}(u),m) = gain_B * inter_power+No+SI_W; %% problem?
                    A(ll,U_l{ll}(u),m)= p(ll,U_l{ll}(u),m) * gain_A;
                    C(ll,U_l{ll}(u),m) = 1*1e6 + 2*z_ini(ll,U_l{ll}(u),m)*sqrt(A(ll,U_l{ll}(u),m)*1e12)-...
                        z_ini(ll,U_l{ll}(u),m)^2*B(ll,U_l{ll}(u),m)*1e6;
                    D(ll,U_l{ll}(u),m) = log(C(ll,U_l{ll}(u),m)*1e-6)/log(2);
                    %                 Lo(ll,m) = Lo(ll,m) + (1-exp(-Psi*p_lo_ini(ll,U_l{ll}(u),m)))...
                    %                     + Psi*exp(-Psi*p_lo_ini(ll,U_l{ll}(u),m))*(p(ll,U_l{ll}(u),m)-p_lo_ini(ll,U_l{ll}(u),m));
                end
            end
        end

        minimize sum(p(:))
        subject to

        for ll = Real_cluster_id
            for u = 1:length(U_l{ll})
                sum(D(ll,U_l{ll}(u),:)) >= mid;
            end
        end
        %     for ll = 1:N_CAP
        %         for m = 1:N_sub
        %             Lo(ll,m) <=1;
        %         end
        %     end

        sum(p(:)) <= Power_CPU_W;

        cvx_end
        if cvx_status(1)=='S'
            lower = mid;
            p_final = p;
            A_final = A;
            B_final = B;
            p_lo_ini = p;
            iter_count = iter_count + 1;
            if iter_count>10
                break
            end
            for m = 1:N_sub
                gain_temp = CC_channel_gain{M_CC(m)};
                for ll = 1:Real_cluster_id
                    gain_temp_1 = gain_temp(ll,:);
                    gain_A = abs(gain_temp_1(ll)).^2;
                    inter_id = setdiff(Real_cluster_id,ll);
                    gain_B = abs(gain_temp_1(inter_id)).^2;
                    inter_power = sum(p(:,:,m),2);
                    inter_power = inter_power(inter_id);
                    inter = gain_B * inter_power;
                    for u = 1:length(U_l{ll})
                        if p(ll,U_l{ll}(u),m)<1e-5
                            continue
                        else
                            z_ini(ll,U_l{ll}(u),m) = sqrt(p(ll,U_l{ll}(u),m) * gain_A)...
                                / (inter+No+SI_W);
                        end
                    end
                end
            end
            if upper-lower<=delta && Power_CPU_W-sum(p(:))>1
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

SINR = zeros(N_CAP,N_User,N_sub);
for m = 1:N_sub
    gain_temp = CC_channel_gain{M_CC(m)};
    for ll = Real_cluster_id
        gain_temp_1 = gain_temp(ll,:);
        gain_A = abs(gain_temp_1(ll)).^2;
        inter_id = setdiff(Real_cluster_id,ll);
        gain_B = abs(gain_temp_1(inter_id)).^2;
        inter_power = sum(p_final(:,:,m),2);
        inter_power = inter_power(inter_id);
        for u = 1:length(U_l{ll})
            interference = gain_B * inter_power+No+SI_W; %% problem?
            signal= p_final(ll,U_l{ll}(u),m) * gain_A;
            SINR(ll,U_l{ll}(u),m) = log2(1 + signal / interference);
        end
    end
end
temp = sum(SINR,3);
temp = temp(find(temp>0));
SE = SubTHz_Bandwidth * min(temp);

end
