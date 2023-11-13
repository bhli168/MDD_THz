function SE = Function_MDD_THz_sub1(H_DL,H_cluster_user,V_RZF_cluster,Cluster_id,G_u,U_l,Number_AP,N_ap,Power_AP_W,No,Sub6G_Bandwidth)

Number_User = length(G_u);
Ini_SINR = zeros(1,Number_User);
AP_user_set = cell(1,Number_User);
AP_user_V = cell(1,Number_User);
H_AP_user = cell(1,Number_User);
for uu = 1:Number_User
    cluster_index = G_u{1,uu};
    for ll = 1:length(cluster_index)
        AP_index = Cluster_id{1,cluster_index(ll)};
        if length(AP_index)==1
            AP_user_set{uu} = [AP_user_set{uu} AP_index];
        else
            AP_user_set{uu} = [AP_user_set{uu} AP_index(2:end)];
        end
        u_index = find(U_l{cluster_index(ll)}==uu);
        AP_user_V{uu} = [AP_user_V{uu};V_RZF_cluster{1,cluster_index(ll)}(:,u_index)]; %%construct RZF precoders from serving APs to each user
    end
    H_AP_user{uu} = cell2mat(H_DL(AP_user_set{uu},uu));
end
%% find upper bound
for uu = 1:Number_User
    AP_serve_index = AP_user_set{1,uu};
    AP_serve_len = length(AP_serve_index);

    %cvx_solver Mosek
    cvx_begin quiet
    variable p_all(Number_AP,Number_User)
    expression SINR
    expression p_opt
    expression p_opt_tri
    %expression gain
    p_opt = p_all(AP_serve_index,uu);
    p_opt_tri = kron(diag(p_opt),eye(N_ap));
    %     for ll = 1:AP_serve_len
    %         SINR(1,ll) = pow_p(norm(H_AP_user{uu}(N_ap*(ll-1)+1:N_ap*ll)' * AP_user_V{uu}(N_ap*(ll-1)+1:N_ap*ll)),2);
    %     end
    %gain = SINR * p_opt / No;
    SINR = real(H_AP_user{uu}' * p_opt_tri * AP_user_V{uu}) / sqrt(No);
    maximize SINR

    subject to
    for ll = 1:AP_serve_len
        pow_p(p_opt(ll) * norm(AP_user_V{uu}(N_ap*(ll-1)+1:N_ap*ll)),2) <= Power_AP_W;
        p_opt(ll) >=0;

    end
    cvx_end

    Ini_SINR(uu) = SINR;
end
p_ini = zeros(Number_AP,Number_User);
for ll = 1:length(U_l)
    user_id = U_l{1,ll};
    ap_id = Cluster_id{1,ll};
    if isempty(user_id)
        continue
    else
        p_ini(ap_id,user_id) = Power_AP_W/length(user_id);
    end
end
p_all_final = p_ini;
[upper,~] = min(Ini_SINR.^2);
lower = 0;
delta = 1;
while upper - lower > delta
    mid = (upper+lower) / 2;

    cvx_begin quiet
    variable p_all(Number_AP,Number_User) nonnegative
    expression I_u(Number_User,Number_User)
    expression S_u(Number_User,1)
    for uu = 1:Number_User
        AP_serve_index = AP_user_set{1,uu};
        AP_serve_len = length(AP_serve_index);
        p_opt_tri = kron(diag(p_all(AP_serve_index,uu)),eye(N_ap));
        S_u(uu) = real(H_AP_user{uu}' * p_opt_tri * AP_user_V{uu});
        %     for ll = 1:AP_serve_len
        %         temp = AP_serve_index(ll);
        %         S_u(uu) =  S_u(uu) + p_all(temp,uu) * (H_AP_user{uu}(N_ap*(ll-1)+1:N_ap*ll)' * AP_user_V{uu}(N_ap*(ll-1)+1:N_ap*ll)); %% p = sqrt(p~)
        %     end
        for uuu = 1:Number_User
            if uu == uuu
                I_u(uu,uuu) = 0;
                continue;
            else
                Inter_AP = AP_user_set{1,uuu};
                H_inter = cell2mat(H_DL(Inter_AP,uu));
                p_inter_opt = kron(diag(p_all(Inter_AP,uuu)),eye(N_ap));
                I_u(uu,uuu) = real(H_inter' * p_inter_opt * AP_user_V{uuu});
            end
        end
    end

    minimize norm(sum(p_all,1))

    subject to
    for uu = 1:Number_User
        S_u(uu) * sqrt(1/mid) >= norm([I_u(uu,:) sqrt(No)]);
    end
    for ll = 1:length(Cluster_id)
        if length(Cluster_id{1,ll})==1
            AP_id = Cluster_id{1,ll};
        else
            AP_id = Cluster_id{1,ll}(2:end);
        end
        User_id = U_l{1,ll};
        if isempty(User_id)
            continue
        end
        for lll = 1:length(AP_id)
            norm(norms(V_RZF_cluster{ll}(N_ap*(lll-1)+1:N_ap*lll,:))...
                * diag(p_all(AP_id(lll),User_id))) <= sqrt(Power_AP_W);
        end

    end

    cvx_end
    if cvx_status(1)=='S'
        lower = mid;
        p_all_final = p_all;
    else
        upper = mid;
    end

end

%% Compute real rate
S_final = zeros(1,Number_User);
I_final = zeros(Number_User,Number_User);
for uu = 1:Number_User
    AP_serve_index = AP_user_set{1,uu};
    p_final_tri = kron(diag(p_all_final(AP_serve_index,uu)),eye(N_ap));
    S_final(uu) = real(H_AP_user{uu}' * p_final_tri * AP_user_V{uu});
    for uuu = 1:Number_User
        if uu == uuu
            I_final(uu,uuu) = 0;
            continue;
        else
            Inter_AP = AP_user_set{1,uuu};
            H_inter = cell2mat(H_DL(Inter_AP,uu));
            p_inter_opt = kron(diag(p_all_final(Inter_AP,uuu)),eye(N_ap));
            I_final(uu,uuu) = real(H_inter' * p_inter_opt * AP_user_V{uuu});
        end
    end
end
for uu = 1:Number_User
    SNR(uu) = abs(S_final(uu))^2 / norm([I_final(uu,:) sqrt(No)])^2;
end
SE = Sub6G_Bandwidth * log2(1+min(SNR));


end