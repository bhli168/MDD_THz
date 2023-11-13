function [P_location,P_height] = PLsetup(Num_P,len,wid,sys)

%% AP MS distribution
while 1
    P_location = len * rand(Num_P,1) + 1i * wid * rand(Num_P,1);

    P_P_PL = abs(repmat(P_location,1,Num_P) - repmat(P_location.',Num_P,1));
    a1 = length(find(P_P_PL>0&P_P_PL<1));

    if a1==0 || strcmp('user',sys)
        break
    end
end

P_height = (6-4) .* rand(Num_P,1) + 4;

end