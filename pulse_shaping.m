function [p_rc] = pulse_shaping(S_t,beta,Ts)

if S_t==Ts/(2*beta)
    p_rc = (pi/4) * sinc(1/(2*beta));
else
    p_rc = sinc(S_t/Ts)*cos(pi*beta*S_t/Ts)/(1-(2*beta*S_t/Ts)^2);



end

