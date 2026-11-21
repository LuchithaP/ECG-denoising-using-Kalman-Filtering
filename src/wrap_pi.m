function y = wrap_pi(x)
%WRAP_PI  Wrap angles to [-pi, pi).
    y = mod(x + pi, 2*pi) - pi;
end
