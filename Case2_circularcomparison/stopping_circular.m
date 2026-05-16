function [value,isterminal,direction] = stopping_circular(t,x)
global rd
value(1)      = ((x(1) - (rd+.1)) && (x(7) - (rd+.1))) ; 
isterminal(1) = 1; % stop the integration(once the condition is met stop the integration)
direction(1)  = 0; % negative direction(as R decreases from positive to zero d=-1;
                  % If R increases from negative to zero d=+1;d=0 implies no need of direction )
% value(2)      = x(7) - (rd+.01); 
% isterminal(2) = 1; % stop the integration(once the condition is met stop the integration)
% direction(2)  = 0; % negative direction(as R decreases from positive to zero d=-1;
%                   % If R increases from negative to zero d=+1;d=0 implies no need of direction )
end