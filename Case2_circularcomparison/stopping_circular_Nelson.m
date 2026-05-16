function [value,isterminal,direction] = stopping_circular_Nelson(t,y)
global rd
value(1)      = y(1) - (rd+.2); 
isterminal(1) = 1; % stop the integration(once the condition is met stop the integration)
direction(1)  = 0; % negative direction(as R decreases from positive to zero d=-1;
                  % If R increases from negative to zero d=+1;d=0 implies no need of direction )
end