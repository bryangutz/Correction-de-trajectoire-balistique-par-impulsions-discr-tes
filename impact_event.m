function [value, isterminal, direction] = impact_event(~, x)
    value = x(3); % z=0
    isterminal = 1; % stop integration
    direction = +1; % detect change of sign from - to +
end