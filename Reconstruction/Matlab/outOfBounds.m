function [outOfBounds] = outOfBounds(x, y)

right = bool(x > MAX_X);

left = bool(x < MIN_X);

top = bool(y > MAX_Y);

bottom = bool(y < MIN_Y);

outOfBounds = bool(right && left && top && bottom);

return;