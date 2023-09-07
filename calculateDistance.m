
function distance = calculateDistance(coord1, coord2)
% Calculate the distance between two 3-D coordinates
% coord1 and coord2 should be arrays of size 1x3

% Calculate the differences in each dimension
dx = coord2(1) - coord1(1);
dy = coord2(2) - coord1(2);
dz = coord2(3) - coord1(3);

% Calculate the distance using the Pythagorean theorem
distance = sqrt(dx^2 + dy^2 + dz^2);
end