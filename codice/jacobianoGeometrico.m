function [J] = jacobianoGeometrico(q,A10,A20,A30,A40)
p0 = [0;0;0];
z0 = [0;0;1];

p1 =A10(1:3,end);
z1 =A10(1:3,3);

p2 =A20(1:3,end);
z2 =A20(1:3,3);

z3 =A30(1:3,3);

p4 =A40(1:3,end);

J1 = [cross(z0,(p4-p0)); z0];
J2 = [cross(z1,(p4-p1)); z1];
J3 = [cross(z2,(p4-p2)); z2];
J4 = [z3; p0];

J = [J1, J2, J3, J4];
end

