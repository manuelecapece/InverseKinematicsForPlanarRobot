function [A10,A20,A30,A40] = untitled(a,q)
%Funzione che calcola la matrice di trasformazione omogena T04 data una
%determinata configurazione dello spazio dei giunti q

A10= [cos(q(1)),-sin(q(1)),0,a*cos(q(1));sin(q(1)),cos(q(1)),0,a*sin(q(1));0,0,1,0;0,0,0,1];
A21= [cos(q(2)),-sin(q(2)),0,a*cos(q(2));sin(q(2)),cos(q(2)),0,a*sin(q(2));0,0,1,0;0,0,0,1];
A20 = A10*A21;
A32= [cos(q(3)),0,sin(q(3)),0;sin(q(3)),0,-cos(q(3)),0;0,1,0,0;0,0,0,1];
A30 = A20*A32;
A43= [1,0,0,0;0,1,0,0;0,0,1,q(4);0,0,0,1];
A40 = A30*A43

end

