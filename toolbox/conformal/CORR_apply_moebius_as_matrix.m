function [T_Z] = CORR_apply_mobius_as_matrix(Mob,Z)

%% Z is Nx1 complex vector

T_Z = (Mob(1,1)*Z+Mob(1,2))./(Mob(2,1)*Z+Mob(2,2));