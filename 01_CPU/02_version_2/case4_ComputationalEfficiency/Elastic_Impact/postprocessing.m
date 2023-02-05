clear,close
set(0,'defaulttextinterpreter','latex')                                   ;%
fslab  = 14; fsleg  = 14; fstit  = 14; fstick = 14;
nstr = 4;

for i = 1:8
D_time = load(['.\data\CPDI2q_time_vectorized_' num2str(i) '.mat']);
D_mpDmeD = load(['.\data\data_' num2str(i) '.mat']);

it_per_sec(i,:) = D_time.ItPs
np(i)           = D_mpDmeD.mpD.n

end

loglog(np,it_per_sec(:,1),'o--')


Julia = [132.80 33.37 26.45 1.82]
MATLAB= it_per_sec([3 4 5 8])
gain  = MATLAB./Julia


