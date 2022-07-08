%%
% 该函数可以将电流指定连续单元的电流强度赋值为0
% variables:
% Amp: 馈电电流幅度
%%
function Amp = makeImZero(Amp,istart,iend,jstart,jend)
for i = istart:1:iend
    for j = jstart:1:jend
        Amp(i,j)=0;
    end
end

