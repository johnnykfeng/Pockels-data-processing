% E_field_integration for each pixel column
E=Output.E_field_biased_corrected;
I0=Calib.I_parallel_nobias_field;%Output.I_parallel_nobias_field;
dimention=size(E);
N=dimention(1);
L=dimention(2);
% find edge
for i=1:N
    Iline=I0(i,:);
    I_signal_rough=Iline(Iline>200);
    Signal_median=median(I_signal_rough);
    Th=Signal_median/2;%threshold for judging a edge
    signal_range=find(Iline>Th);
    Edge_left(i)=signal_range(1);
    Edge_right(i)=signal_range(end);
end
% find outliers
Edge_left_median=median(Edge_left);
Edge_right_median=median(Edge_right);
T=Edge_right-Edge_left;% thickness of sensor in unite of pixel
T_median=median(T);
outlier=find((abs(T-T_median)/T_median>0.2)==1);
for i=1:N
    Eline=E(i,Edge_left(i):Edge_right(i));
    x=(0:T(i))/T(i)*2E-3;
    E_int(i)=trapz(x,Eline);
end
figure
plot(E_int)