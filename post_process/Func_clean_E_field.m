% clean up the electric field map by removing the region out of the sensor
% and the outliers inside the sensor
function [E_clean,Edge_cathode,Edge_anode]=Func_clean_E_field(Output,Calib)
%% define signal range
% Threth_I0=200;
version=str2double(Output.version);
if version<2
    I0=Output.I_parallel_nobias_field;
    I0_cross=Output.I_crossed_nobias_field;
else
    I0=Calib.I_parallel_nobias_field;
    I0_cross=Calib.I_crossed_nobias_field;
end

dimention=size(Output.I_crossed_biased_field);
N=dimention(1);
L=dimention(2);
%% define threthold of signal in I0
[counts,edges] = histcounts(I0,50);
I0_peak=edges(counts==max(counts(10:end)));
Threth_I0=I0_peak/4;
%----------
for i=1:N
    Iline=I0(i,:);
    I_signal_rough=Iline(Iline>Threth_I0);
    Signal_median=median(I_signal_rough);
    Th=Signal_median/2;%threshold for judging a edge
    signal_range=find(Iline>Th);
    if ~isempty(signal_range)
        Edge_cathode(i)=signal_range(1);
        Edge_anode(i)=signal_range(end);
        Iline(1:Edge_cathode(i)-1)=NaN;
        Iline(Edge_anode(i):L)=NaN;
        I0(i,:)=Iline;
    else
        Edge_cathode(i)=1;
        Edge_anode(i)=L;
        Iline=NaN;
        I0(i,:)=Iline;
    end
end

%% re-calculate E-field with different procedure for outliers who satisfis I_over_I0>1
Nom3(:,:,1)=(Output.I_crossed_biased_field(:,:,1)-I0_cross);
A3=Nom3(:,:,1)<0;
Nom3(A3)=0;
I_over_I0=(Nom3(:,:,1)./I0);
I_over_I0_filtered=I_over_I0;
%--deal with outliers-----
outliers=(I_over_I0>1);
[x,y]=find(I_over_I0>1);
for i=1:length(x)
    for j=1:3
        for k=1:3
            X=(x(i)+j-2);
            Y=(y(i)+k-2);
            if ((X==0)||(X==N+1)||(Y==0)||(Y==L+1))
                neighbor(j,k)=NaN;
            else
                neighbor(j,k)=I_over_I0(X,Y);
            end
        end
    end
    neighbor(neighbor>1)=NaN;
    neighbor_mean=mean(mean(neighbor,'omitnan'));
    if isnan(neighbor_mean)
        neighbor_mean=1;
    end
    I_over_I0_filtered(x(i),y(i))=neighbor_mean;
end
%-------------------------
E_clean=(1/Calib.alpha) * (asin(sqrt(I_over_I0_filtered)));
bias_string=num2str(Output.bias);
flux_string=num2str(Output.flux);
figure
imagesc(E_clean)
title(['E-field sensor-' Output.sensor_name '@' bias_string 'V,' flux_string 'mA'])
axes1=gca;
set(axes1,'DataAspectRatio',[1 1 1],'Layer','top');
set(axes1,'View',[-90 90]);
colorbar
box


end
