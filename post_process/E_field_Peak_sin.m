% to judge whether the calculated E-field reaches the peak of sin()
% function [not done]

%% define signal range
version=str2double(Output.version);
if version<2
    I0=Output.I_parallel_nobias_field;
else
    I0=Calib.I_parallel_nobias_field;
end

dimention=size(Output.I_crossed_biased_field);
N=dimention(1);
L=dimention(2);

for i=1:N
    Iline=I0(i,:);
    I_signal_rough=Iline(Iline>200);
    Signal_median=median(I_signal_rough);
    Th=Signal_median/2;%threshold for judging a edge
    signal_range=find(Iline>Th);
    if ~isempty(signal_range)
        Edge_left(i)=signal_range(1);
        Edge_right(i)=signal_range(end);
        Iline(1:Edge_left(i)-1)=NaN;
        Iline(Edge_right(i):L)=NaN;
        I0(i,:)=Iline;
    else
        Edge_left(i)=1;
        Edge_right(i)=L;
        Iline=NaN;
        I0(i,:)=Iline;
    end
end

%% re-calculate E-field with different procedure for outliers who satisfis I_over_I0>1
Nom3(:,:,1)=(Output.I_crossed_biased_field(:,:,1)-Calib.I_crossed_nobias_field);
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
    
% 
%     neighbor(1,1)=I_over_I0((x(i)-1),(y(i)-1));
%     
%     neighbor(1,2)=I_over_I0((x(i)),(y(i)-1));
%     
%     neighbor(1,3)=I_over_I0((x(i)+1),(y(i)-1));
%     
%         
%     neighbor(2,1)=I_over_I0((x(i)-1),(y(i)));
%     
%     neighbor(2,2)=I_over_I0((x(i)),(y(i)));
%     
%     neighbor(2,3)=I_over_I0((x(i)+1),(y(i)));
%     
%         
%     neighbor(3,1)=I_over_I0((x(i)-1),(y(i)+1));
%     
%     neighbor(3,2)=I_over_I0((x(i)),(y(i)+1));
%     
%     neighbor(3,3)=I_over_I0((x(i)+1),(y(i)+1));
end
%-------------------------
E_biased_corrected=(1/Calib.alpha) * (asin(sqrt(I_over_I0_filtered)));

imagesc(E_biased_corrected)
% title(['E-field sensor-' Output.sensor_name '@' bias_string 'V,' flux_string 'mA'])
axes1=gca;
set(axes1,'DataAspectRatio',[1 1 1],'Layer','top');
set(axes1,'View',[-90 90]);
colorbar
box



E_max=(1/Calib.alpha) * (asin(1));