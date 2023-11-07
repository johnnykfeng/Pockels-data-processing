function [Output]=Read_Pockels_data_group_extraction(sensor_name,bias,X_ray_tube_current)
%
Output.sensor_name=sensor_name;
Output.bias=bias;
Output.flux=X_ray_tube_current;
bias_string=num2str(bias);
flux_string=num2str(X_ray_tube_current);
waitfor(msgbox(['Select the crossed data under ' bias_string 'V, ' flux_string 'mA']));
[file_biased_crossed,path_biased_crossed] = uigetfile('*.txt');


buffer3 = fileread([path_biased_crossed  file_biased_crossed]);
data3 = textscan(buffer3,'%f %f %f', 'headerlines', 19) ;
I_crossed_biased_matrix(:,1)=data3{1,1};
I_crossed_biased_matrix(:,2)=data3{1,2};
I_crossed_biased_matrix(:,3)=data3{1,3};




% %%% Reshapes the 3rd column into a 2D matrix
% I_parallel_nobias_field=reshape(I_parallel_nobias_matrix(:,3),[696,520]);
% I_crossed_nobias_field=reshape(I_crossed_nobias_matrix(:,3),[696,520]);
I_crossed_biased_field(:,:,1)=reshape(I_crossed_biased_matrix(:,3,1),[696,520]);


% Output.I_parallel_nobias_field=(single(I_parallel_nobias_field));
% Output.I_crossed_nobias_field=(single(I_crossed_nobias_field));
Output.I_crossed_biased_field(:,:,1)=(single(I_crossed_biased_field(:,:,1)));


end