function [Output]=Read_Pockels_data


prompt = {'Enter sensor name','Enter bias:','Xray tube current (enter 0 if no Xray)'};
dlgtitle = 'Test information ';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
sensor_name=char(answer{1,1});

% 
Output.sensor_name=sensor_name;
Output.bias=answer{2,1};
Output.flux=answer{3,1};

waitfor(msgbox(['Select the 0V, parallel data']));
[file_0V_parallel,path_0V_parallel] = uigetfile('*.txt');
Output.data_directory=path_0V_parallel;

buffer1 = fileread(char([path_0V_parallel  file_0V_parallel]));
data1 = textscan(buffer1,'%f %f %f', 'headerlines', 19) ;
I_parallel_nobias_matrix(:,1)=data1{1,1};
I_parallel_nobias_matrix(:,2)=data1{1,2};
I_parallel_nobias_matrix(:,3)=data1{1,3};
%[I_parallel_nobias_matrix(:,1) I_parallel_nobias_matrix(:,2) I_parallel_nobias_matrix(:,3)] = textread([path_0V_parallel  file_0V_parallel], '%f %f %f', 'headerlines', 19) ;

waitfor(msgbox(['Select the 0V, crossed data']));
[file_0V_crossed,path_0V_crossed] = uigetfile('*.txt');

buffer2 = fileread(char([path_0V_crossed  file_0V_crossed]));
data2 = textscan(buffer2,'%f %f %f', 'headerlines', 19) ;
I_crossed_nobias_matrix(:,1)=data2{1,1};
I_crossed_nobias_matrix(:,2)=data2{1,2};
I_crossed_nobias_matrix(:,3)=data2{1,3};
%[I_crossed_nobias_matrix(:,1) I_crossed_nobias_matrix(:,2) I_crossed_nobias_matrix(:,3)] = textread([path_0V_crossed  file_0V_crossed], '%f %f %f', 'headerlines', 19) ;


waitfor(msgbox(['Select the crossed data under this bias and flux']));
[file_biased_crossed,path_biased_crossed] = uigetfile('*.txt');
%[file_biased_crossed,path_biased_crossed] = uigetfile('*.txt');

buffer3 = fileread([path_biased_crossed  file_biased_crossed]);
data3 = textscan(buffer3,'%f %f %f', 'headerlines', 19) ;
I_crossed_biased_matrix(:,1)=data3{1,1};
I_crossed_biased_matrix(:,2)=data3{1,2};
I_crossed_biased_matrix(:,3)=data3{1,3};
%[I_crossed_biased_matrix(:,1,1) I_crossed_biased_matrix(:,2,1) I_crossed_biased_matrix(:,3,1)] = textread([path_biased_crossed  file_biased_crossed], '%f %f %f', 'headerlines', 19) ;

% choice = menu('Is there more biased input you would like to analyze?','Yes','No');
% if choice==2 | choice==0
%    break;
% end
% n=n+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  
  %%% Reshapes the 3rd column into a 2D matrix
        I_parallel_nobias_field=reshape(I_parallel_nobias_matrix(:,3),[696,520]);
        I_crossed_nobias_field=reshape(I_crossed_nobias_matrix(:,3),[696,520]);
        I_crossed_biased_field(:,:,1)=reshape(I_crossed_biased_matrix(:,3,1),[696,520]);
   
   
        Output.I_parallel_nobias_field=(single(I_parallel_nobias_field));
        Output.I_crossed_nobias_field=(single(I_crossed_nobias_field));
        Output.I_crossed_biased_field(:,:,1)=(single(I_crossed_biased_field(:,:,1)));
     
end