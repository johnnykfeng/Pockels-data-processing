function [Output]=Calculate_Efield(Output, alpha)

%%Remove background light with crossed polarizers
Nom3(:,:,1)=(Output.I_crossed_biased_field(:,:,1)-Output.I_crossed_nobias_field);
A3=Nom3(:,:,1)<0;
Nom3(A3)=0;

%       %%small angle approximation
%         E_biased(:,:,1)=(1/alpha) * (sqrt(Nom3(:,:,1)./I_parallel_nobias_field));
%
%
%         %%% Cross section of field calculation
%         E_cross_section_average(:,:,1)=mean(E_biased,'omitnan');
%
%
%          figure
%
%         plot([1:520],E_cross_section_average(:,:,1));
%         title(['Pockels average profile - sensor=' detector_num 'bias=' bias, ', Xray tube current=' flux ' , ' time])
%

%%%%%%Correcting for sin (Not using  small angle approximation)
I_over_I0=(Nom3(:,:,1)./Output.I_parallel_nobias_field);
I_over_I0(I_over_I0>1)=1;
E_biased_corrected=(1/alpha) * (asin(sqrt(I_over_I0)));

Output.E_field_biased_corrected(:,:,1)=E_biased_corrected;
Output.E_cross_section_average_corrected(:,:,1)=mean(E_biased_corrected,'omitnan');

end