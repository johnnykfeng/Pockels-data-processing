function [Output]=Calculate_Efield(Output, Calib)

%%Remove background light with crossed polarizers
Nom3(:,:,1)=(Output.I_crossed_biased_field(:,:,1)-Calib.I_crossed_nobias_field);
A3=Nom3(:,:,1)<0;
Nom3(A3)=0;

%%%%%%Correcting for sin (Not using  small angle approximation)
I_over_I0=(Nom3(:,:,1)./Calib.I_parallel_nobias_field);
I_over_I0(I_over_I0>1)=1;
% I_over_I0(I_over_I0>1)=NaN;
E_biased_corrected=(1/Calib.alpha) * (asin(sqrt(I_over_I0)));

Output.I_over_I0=I_over_I0;
Output.E_field_biased_corrected(:,:,1)=E_biased_corrected;
Output.E_cross_section_average_corrected(:,:,1)=mean(E_biased_corrected,'omitnan');

end