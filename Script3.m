%This lets you plot your 4D data to compare parameter distributions
%load('chi2_values_ADH4.mat')
%load('chi2_values_SH9.mat')
%or load whatever data you want to compare
for i = 1:length(elongation_input)
    figure(i)
    Red = chi2_values(:,:,:,i) > 0.7; %This thresholds your parameters, variable stringency
    Yellow = 0.2<chi2_values(:,:,:,i) & chi2_values(:,:,:,i) <0.7;
    Any = chi2_values(:,:,:,i); %This will find any match even if very small. 
    R = find(Any);
    [S,T,V] = ind2sub(size(Any),R);
    C = find(Red);
    N = find(Yellow);
    [D,E,F] = ind2sub(size(Red),C);
    [J,K,L] = ind2sub(size(Yellow),N);
    figure(i)
    scatter3(D,E,F,25,'bs', 'LineWidth', 4)
    %hold on
    scatter3(J,K,L,25,'bs', 'LineWidth', 4)
    hold on
    %scatter3(S,T,V,25,'bs' ,'LineWidth', 4) %Don't normally plot as is
    %very non-specific and picks up on any time similarity
    hold on

    Red2 = chi2_values_cytoplasm(:,:,:,i) > 0.7;
    Yellow2 = 0.02<chi2_values_cytoplasm(:,:,:,i) & chi2_values_cytoplasm(:,:,:,i) <0.7;
    Any2 = chi2_values_cytoplasm(:,:,:,i);
    R2 = find(Any2);
    [S2,T2,V2] = ind2sub(size(Any2),R2);
    C2 = find(Red2);
    N2 = find(Yellow2);
    [D2,E2,F2] = ind2sub(size(Red2),C2);
    [J2,K2,L2] = ind2sub(size(Yellow2),N2);
    figure(i) %(i+length(elongation_input)) 
    scatter3(D2,E2,F2,18,'co', 'MarkerFaceColor', 'c')
    hold on
    scatter3(J2,K2,L2,18,'co', 'MarkerFaceColor', 'c')
  %  hold on
    %scatter3(S2,T2,V2,18,'go', 'MarkerFaceColor', 'g')
    
    title(sprintf('Elongation Rate = %g kb per min', elongation_input(i)), 'FontName', 'Cambria Math', 'FontSize', 18 )
    xlabel('Av. time off (min)', 'FontName', 'Cambria Math', 'FontSize', 18)
    set(gca,'XTick', [1:2:length(av_time_off_input)])
    set(gca,'XTickLabel', av_time_off_input(1:2:end), 'FontName', 'Cambria Math', 'FontSize', 18)

    ylabel('Av. time on (min)', 'FontName', 'Cambria Math', 'FontSize', 18)
    set(gca,'YTick', [1:2:length(av_time_on_input)])
    set(gca,'YTickLabel', av_time_on_input(1:2:end), 'FontName', 'Cambria Math', 'FontSize', 18)
    zlabel('c (transcripts per min)', 'FontName', 'Cambria Math', 'FontSize', 18)
    %legend('GAL1:ADH1T nuclear transcripts', 'TATA-box mutant nuclear transcripts')
   set(gca,'ZTick', [1:2:length(c_input)])
    set(gca,'ZTickLabel', c_input(1:2:end), 'FontName', 'Cambria Math', 'FontSize', 18)
end
