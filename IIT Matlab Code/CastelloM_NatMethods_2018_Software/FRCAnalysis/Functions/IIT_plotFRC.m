function IIT_plotFRC(FRCresult)
%IIT_FRC_plotFRC Plot the results of the FRC analysis
 
    theta_sess = round(360/(2*pi) * FRCresult.Theta);
    if ~isnan(FRCresult.CutOff.fixedThreshold_smallAngles)
        if isempty(FRCresult.Title)
            T = 'FRC analysis_smallAngles';
        else
            T = ['FRC analysis_smallAngles:',' ', FRCresult.Title];
        end
        figure('Name',T);
        clear title
        plot(FRCresult.Scale, FRCresult.sFRC.smallAngles, '.-b', 'MarkerSize',10);
        hold on
        plot(FRCresult.Scale, FRCresult.FixedThreshold, ':r');
        plot([FRCresult.Scale(FRCresult.CutOff.fixedThreshold_smallAngles) FRCresult.Scale(FRCresult.CutOff.fixedThreshold_smallAngles) ],[ 0 1 ],'--r');
        
        txt = strcat('Resolution = ', num2str(round (FRCresult.Resolution.fixedThreshold_smallAngles)), ' nm. Cutoff frequency = ', ...
            num2str(FRCresult.Scale(FRCresult.CutOff.fixedThreshold_smallAngles)), ' \mum^{-1}' ) ;
        h = text(FRCresult.Scale(FRCresult.CutOff.fixedThreshold_smallAngles)+ 0.5 ,0.95 ,txt,'fontsize',12);
        set(h,'Rotation',-90);
        
        xlabel('Spatial Frequency (1/\mu)','fontsize',14);
        ylabel('FRC','fontsize',14);
        
        if (FRCresult.Theta == 0 || FRCresult.Theta == pi/2)
            t = strcat('FRC analysis, \forall \theta');
        else
            t = strcat('FRC analysis, ', num2str(-theta_sess) ,'\circ < \theta < ', num2str(theta_sess),'\circ');
        end
        title(t ,'fontSize', 15 );
    end
    
    if ~isnan(FRCresult.CutOff.fixedThreshold_largeAngles)
        if isempty(FRCresult.Title)
            T = 'FRC analysis_largeAngles';
        else
            T = ['FRC analysis_largeAngles:',' ', FRCresult.Title];
        end
        
        figure('Name',T);
        clear title
        plot(FRCresult.Scale, FRCresult.sFRC.largeAngles, '.-b', 'MarkerSize',10);
        hold on
        plot(FRCresult.Scale, FRCresult.FixedThreshold, ':r');
        plot([FRCresult.Scale(FRCresult.CutOff.fixedThreshold_largeAngles) FRCresult.Scale(FRCresult.CutOff.fixedThreshold_largeAngles) ],[ 0 1 ],'--r');
        
        txt = strcat('Resolution = ', num2str(round (FRCresult.Resolution.fixedThreshold_largeAngles)), ' nm. Cutoff frequency = ', ...
            num2str(FRCresult.Scale(FRCresult.CutOff.fixedThreshold_largeAngles)), ' \mum^{-1}' ) ;
        h = text(FRCresult.Scale(FRCresult.CutOff.fixedThreshold_largeAngles)+ 0.5 ,0.95 ,txt,'fontsize',12');
        set(h,'Rotation',-90);
        
        xlabel('Spatial Frequency (1/\mu)','fontsize',14);
        ylabel('FRC','fontsize',14);
        
        if (FRCresult.Theta == 0 || FRCresult.Theta == pi/2)
            t = strcat('FRC analysis, \forall \theta');
        else
            t = strcat('FRC analysis, \theta < ', num2str(-theta_sess) ,'\circ U \theta > ', num2str(theta_sess),'\circ' );
        end
        title(t ,'fontSize', 15 );
    end

end

