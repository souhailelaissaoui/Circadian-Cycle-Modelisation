function bio341_mini_projet1()
% 14/11/2017

clear all
close all

what_to_solve = [1 0] ;
%what_to_solve(1) : Partie A
%what_to_solve(2)  : Partie B


%% PARTIE A
if(what_to_solve(1))
    
    %% SIMULATION OF THE MODEL
    
    %% Question 4
    %trajectoires  - Euler
    % dt1=0.01;
    % n1=10000/dt1;
    % x=zeros(1,n1);
    % y=zeros(1,n1);
    % z=zeros(1,n1);
    % x(1)=0; %x0
    % y(1)=1; %y0
    % z(1)=5; %z0
    %
    % for i=1:n1-1
    %    x(i+1)=x(i)+dt1*(v1 * K1^4/(K1^4+z(i)^4)- v2 * x(i)/(K2+x(i))) ;
    %    y(i+1)=y(i)+dt1*(k3 * x(i)- v4 * y(i)/(K4+y(i)));
    %    z(i+1)=z(i)+dt1*(k5 * y(i)- v6 * z(i)/(K6+z(i)));
    % end
    %
    % figure
    % subplot(3,1,1)
    % plot(x) %idealement, il faut faire plot(t,x)
    %
    % subplot(3,1,2)
    % plot(y) %plutot plot(t,y)
    %
    % subplot(3,1,3)
    % plot(z) %plutot plot(t,z)
    %
    % figure
    % plot3(x,y,z); hold on
    
    % ODE45
    dt = 0.1 ; %timestep
    tspan = 0 : dt : 1000; %10 000 heures = 416.67 jours
    condition_ini = [0 ; 1 ; 5] ;
    
    v1=0.7 ;
    [t, sol] = ode45(@(t,x) equation(t,x,v1), tspan, condition_ini);
    x=sol(:,1);
    y=sol(:,2);
    z=sol(:,3);
   
    figure
    subplot(3,1,1)
    plot(t/24,x)
    ylabel('x (nM)')
    
    subplot(3,1,2)
    plot(t/24,y)
    ylabel('y (nM)')
    
    subplot(3,1,3)
    plot(t/24,z)
    ylabel('z (nM)')
    xlabel('time [day]')
    
    figure
    plot3(x,y,z)
    xlabel('x') ; ylabel('y') ;zlabel('z')
    
    
    
    %% Question 5
    
    figure ;
    subplot(3,1,1)
    facteur_fin = 0.8 ; % on prend les valeurs de 80% à 100% de la fin
    x_end = x(length(x)*facteur_fin:end) ;
    [ppx, f] = periodogram(x_end-mean(x),[],[],1/dt) ;
    plot(1./f, ppx)
    ylabel('ppx')
    xlim([0 100])
    text(22.76, 0.16,'\leftarrow T=22.76');
    subplot(3,1,2)
    y_end = y(length(y)*facteur_fin:end) ;
    [ppy, f] = periodogram(y_end-mean(y),[],[],1/dt) ;
    plot(1./f, ppy)
    ylabel('ppy')
    xlim([0 100])
    text(22.76, 0.74,'\leftarrow T=22.76');
    subplot(3,1,3)
    z_end = z(length(z)*facteur_fin:end) ;
    [ppz, f] = periodogram(z_end-mean(z),[],[],1/dt) ;
    plot(1./f, ppz)
    ylabel('ppz')
    xlabel('T period [h]')
    xlim([0 100])
    text(22.76, 4.676,'\leftarrow T=22.76');
    
    %% Question 6
    
    v1=linspace(0.001, 5, 5) ; % !!!! changer le 10 en un nombre plus grand
    x_min=zeros(1,length(v1));
    x_max=zeros(1,length(v1));
    period=zeros(1,length(v1));
    pf_x=zeros(1,length(v1));
    pf_y=zeros(1,length(v1));
    pf_z=zeros(1,length(v1));
    
    figure
    
    for j=1:length(v1)
        
        [~, sol] = ode45(@(t,x) equation(t,x,v1(j)), tspan, condition_ini);
        x=sol(:,1);
        y=sol(:,2);
        z=sol(:,3);
        
        plot3(x,y,z) ; hold on
        xlabel('x') ; ylabel('y') ;zlabel('z')
        
        %     x=zeros(1,n1);
        %     y=zeros(1,n1);
        %     z=zeros(1,n1);
        %     x(1)=0; %x0
        %     y(1)=1; %y0
        %     z(1)=5; %z0
        %
        %     for i=1:n1-1
        %        x(i+1)=x(i)+dt1*(v1(j) * K1^4/(K1^4+z(i)^4)- v2 * x(i)/(K2+x(i))) ;
        %        y(i+1)=y(i)+dt1*(k3 * x(i)- v4 * y(i)/(K4+y(i)));
        %        z(i+1)=z(i)+dt1*(k5 * y(i)- v6 * z(i)/(K6+z(i)));
        %     end
        %
        %     plot3(x,y,z); hold on
        
        % question 7
        x_min(j)=min(x(length(x)*facteur_fin : end)) ;
        x_max(j)=max(x(length(x)*facteur_fin : end)) ;
        
        %question 8 - 9
        if v1(j)>=0.3 && v1(j)<=4
            %question 8
            x_end = x(length(x)*facteur_fin:end) ;
            [ppx, f] = periodogram(x_end-mean(x),[],[],1/dt) ;
            [~,index] = max(ppx) ;
            period(j) = 1./f(index) ;
            
            %question 9
            [sol_pf] = fsolve(@(x) solve_points_fixes(x,v1(j)),[1;1;2]);
            pf_x(j) = sol_pf(1) ;
            pf_y(j) = sol_pf(2) ;
            pf_z(j) = sol_pf(3);
        end
    end
    
    %% Question 7 (voir question 6)
    figure
    subplot(1,2,1)
    plot(v1,x_min,'-o', v1, x_max,'-x') ; hold on
    xlabel('v1')
    
    %% Question 8 (voir question 6)
    
    % v1=linspace(0.3, 1.3, 20) ;
    % period=zeros(1,length(v1));
    %
    % for j=1:length(v1)
    %     x=zeros(1,n1);
    %     y=zeros(1,n1);
    %     z=zeros(1,n1);
    %     x(1)=0; %x0
    %     y(1)=1; %y0
    %     z(1)=5; %z0
    %
    %     for i=1:n1-1
    %        x(i+1)=x(i)+dt1*(v1(j) * K1^4/(K1^4+z(i)^4)- v2 * x(i)/(K2+x(i))) ;
    %        y(i+1)=y(i)+dt1*(k3 * x(i)- v4 * y(i)/(K4+y(i)));
    %        z(i+1)=z(i)+dt1*(k5 * y(i)- v6 * z(i)/(K6+z(i)));
    %     end
    %     [ppx, f] = periodogram(x(50000:length(x))-mean(x),[],[],1/dt1) ;
    %     [maxi,index] = max(ppx) ;
    %     period(j) = 1./f(index) ;
    % end
    subplot(1,2,2)
    plot(v1, period, '-o')
    xlabel('v1')
    ylabel('period [h]')
    xlim([0.3 4])
    
    %% BEHAVIOR NEAR THE HOPF BIFURCATION
    
    %% Question 9
    subplot(1,2,1)
    plot(v1, pf_x,'*') ;
    legend('x_{min}', 'x_{max}','x_{PF}')
    
end
if(what_to_solve(2))
    dt = 0.1 ; %timestep
    tspan = 0 : dt : 10000; %10 000 heures = 416.67 jours
    condition_ini = [0 ; 1 ; 5] ;
    %% Question 10
    figure
    
    v1 = 0.5 ;
    
    [~, sol] = ode45(@(t,x) equation(t,x,v1), tspan, condition_ini);
    x=sol(:,1);
    y=sol(:,2);
    z=sol(:,3);

    plot3(x,y,z) ; hold on
    xlabel('x') ; ylabel('y') ;zlabel('z');
    
    [sol_pf] = fsolve(@(x) solve_points_fixes(x,v1),[1;1;2]);
    pf10_x = sol_pf(1) ;
    pf10_y = sol_pf(2) ;
    pf10_z = sol_pf(3);
    
    figure
    plot3(x,y,z) ; hold on
    xlim([pf10_x-0.1, pf10_x+0.1])
    ylim([pf10_y-0.1, pf10_y+0.1])
    zlim([pf10_z-0.1, pf10_z+0.1])
    xlabel('x') ; ylabel('y') ;zlabel('z')
    plot3(pf10_x, pf10_y, pf10_z, '-o', 'LineWidth',2) ; hold on
    
    figure
    plot(x,y)
    Jaco = jacobien(pf10_x,  pf10_y, pf10_z, v1)
    [V, D] = eig(Jaco) 
    
end

    function dxdt = equation(~,x,v1)
        %equations of the system
        %x(1) = x ; x(2) = y ; x(3) = z
        
        % Parameters
        %v1 est en parametre de la fonction car peut varier en fonciton de la
        %question
        v2 = 0.35 ;
        v4 = v2 ;
        v6 = v2 ;
        
        K1 = 1 ;
        K2 = K1 ;
        K4 = K1 ;
        K6 = K1 ;
        
        k3 = 0.7 ;
        k5 = k3 ;
        
        dx = v1 * K1^4/(K1^4+x(3)^4)- v2 * x(1)/(K2+x(1));
        dy = k3 * x(1)- v4 * x(2)/(K4+x(2)) ;
        dz = k5 * x(2)- v6 * x(3)/(K6+x(3)) ;
        
        dxdt=[dx; dy; dz];
    end

    function F = solve_points_fixes(x,v1)
        %equations of the system
        %x(1) = x ; x(2) = y ; x(3) = z
        
        % Parameters
        v2 = 0.35 ;
        v4 = v2 ;
        v6 = v2 ;
        
        K1 = 1 ;
        K2 = K1 ;
        K4 = K1 ;
        K6 = K1 ;
        
        k3 = 0.7 ;
        k5 = k3 ;
        
        F = [ v1 * K1^4/(K1^4+x(3)^4)- v2 * x(1)/(K2+x(1));
            k3 * x(1)- v4 * x(2)/(K4+x(2)) ;
            k5 * x(2)- v6 * x(3)/(K6+x(3))] ;
    end

    function J = jacobien(x,y,z, v1)
        v2 = 0.35 ;
        v4 = v2 ;
        v6 = v2 ;
        
        K1 = 1 ;
        K2 = K1 ;
        K4 = K1 ;
        K6 = K1 ;
        
        k3 = 0.7 ;
        k5 = k3 ;
        J = [-v2*K2/(K2+x)^2 ,   0,       -v1*K1^4*4*z^3/(K1^4+ z^4)^2 ;
            k3,             -v4*K4/(K4+y)^2,        0;
            0,                 k5,             -v6*K6/(K6+z)^2];
    end

end
