clc;
clear all;
clf;

d1 = 1;
d2 = 1;
p = 4;
m = 4;

element_type = 'D2TR3N';

inculusion = 0 ;

R = 0.2;  


if isequal(element_type, 'D2TR3N')

    NPE = 3;
    if inculusion

        NoN = 2*(p+1)*(2*m+1) + 2*(p-1)*(2*m+1) + (p-1)^2;
        NoE = 16*p*m+(p^2)*2;
    
    else
        NoN = 2*(p+1)*(m+1) + 2*(p-1)*(m+1);
        NoE = p*m*8;
    end


elseif isequal(element_type, 'D2TR6N')

    NPE = 6;
    if inculusion

        NoN = 32*p*m + (2*p-1)^2;
        NoE = 16*p*m+2*(p^2);
    
    else
        NoN = (2*(p+1)*(m+1) + 2*(p-1)*(m+1))*2 + 8 * p * m;
        NoE = 8*p*m;
    end


elseif isequal(element_type, 'D2QU8N')

    NPE = 8;
    if inculusion

        NoN = 32*p*m + (2*p-1)^2;
        NoE = 8*p*m+p^2;
    
    else
        NoN = (2*(p+1)*(m+1) + 2*(p-1)*(m+1))*2 + 8 * p * m;
        NoE = 4*p*m;
    end


elseif isequal(element_type, 'D2QU9N')

    NPE = 9;
    if inculusion

        NoN = 32*p*m + (2*p-1)^2;
        NoE = 8*p*m+p^2;
    
    else
        NoN = (2*(p+1)*(m+1) + 2*(p-1)*(m+1))*2 + 8 * p * m;
        NoE = 4*p*m;
    end


elseif isequal(element_type, 'D2QU4N')

    NPE = 4;
    if inculusion

        NoN = 2*(p+1)*(2*m+1) + 2*(p-1)*(2*m+1) + (p-1)^2;
        NoE = 8*p*m+(p^2);
    
    else
        NoN = 2*(p+1)*(m+1) + 2*(p-1)*(m+1);
        NoE = 4*p*m;
    end

end

[NL, EL] = void_mesh(d1, d2, p, m, R, element_type, inculusion);


switch element_type

    case 'D2TR3N'
        
        for i = 1:NoE
            hold on;
            plot([NL(EL(i,1),1), NL(EL(i,2),1)], [NL(EL(i,1),2), NL(EL(i,2),2)],'m');
            plot([NL(EL(i,2),1), NL(EL(i,3),1)], [NL(EL(i,2),2), NL(EL(i,3),2)],'m');
            plot([NL(EL(i,3),1), NL(EL(i,1),1)], [NL(EL(i,3),2), NL(EL(i,1),2)],'m');
            x = (NL(EL(i,1),1) + NL(EL(i,2),1) + NL(EL(i,3),1))/3;
            y = (NL(EL(i,1),2) + NL(EL(i,2),2) + NL(EL(i,3),2))/3;
            text(x,y,num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment','center')
        end
        
        for i = 1:NoN
            hold on;
            plot(NL(i,1),NL(i,2),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1])
            text(NL(i,1),NL(i,2), num2str(i), 'Color','w','FontSize',8,'HorizontalAlignment','center')
        end
        
        axis equal


    case 'D2TR6N'
        
        for i = 1:NoE

            hold on;

            plot([NL(EL(i,1),1), NL(EL(i,4),1)], [NL(EL(i,1),2), NL(EL(i,4),2)],'m');
            plot([NL(EL(i,4),1), NL(EL(i,2),1)], [NL(EL(i,4),2), NL(EL(i,2),2)],'m');
            plot([NL(EL(i,2),1), NL(EL(i,5),1)], [NL(EL(i,2),2), NL(EL(i,5),2)],'m');
            plot([NL(EL(i,5),1), NL(EL(i,3),1)], [NL(EL(i,5),2), NL(EL(i,3),2)],'m');
            plot([NL(EL(i,3),1), NL(EL(i,6),1)], [NL(EL(i,3),2), NL(EL(i,6),2)],'m');
            plot([NL(EL(i,6),1), NL(EL(i,1),1)], [NL(EL(i,6),2), NL(EL(i,1),2)],'m');

            x = (NL(EL(i,6),1) + NL(EL(i,2),1)+ NL(EL(i,1),1) + NL(EL(i,3),1)+ NL(EL(i,4),1) + NL(EL(i,5),1))/6;
            y = (NL(EL(i,6),2) + NL(EL(i,2),2)+ NL(EL(i,1),2) + NL(EL(i,3),2)+ NL(EL(i,4),2) + NL(EL(i,5),2))/6;

            text(x,y,num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment','center')
        end
        NoN = size(NL,1);
        for i = 1:NoN
            hold on;
            plot(NL(i,1),NL(i,2),'o','MarkerSize',8, 'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1])
            text(NL(i,1),NL(i,2), num2str(i), 'Color','w','FontSize',8,'HorizontalAlignment','center')
        end
        
        axis equal


        
    case 'D2QU4N'
        
        for i = 1:NoE
            hold on;
            plot([NL(EL(i,1),1), NL(EL(i,2),1)], [NL(EL(i,1),2), NL(EL(i,2),2)],'m');
            plot([NL(EL(i,2),1), NL(EL(i,3),1)], [NL(EL(i,2),2), NL(EL(i,3),2)],'m');
            plot([NL(EL(i,3),1), NL(EL(i,4),1)], [NL(EL(i,3),2), NL(EL(i,4),2)],'m');
            plot([NL(EL(i,4),1), NL(EL(i,1),1)], [NL(EL(i,4),2), NL(EL(i,1),2)],'m');
            x = (NL(EL(i,1),1) + NL(EL(i,2),1) + NL(EL(i,3),1) + NL(EL(i,4),1))/4;
            y = (NL(EL(i,1),2) + NL(EL(i,2),2) + NL(EL(i,3),2) + NL(EL(i,4),2))/4;
            text(x,y,num2str(i), 'Color', 'k', 'FontSize', 12, 'HorizontalAlignment','center')
        end
        
        
        for i = 1:NoN
            hold on;
            plot(NL(i,1),NL(i,2),'o','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1])
            text(NL(i,1),NL(i,2),num2str(i),'Color','w','FontSize',12,'HorizontalAlignment','center')
        end
        
        axis equal


    case 'D2QU9N'
        NoN = size(NL,1);

        for i = 1:NoN
            hold on;
            
            plot(NL(i,1),NL(i,2),'o','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1])
            
            text(NL(i,1),NL(i,2),num2str(i),'Color','w','FontSize',10,'HorizontalAlignment','center')
        end
        
        
        for i = 1:NoE
            hold on;
            
            plot([NL(EL(i,1),1), NL(EL(i,5),1)], [NL(EL(i,1),2), NL(EL(i,5),2)],'m');
            plot([NL(EL(i,5),1), NL(EL(i,2),1)], [NL(EL(i,5),2), NL(EL(i,2),2)],'m');
            plot([NL(EL(i,2),1), NL(EL(i,6),1)], [NL(EL(i,2),2), NL(EL(i,6),2)],'m');
            plot([NL(EL(i,6),1), NL(EL(i,3),1)], [NL(EL(i,6),2), NL(EL(i,3),2)],'m');
            plot([NL(EL(i,3),1), NL(EL(i,7),1)], [NL(EL(i,3),2), NL(EL(i,7),2)],'m');
            plot([NL(EL(i,7),1), NL(EL(i,4),1)], [NL(EL(i,7),2), NL(EL(i,4),2)],'m');
            plot([NL(EL(i,4),1), NL(EL(i,8),1)], [NL(EL(i,4),2), NL(EL(i,8),2)],'m');
            plot([NL(EL(i,8),1), NL(EL(i,1),1)], [NL(EL(i,8),2), NL(EL(i,1),2)],'m');
            
            x = (NL(EL(i,1),1) + NL(EL(i,2),1) + NL(EL(i,3),1) + NL(EL(i,4),1))/4;
            y = (NL(EL(i,1),2) + NL(EL(i,2),2) + NL(EL(i,3),2) + NL(EL(i,4),2))/4;
            
            text(x,y,num2str(i), 'Color', 'r', 'FontSize', 12, 'HorizontalAlignment','center')
        end
        
        

        axis equal


    case 'D2QU8N'
        
        for i = 1:NoE
            hold on;
            
            plot([NL(EL(i,1),1), NL(EL(i,5),1)], [NL(EL(i,1),2), NL(EL(i,5),2)],'m');
            plot([NL(EL(i,5),1), NL(EL(i,2),1)], [NL(EL(i,5),2), NL(EL(i,2),2)],'m');
            plot([NL(EL(i,2),1), NL(EL(i,6),1)], [NL(EL(i,2),2), NL(EL(i,6),2)],'m');
            plot([NL(EL(i,6),1), NL(EL(i,3),1)], [NL(EL(i,6),2), NL(EL(i,3),2)],'m');
            plot([NL(EL(i,3),1), NL(EL(i,7),1)], [NL(EL(i,3),2), NL(EL(i,7),2)],'m');
            plot([NL(EL(i,7),1), NL(EL(i,4),1)], [NL(EL(i,7),2), NL(EL(i,4),2)],'m');
            plot([NL(EL(i,4),1), NL(EL(i,8),1)], [NL(EL(i,4),2), NL(EL(i,8),2)],'m');
            plot([NL(EL(i,8),1), NL(EL(i,1),1)], [NL(EL(i,8),2), NL(EL(i,1),2)],'m');
            
            x = (NL(EL(i,1),1) + NL(EL(i,2),1) + NL(EL(i,3),1) + NL(EL(i,4),1))/4;
            y = (NL(EL(i,1),2) + NL(EL(i,2),2) + NL(EL(i,3),2) + NL(EL(i,4),2))/4;
            
            text(x,y,num2str(i), 'Color', 'r', 'FontSize', 12, 'HorizontalAlignment','center')
        end
        
        NoN = size(NL,1);
        
        for i = 1:NoN
            hold on;
            if ismember(i,EL)
            
                plot(NL(i,1),NL(i,2),'o','MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1])
                
                text(NL(i,1),NL(i,2),num2str(i),'Color','w','FontSize',12,'HorizontalAlignment','center')
            end
        end
        
        axis equal
        
end