function [NL, EL] = void_square_mesh(d1, d2, p, m, c, element_type)

   inc = 0;

    q = [0 0; d1 0; 0 d2; d1 d2];
    
    PD = 2;
    
    NoN = 2*(p+1)*(m+1) + 2*(p-1)*(m+1);
    
    NoE = 4*p*m;
    
    NPE = 4;
    
    NL = zeros(NoN, PD);
    
    a = d1/p;
    b = d2/p;
    %region 1
    
    coor11 = zeros((p+1)*(m+1), PD);
    
    for i = 1:p+1
        
        coor11(i,1) = q(1,1) + (i-1)*a;
        coor11(i,2) = q(1,2);
        
    end
    
    for i = 1:p+1
        
        coor11(m*(p+1)+i,1) = ((d1-c)/2)+(i-1)*(c/p);
        coor11(m*(p+1)+i,2) = ((d2-c)/2);
        
    end
    for i = 1:m-1
        
        for j = 1:p+1
            
            dx = (coor11(m*(p+1)+j,1) - coor11(j,1))/m;
            dy = (coor11(m*(p+1)+j,2) - coor11(j,2))/m;
            
            coor11(i*(p+1)+j,1) = coor11(((i-1)*(p+1))+j,1) + dx;
            coor11(i*(p+1)+j,2) = coor11(((i-1)*(p+1))+j,2) + dy;
            
        end
        
    end
    
    % for i = 1:(p+1)*(m+1)
    %     hold on;
    %     plot(coor11(i,1),coor11(i,2),'o','MarkerSize',30, 'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1])
    %     text(coor11(i,1),coor11(i,2), num2str(i), 'Color','w','FontSize',20,'HorizontalAlignment','center')
    % end
    
    
    
    
    %region 2
    
    coor22 = zeros((p+1)*(m+1), PD);
    
    for i = 1:p+1
        
        coor22(i,1) = q(3,1) + (i-1)*a;
        coor22(i,2) = q(3,2);
        
    end
    
    for i = 1:p+1
        
        coor22(m*(p+1)+i,1) = ((d1-c)/2)+(i-1)*(c/p);
        coor22(m*(p+1)+i,2) = ((d2+c)/2);
        
    end
    
    for i = 1:m-1
        
        for j = 1:p+1
            
            dx = (coor22(m*(p+1)+j,1) - coor22(j,1))/m;
            dy = (coor22(m*(p+1)+j,2) - coor22(j,2))/m;
            
            coor22(i*(p+1)+j,1) = coor22(((i-1)*(p+1))+j,1) + dx;
            coor22(i*(p+1)+j,2) = coor22(((i-1)*(p+1))+j,2) + dy;
            
        end
        
    end
    % for i = 1:(p+1)*(m+1)
    %     hold on;
    %     plot(coor22(i,1),coor22(i,2),'o','MarkerSize',30, 'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1])
    %     text(coor22(i,1),coor22(i,2), num2str(i), 'Color','w','FontSize',20,'HorizontalAlignment','center')
    % end
    
    
    %Region 3
    
    coor33 = zeros((p-1)*(m+1), PD);
    
    for i = 1:p-1
        
        coor33(i,1) = q(1,1);
        coor33(i,2) = q(1,2) + (i)*b;
        
    end
    
    for i = 1:p-1
        
        coor33(m*(p-1)+i,1) = (d1-c)/2;
        coor33(m*(p-1)+i,2) = ((d2-c)/2)+(i)*(c/p);
        
    end
    
    for i = 1:m-1
        
        for j = 1:p-1
            
            dx = (coor33(m*(p-1)+j,1) - coor33(j,1))/m;
            dy = (coor33(m*(p-1)+j,2) - coor33(j,2))/m;
            
            coor33(i*(p-1)+j,1) = coor33(((i-1)*(p-1))+j,1) + dx;
            coor33(i*(p-1)+j,2) = coor33(((i-1)*(p-1))+j,2) + dy;
            
        end
    
        
    end
    % for i = 1:(p-1)*(m+1)
    %     hold on;
    %     plot(coor33(i,1),coor33(i,2),'o','MarkerSize',30, 'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1])
    %     text(coor33(i,1),coor33(i,2), num2str(i), 'Color','w','FontSize',20,'HorizontalAlignment','center')
    % end
    
    
    
    %region4
    
    
    
    coor44 = zeros((p-1)*(m+1), PD);
    
    for i = 1:p-1
        
        coor44(i,1) = q(2,1);
        coor44(i,2) = q(2,2) + (i)*b;
        
    end
    
    for i = 1:p-1
        
        coor44(m*(p-1)+i,1) = (d1+c)/2;
        coor44(m*(p-1)+i,2) = ((d2-c)/2)+(i)*(c/p);
        
    end
    
    for i = 1:m-1
        
        for j = 1:p-1
            
            dx = (coor44(m*(p-1)+j,1) - coor44(j,1))/m;
            dy = (coor44(m*(p-1)+j,2) - coor44(j,2))/m;
            
            coor44(i*(p-1)+j,1) = coor44(((i-1)*(p-1))+j,1) + dx;
            coor44(i*(p-1)+j,2) = coor44(((i-1)*(p-1))+j,2) + dy;
            
        end
        
    end
    % for i = 1:(p+1)*(m+1)
    %     hold on;
    %     plot(coor11(i,1),coor11(i,2),'o','MarkerSize',30, 'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1])
    %     text(coor11(i,1),coor11(i,2), num2str(i), 'Color','w','FontSize',20,'HorizontalAlignment','center')
    % end
    % axis equal
    % for i = 1:(p+1)*(m+1)
    %     hold on;
    %     plot(coor22(i,1),coor22(i,2),'o','MarkerSize',30, 'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1])
    %     text(coor22(i,1),coor22(i,2), num2str(i), 'Color','w','FontSize',20,'HorizontalAlignment','center')
    % end
    % axis equal
    % for i = 1:(p-1)*(m+1)
    %     hold on;
    %     plot(coor33(i,1),coor33(i,2),'o','MarkerSize',30, 'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1])
    %     text(coor33(i,1),coor33(i,2), num2str(i), 'Color','w','FontSize',20,'HorizontalAlignment','center')
    % end
    % axis equal
    % for i = 1:(p-1)*(m+1)
    %     hold on;
    %     plot(coor44(i,1),coor44(i,2),'o','MarkerSize',30, 'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1])
    %     text(coor44(i,1),coor44(i,2), num2str(i), 'Color','w','FontSize',20,'HorizontalAlignment','center')
    % end
    % 
    % axis equal
    for i = 1:m+1
        
        NL((i-1)*4*p+1:(i)*4*p,:) = [coor11((i-1)*(p+1)+1:(i)*(p+1),:);
                                     coor44((i-1)*(p-1)+1:(i)*(p-1),:);
                                     flipud(coor22((i-1)*(p+1)+1:(i)*(p+1),:));
                                     flipud(coor33((i-1)*(p-1)+1:(i)*(p-1),:))];
        
    end
%     for i = 1:NoN
%         hold on;
%         plot(NL(i,1),NL(i,2),'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1])
%         text(NL(i,1),NL(i,2),num2str(i),'Color','w','FontSize',20,'HorizontalAlignment','center')
%     end
%     
%     axis equal
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Elements   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


EL = zeros(NoE, NPE);

if isequal(element_type, 'D2QU4N')

    for i = 1:m*(1+inc) 
        
        for j = 1:4*p
            
            if j == 1
                
                EL((i-1)*(4*p)+j,1) = (i-1)*(4*p) + 1;
                EL((i-1)*(4*p)+j,2) = EL((i-1)*(4*p)+j,1) + 1;
                EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j,1) + 4*p;
                EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j,4) + 1;
                
                
            elseif j == 4*p
                
                EL((i-1)*(4*p)+j,1) = (i)*(4*p);
                EL((i-1)*(4*p)+j,2) = (i-1)*(4*p) + 1;
                EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j,1) + 1;
                EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j,1) + 4*p;
                
            else
                
                EL((i-1)*(4*p)+j,1) = EL((i-1)*(4*p)+j-1,2);
                EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j-1,3);
                EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j,4) + 1;
                EL((i-1)*(4*p)+j,2) = EL((i-1)*(4*p)+j,1) + 1;
                
            end
        end
    end
    
    
    
    
    
    if inc
        
        n=8*p*m; % Number of elements up to now
        NN = 2*(p+1)*(2*m+1) + 2*(p-1)*(2*m+1) - p*4 + 1; % Node number of the down-left corner
        
        for i = 1:p
            for j = 1:p
                
                if (i == 1)&&(j == 1)
                    EL((i-1)*p+j+n,1) =                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          NN;
                    EL((i-1)*p+j+n,2) = NN+1;
                    EL((i-1)*p+j+n,4) = NN+4*p-1;
                    EL((i-1)*p+j+n,3) = NN+4*p;

                elseif (i == 1)&&(j == p)
                    EL((i-1)*p+j+n,1) = NN+p-1;
                    EL((i-1)*p+j+n,2) = NN+p;
                    EL((i-1)*p+j+n,4) = NN+4*p+(p-2);
                    EL((i-1)*p+j+n,3) = NN+p+1;

                elseif (j == 1)&&(i == p)
                    EL((i-1)*p+j+n,1) = NN+3*p+1;
                    EL((i-1)*p+j+n,2) = NoN-(p-2);
                    EL((i-1)*p+j+n,4) = NN+3*p;
                    EL((i-1)*p+j+n,3) = NN+3*p-1;

                elseif (j == p)&&(i == p) 
                    EL((i-1)*p+j+n,1) = NoN;
                    EL((i-1)*p+j+n,2) = NN+2*p-1;
                    EL((i-1)*p+j+n,4) = NN+2*p+1;
                    EL((i-1)*p+j+n,3) = NN+2*p;
                
                elseif i==1
                    EL((i-1)*p+j+n,1) = NN+j-1;
                    EL((i-1)*p+j+n,2) = NN+j;
                    EL((i-1)*p+j+n,4) = NN+4*p-2+j;
                    EL((i-1)*p+j+n,3) = NN+4*p-1+j;

                elseif i==p
                    EL((i-1)*p+j+n,1) = NoN-(p-2)-2+j;
                    EL((i-1)*p+j+n,2) = NoN-(p-2)-1+j;
                    EL((i-1)*p+j+n,4) = NN+3*p+1-j;
                    EL((i-1)*p+j+n,3) = NN+3*p-j;

                elseif j == 1
                    EL((i-1)*p+j+n,1) = NN+4*p+1-i;
                    EL((i-1)*p+j+n,2) = NN+4*p+(p-1)*(i-2);
                    EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n,1) - 1;
                    EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2) + (p-1);

                elseif j == p
                    EL((i-1)*p+j+n,1) = NN+4*p+(p-2)+(p-1)*(i-2);
                    EL((i-1)*p+j+n,2) = NN+p-1+i;
                    EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n,1) + (p-1);
                    EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2) + 1;

                else
                    
                    EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,2);
                    EL((i-1)*p+j+n,4) = EL((i-1)*p+j-1+n,3);
                    EL((i-1)*p+j+n,2) = EL((i-1)*p+j+n,1) + 1;
                    EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,4) + 1;

                end
                
            end
        end
        
    end
        
    end



if isequal(element_type, 'D2QU9N')

    if inc
        NL_Temp = NL(1:4*p*(2*m+1),:);
        [NL_Temp,NoN] = NL_Extender(NL_Temp,p,m,inc);
        NL_Temp (8*p*(4*m+1):(8*p*(4*m+1))+(p-1)^2,:) = NL(4*p*(2*m+1):8*p*m+(p+1)^2,:);
        NL = NL_Temp;
        NL(8*p*(4*m+1),:) = (NL(8*p*(4*m+1)-1,:) + NL(8*p*(4*m+1)-(8*p)+1,:))/2;

        for i = 1:((p-1)*4) 
            if (i < p) 
                NL(8*p*(4*m+1)+(p-1)^2+i,:) = (NL(8*p*(4*m+1)-(8*p)+2*i+1,:) + NL(8*p*(4*m+1)+i,:))/2;
            elseif (i < 3*(p-1)+1) && (i>p-1)
                NL(8*p*(4*m+1)+(p-1)^2+i,:) = (NL(8*p*(4*m+1)+i-p+1,:) + NL(8*p*(4*m+1)+i,:))/2;
            else
                NL(8*p*(4*m+1)+(p-1)^2+i,:) = (NL(8*p*(4*m+1)-(p-1)+i,:) + NL(8*p*(4*m+1)-((p-1)*2)-1-(2*(i-(p-1)^2)),:))/2;
            end
        end
        for i = 1:p 
            for j = 1:(2*p-1)
                NL(8*p*(4*m+1)+(p-1)^2+((2*p-1)*(i-1))+j+((p-1)*4),1) = ((8-((2*i)-1))*(NL(8*p*(4*m+1)-j+1,1)) + (((2*i)-1)*NL(8*p*(4*m+1)-8*p+1+2*p+j,1)))/8;
                NL(8*p*(4*m+1)+(p-1)^2+((2*p-1)*(i-1))+j+((p-1)*4),2) = (NL(8*p*(4*m+1)-j+1,2) + NL(8*p*(4*m+1)-8*p+1+2*p+j,2))/2;
            end
           
        end

        NL_Temp = NL((4*m+1)*8*p+1:size(NL,1),:);
        NL_Temp = sortrows(NL_Temp);
        NL((4*m+1)*8*p+1:size(NL,1),:) = NL_Temp;




    else
        [NL,NoN] = NL_Extender(NL,p,m,inc);

    end

    for i = 1:m*(1+inc) 
        
        for j = 1:(4*p)
            
            if j == 1
                
                EL(2*(i-1)*(2*p)+j,1) = 4*(i-1)*4*p + j;
                EL(2*(i-1)*(2*p)+j,2) = EL(2*(i-1)*(2*p)+j,1) + 2;
                EL(2*(i-1)*(2*p)+j,4) = EL(2*(i-1)*(2*p)+j,1) + 16*p;
                EL(2*(i-1)*(2*p)+j,3) = EL(2*(i-1)*(2*p)+j,4) + 2;
                EL(2*(i-1)*(2*p)+j,5) = EL(2*(i-1)*(2*p)+j,1) + 1;
                EL(2*(i-1)*(2*p)+j,6) = EL(2*(i-1)*(2*p)+j,2) + 8*p;
                EL(2*(i-1)*(2*p)+j,7) = EL(2*(i-1)*(2*p)+j,4) + 1;
                EL(2*(i-1)*(2*p)+j,8) = EL(2*(i-1)*(2*p)+j,1) + 8*p;
                EL(2*(i-1)*(2*p)+j,9) = EL(2*(i-1)*(2*p)+j,8) + 1;
                
            elseif j == (4*p)
                
                EL(2*(i-1)*(2*p)+j,1) = (8*p*((2*i)-1))-1;
                EL(2*(i-1)*(2*p)+j,2) = (8*p*((2*i)-1)) + 1 - (8*p);
                EL(2*(i-1)*(2*p)+j,3) = EL(2*(i-1)*(2*p)+j,2) + 16*p;
                EL(2*(i-1)*(2*p)+j,4) = EL(2*(i-1)*(2*p)+j,1) + 16*p;
                EL(2*(i-1)*(2*p)+j,5) = (EL(2*(i-1)*(2*p)+j,1)) + 1;
                EL(2*(i-1)*(2*p)+j,6) = (EL(2*(i-1)*(2*p)+j,5)) + 1;
                EL(2*(i-1)*(2*p)+j,7) = (EL(2*(i-1)*(2*p)+j,4)) + 1;
                EL(2*(i-1)*(2*p)+j,8) = (EL(2*(i-1)*(2*p)+j,1)) + 8*p;
                EL(2*(i-1)*(2*p)+j,9) = (EL(2*(i-1)*(2*p)+j,8)) + 1;
                
            else
                
                EL(2*(i-1)*(2*p)+j,1) = EL(2*(i-1)*(2*p)+j-1,2);
                EL(2*(i-1)*(2*p)+j,4) = EL(2*(i-1)*(2*p)+j-1,3);
                EL(2*(i-1)*(2*p)+j,2) = EL(2*(i-1)*(2*p)+j,1) + 2;
                EL(2*(i-1)*(2*p)+j,3) = EL(2*(i-1)*(2*p)+j,4) + 2;
                EL(2*(i-1)*(2*p)+j,5) = EL(2*(i-1)*(2*p)+j,1) + 1;
                EL(2*(i-1)*(2*p)+j,6) = EL(2*(i-1)*(2*p)+j,2) + 8*p;
                EL(2*(i-1)*(2*p)+j,7) = EL(2*(i-1)*(2*p)+j,4) + 1;
                EL(2*(i-1)*(2*p)+j,8) = EL(2*(i-1)*(2*p)+j,1) + 8*p;
                EL(2*(i-1)*(2*p)+j,9) = EL(2*(i-1)*(2*p)+j,8) + 1;
                
            end
        end
    end




    if inc
        n=8*p*m; % Number of elements up to now
        NN = 32*p*m; % Node number of the down-left corner 512
        H = (4*m+1)*8*p; %544
        
        for i = 1:p
            for j = 1:p
                
                if (i == 1)
                    if j == 1
                        EL((i-1)*p+j+n,1) = NN+1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = NN+3;
                        EL((i-1)*p+j+n,4) = H-1;
                        EL((i-1)*p+j+n,3) = H+(p-1)^2;
                        EL((i-1)*p+j+n,5) = NN+2;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = H+2;
                        EL((i-1)*p+j+n,8) = H;
                        EL((i-1)*p+j+n,9) = H+1;

                    elseif j == p
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n,1)-2;
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,4)-2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,2)+1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,4)-1;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,1)-1;
                        EL((i-1)*p+j+n,9) = EL((i-1)*p+j+n,5)+1;


                    else
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n,1)-2;
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2)+2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,5)+2;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,4)+1;
                        EL((i-1)*p+j+n,9) = EL((i-1)*p+j+n,7)-1;
                    end

                elseif (i == p)
                    if j == 1
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j+n-p,2);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j+n,1)+2;
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n-p,3);
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2)+2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j+n,1)+1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n-p,3)+(p-1)^2-2;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n-p,6);
                        EL((i-1)*p+j+n,9) = EL((i-1)*p+j+n,7)-1;

                    elseif j == p
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n-p,3);
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2)+2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,2)+1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,4)-1;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,1)+1;
                        EL((i-1)*p+j+n,9) = EL((i-1)*p+j+n,5)+1;


                    else
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n-p,3);
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2)+2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,5)+2;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,1)+1;
                        EL((i-1)*p+j+n,9) = EL((i-1)*p+j+n,5)+1;
                    end


                else
                    if j == 1
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j+n-p,2);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j+n,1)+2;
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n-p,3);
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n-p,3)+2*((p-1)^2-2);
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j+n,1)+1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n-p,3)+(p-1)^2-2;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n-p,6);
                        EL((i-1)*p+j+n,9) = EL((i-1)*p+j+n,7)-1;

                    elseif j == p
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n-p,3);
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,4)-2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,2)+1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,4)-1;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,1)+1;
                        EL((i-1)*p+j+n,9) = EL((i-1)*p+j+n,5)+1;
                    else
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n,1)+2;
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2)+2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,5)+2;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,4)+1;
                        EL((i-1)*p+j+n,9) = EL((i-1)*p+j+n,7)-1;
                    end

                end
            end
        end
        
    end
end


if isequal(element_type, 'D2QU8N')

    if inc
        NL_Temp = NL(1:4*p*(2*m+1),:);
        [NL_Temp,NoN] = NL_Extender(NL_Temp,p,m,inc);
        NL_Temp (8*p*(4*m+1):(8*p*(4*m+1))+(p-1)^2,:) = NL(4*p*(2*m+1):8*p*m+(p+1)^2,:);
        NL = NL_Temp;
        NL(8*p*(4*m+1),:) = (NL(8*p*(4*m+1)-1,:) + NL(8*p*(4*m+1)-(8*p)+1,:))/2;

        for i = 1:((p-1)*4) 
            if (i < p) 
                NL(8*p*(4*m+1)+(p-1)^2+i,:) = (NL(8*p*(4*m+1)-(8*p)+2*i+1,:) + NL(8*p*(4*m+1)+i,:))/2;
            elseif (i < 3*(p-1)+1) && (i>p-1)
                NL(8*p*(4*m+1)+(p-1)^2+i,:) = (NL(8*p*(4*m+1)+i-p+1,:) + NL(8*p*(4*m+1)+i,:))/2;
            else
                NL(8*p*(4*m+1)+(p-1)^2+i,:) = (NL(8*p*(4*m+1)-(p-1)+i,:) + NL(8*p*(4*m+1)-((p-1)*2)-1-(2*(i-(p-1)^2)),:))/2;
            end
        end
        for i = 1:p 
            for j = 1:(2*p-1)
                NL(8*p*(4*m+1)+(p-1)^2+((2*p-1)*(i-1))+j+((p-1)*4),1) = ((8-((2*i)-1))*(NL(8*p*(4*m+1)-j+1,1)) + (((2*i)-1)*NL(8*p*(4*m+1)-8*p+1+2*p+j,1)))/8;
                NL(8*p*(4*m+1)+(p-1)^2+((2*p-1)*(i-1))+j+((p-1)*4),2) = (NL(8*p*(4*m+1)-j+1,2) + NL(8*p*(4*m+1)-8*p+1+2*p+j,2))/2;
            end
           
        end

        NL_Temp = NL((4*m+1)*8*p+1:size(NL,1),:);
        NL_Temp = sortrows(NL_Temp);
        NL((4*m+1)*8*p+1:size(NL,1),:) = NL_Temp;




    else
        [NL,NoN] = NL_Extender(NL,p,m,inc);

    end

    for i = 1:m*(1+inc) 
        
        for j = 1:(4*p)
            
            if j == 1
                
                EL(2*(i-1)*(2*p)+j,1) = 4*(i-1)*4*p + j;
                EL(2*(i-1)*(2*p)+j,2) = EL(2*(i-1)*(2*p)+j,1) + 2;
                EL(2*(i-1)*(2*p)+j,4) = EL(2*(i-1)*(2*p)+j,1) + 16*p;
                EL(2*(i-1)*(2*p)+j,3) = EL(2*(i-1)*(2*p)+j,4) + 2;
                EL(2*(i-1)*(2*p)+j,5) = EL(2*(i-1)*(2*p)+j,1) + 1;
                EL(2*(i-1)*(2*p)+j,6) = EL(2*(i-1)*(2*p)+j,2) + 8*p;
                EL(2*(i-1)*(2*p)+j,7) = EL(2*(i-1)*(2*p)+j,4) + 1;
                EL(2*(i-1)*(2*p)+j,8) = EL(2*(i-1)*(2*p)+j,1) + 8*p;
                
                
            elseif j == (4*p)
                
                EL(2*(i-1)*(2*p)+j,1) = (8*p*((2*i)-1))-1;
                EL(2*(i-1)*(2*p)+j,2) = (8*p*((2*i)-1)) + 1 - (8*p);
                EL(2*(i-1)*(2*p)+j,3) = EL(2*(i-1)*(2*p)+j,2) + 16*p;
                EL(2*(i-1)*(2*p)+j,4) = EL(2*(i-1)*(2*p)+j,1) + 16*p;
                EL(2*(i-1)*(2*p)+j,5) = (EL(2*(i-1)*(2*p)+j,1)) + 1;
                EL(2*(i-1)*(2*p)+j,6) = (EL(2*(i-1)*(2*p)+j,5)) + 1;
                EL(2*(i-1)*(2*p)+j,7) = (EL(2*(i-1)*(2*p)+j,4)) + 1;
                EL(2*(i-1)*(2*p)+j,8) = (EL(2*(i-1)*(2*p)+j,1)) + 8*p;
                
                
            else
                
                EL(2*(i-1)*(2*p)+j,1) = EL(2*(i-1)*(2*p)+j-1,2);
                EL(2*(i-1)*(2*p)+j,4) = EL(2*(i-1)*(2*p)+j-1,3);
                EL(2*(i-1)*(2*p)+j,2) = EL(2*(i-1)*(2*p)+j,1) + 2;
                EL(2*(i-1)*(2*p)+j,3) = EL(2*(i-1)*(2*p)+j,4) + 2;
                EL(2*(i-1)*(2*p)+j,5) = EL(2*(i-1)*(2*p)+j,1) + 1;
                EL(2*(i-1)*(2*p)+j,6) = EL(2*(i-1)*(2*p)+j,2) + 8*p;
                EL(2*(i-1)*(2*p)+j,7) = EL(2*(i-1)*(2*p)+j,4) + 1;
                EL(2*(i-1)*(2*p)+j,8) = EL(2*(i-1)*(2*p)+j,1) + 8*p;
               
                
            end
        end
    end
    
    
    
    
    if inc
        n=8*p*m; % Number of elements up to now
        NN = 32*p*m; % Node number of the down-left corner 512
        H = (4*m+1)*8*p; %544
        
        for i = 1:p
            for j = 1:p
                
                if (i == 1)
                    if j == 1
                        EL((i-1)*p+j+n,1) = NN+1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = NN+3;
                        EL((i-1)*p+j+n,4) = H-1;
                        EL((i-1)*p+j+n,3) = H+(p-1)^2;
                        EL((i-1)*p+j+n,5) = NN+2;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = H+2;
                        EL((i-1)*p+j+n,8) = H;
                        

                    elseif j == p
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n,1)-2;
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,4)-2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,2)+1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,4)-1;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,1)-1;
                        


                    else
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n,1)-2;
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2)+2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,5)+2;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,4)+1;
                        
                    end

                elseif (i == p)
                    if j == 1
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j+n-p,2);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j+n,1)+2;
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n-p,3);
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2)+2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j+n,1)+1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n-p,3)+(p-1)^2-2;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n-p,6);
                        

                    elseif j == p
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n-p,3);
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2)+2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,2)+1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,4)-1;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,1)+1;
                        


                    else
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n-p,3);
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2)+2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,5)+2;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,1)+1;
                        
                    end


                else
                    if j == 1
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j+n-p,2);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j+n,1)+2;
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n-p,3);
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n-p,3)+2*((p-1)^2-2);
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j+n,1)+1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n-p,3)+(p-1)^2-2;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n-p,6);
                        

                    elseif j == p
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n-p,3);
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,4)-2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,2)+1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,4)-1;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,1)+1;
                        
                    else
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n,1)+2;
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2)+2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,5)+2;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,4)+1;
                        
                    end

                end
            end
        end
        
    end
end


if isequal(element_type, 'D2TR6N')
    if inc
        NL_Temp = NL(1:4*p*(2*m+1),:);
        [NL_Temp,NoN] = NL_Extender(NL_Temp,p,m,inc);
        NL_Temp (8*p*(4*m+1):(8*p*(4*m+1))+(p-1)^2,:) = NL(4*p*(2*m+1):8*p*m+(p+1)^2,:);
        NL = NL_Temp;
        NL(8*p*(4*m+1),:) = (NL(8*p*(4*m+1)-1,:) + NL(8*p*(4*m+1)-(8*p)+1,:))/2;

        for i = 1:((p-1)*4) 
            if (i < p) 
                NL(8*p*(4*m+1)+(p-1)^2+i,:) = (NL(8*p*(4*m+1)-(8*p)+2*i+1,:) + NL(8*p*(4*m+1)+i,:))/2;
            elseif (i < 3*(p-1)+1) && (i>p-1)
                NL(8*p*(4*m+1)+(p-1)^2+i,:) = (NL(8*p*(4*m+1)+i-p+1,:) + NL(8*p*(4*m+1)+i,:))/2;
            else
                NL(8*p*(4*m+1)+(p-1)^2+i,:) = (NL(8*p*(4*m+1)-(p-1)+i,:) + NL(8*p*(4*m+1)-((p-1)*2)-1-(2*(i-(p-1)^2)),:))/2;
            end
        end
        for i = 1:p 
            for j = 1:(2*p-1)
                NL(8*p*(4*m+1)+(p-1)^2+((2*p-1)*(i-1))+j+((p-1)*4),1) = ((8-((2*i)-1))*(NL(8*p*(4*m+1)-j+1,1)) + (((2*i)-1)*NL(8*p*(4*m+1)-8*p+1+2*p+j,1)))/8;
                NL(8*p*(4*m+1)+(p-1)^2+((2*p-1)*(i-1))+j+((p-1)*4),2) = (NL(8*p*(4*m+1)-j+1,2) + NL(8*p*(4*m+1)-8*p+1+2*p+j,2))/2;
            end
           
        end

        NL_Temp = NL((4*m+1)*8*p+1:size(NL,1),:);
        NL_Temp = sortrows(NL_Temp);
        NL((4*m+1)*8*p+1:size(NL,1),:) = NL_Temp;




    else
        [NL,NoN] = NL_Extender(NL,p,m,inc);

    end

    for i = 1:m*(1+inc) 
        
        for j = 1:(4*p)
            
            if j == 1
                
                EL(2*(i-1)*(2*p)+j,1) = 4*(i-1)*4*p + j;
                EL(2*(i-1)*(2*p)+j,2) = EL(2*(i-1)*(2*p)+j,1) + 2;
                EL(2*(i-1)*(2*p)+j,4) = EL(2*(i-1)*(2*p)+j,1) + 16*p;
                EL(2*(i-1)*(2*p)+j,3) = EL(2*(i-1)*(2*p)+j,4) + 2;
                EL(2*(i-1)*(2*p)+j,5) = EL(2*(i-1)*(2*p)+j,1) + 1;
                EL(2*(i-1)*(2*p)+j,6) = EL(2*(i-1)*(2*p)+j,2) + 8*p;
                EL(2*(i-1)*(2*p)+j,7) = EL(2*(i-1)*(2*p)+j,4) + 1;
                EL(2*(i-1)*(2*p)+j,8) = EL(2*(i-1)*(2*p)+j,1) + 8*p;
                EL(2*(i-1)*(2*p)+j,9) = EL(2*(i-1)*(2*p)+j,8) + 1;
                
            elseif j == (4*p)
                
                EL(2*(i-1)*(2*p)+j,1) = (8*p*((2*i)-1))-1;
                EL(2*(i-1)*(2*p)+j,2) = (8*p*((2*i)-1)) + 1 - (8*p);
                EL(2*(i-1)*(2*p)+j,3) = EL(2*(i-1)*(2*p)+j,2) + 16*p;
                EL(2*(i-1)*(2*p)+j,4) = EL(2*(i-1)*(2*p)+j,1) + 16*p;
                EL(2*(i-1)*(2*p)+j,5) = (EL(2*(i-1)*(2*p)+j,1)) + 1;
                EL(2*(i-1)*(2*p)+j,6) = (EL(2*(i-1)*(2*p)+j,5)) + 1;
                EL(2*(i-1)*(2*p)+j,7) = (EL(2*(i-1)*(2*p)+j,4)) + 1;
                EL(2*(i-1)*(2*p)+j,8) = (EL(2*(i-1)*(2*p)+j,1)) + 8*p;
                EL(2*(i-1)*(2*p)+j,9) = (EL(2*(i-1)*(2*p)+j,8)) + 1;
                
            else
                
                EL(2*(i-1)*(2*p)+j,1) = EL(2*(i-1)*(2*p)+j-1,2);
                EL(2*(i-1)*(2*p)+j,4) = EL(2*(i-1)*(2*p)+j-1,3);
                EL(2*(i-1)*(2*p)+j,2) = EL(2*(i-1)*(2*p)+j,1) + 2;
                EL(2*(i-1)*(2*p)+j,3) = EL(2*(i-1)*(2*p)+j,4) + 2;
                EL(2*(i-1)*(2*p)+j,5) = EL(2*(i-1)*(2*p)+j,1) + 1;
                EL(2*(i-1)*(2*p)+j,6) = EL(2*(i-1)*(2*p)+j,2) + 8*p;
                EL(2*(i-1)*(2*p)+j,7) = EL(2*(i-1)*(2*p)+j,4) + 1;
                EL(2*(i-1)*(2*p)+j,8) = EL(2*(i-1)*(2*p)+j,1) + 8*p;
                EL(2*(i-1)*(2*p)+j,9) = EL(2*(i-1)*(2*p)+j,8) + 1;
                
            end
        end
    end




    if inc
        n=8*p*m; % Number of elements up to now
        NN = 32*p*m; % Node number of the down-left corner 512
        H = (4*m+1)*8*p; %544
        
        for i = 1:p
            for j = 1:p
                
                if (i == 1)
                    if j == 1
                        EL((i-1)*p+j+n,1) = NN+1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = NN+3;
                        EL((i-1)*p+j+n,4) = H-1;
                        EL((i-1)*p+j+n,3) = H+(p-1)^2;
                        EL((i-1)*p+j+n,5) = NN+2;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = H+2;
                        EL((i-1)*p+j+n,8) = H;
                        EL((i-1)*p+j+n,9) = H+1;

                    elseif j == p
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n,1)-2;
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,4)-2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,2)+1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,4)-1;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,1)-1;
                        EL((i-1)*p+j+n,9) = EL((i-1)*p+j+n,5)+1;


                    else
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n,1)-2;
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2)+2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,5)+2;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,4)+1;
                        EL((i-1)*p+j+n,9) = EL((i-1)*p+j+n,7)-1;
                    end

                elseif (i == p)
                    if j == 1
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j+n-p,2);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j+n,1)+2;
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n-p,3);
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2)+2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j+n,1)+1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n-p,3)+(p-1)^2-2;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n-p,6);
                        EL((i-1)*p+j+n,9) = EL((i-1)*p+j+n,7)-1;

                    elseif j == p
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n-p,3);
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2)+2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,2)+1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,4)-1;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,1)+1;
                        EL((i-1)*p+j+n,9) = EL((i-1)*p+j+n,5)+1;


                    else
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n-p,3);
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2)+2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,5)+2;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,1)+1;
                        EL((i-1)*p+j+n,9) = EL((i-1)*p+j+n,5)+1;
                    end


                else
                    if j == 1
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j+n-p,2);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j+n,1)+2;
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n-p,3);
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n-p,3)+2*((p-1)^2-2);
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j+n,1)+1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n-p,3)+(p-1)^2-2;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n-p,6);
                        EL((i-1)*p+j+n,9) = EL((i-1)*p+j+n,7)-1;

                    elseif j == p
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n-p,3);
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,4)-2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,2)+1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,4)-1;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,1)+1;
                        EL((i-1)*p+j+n,9) = EL((i-1)*p+j+n,5)+1;
                    else
                        EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,4);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,2) = EL((i-1)*p+j-1+n,3);
                        EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n,1)+2;
                        EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2)+2;
                        EL((i-1)*p+j+n,5) = EL((i-1)*p+j-1+n,7);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        NN;
                        EL((i-1)*p+j+n,6) = EL((i-1)*p+j+n,3)-1;
                        EL((i-1)*p+j+n,7) = EL((i-1)*p+j+n,5)+2;
                        EL((i-1)*p+j+n,8) = EL((i-1)*p+j+n,4)+1;
                        EL((i-1)*p+j+n,9) = EL((i-1)*p+j+n,7)-1;
                    end

                end
            end
        end
        
    end

    NPE_new = 6;
    EL_new = zeros(NoE, NPE_new);
    
    for i = 1:NoE
        
        EL_new(2*(i-1)+1,1) = EL(i,1); 
        EL_new(2*(i-1)+1,2) = EL(i,2);
        EL_new(2*(i-1)+1,3) = EL(i,3);
        EL_new(2*(i-1)+1,4) = EL(i,5); 
        EL_new(2*(i-1)+1,5) = EL(i,6);
        EL_new(2*(i-1)+1,6) = EL(i,9);
        
        EL_new(2*(i-1)+2,1) = EL(i,1); 
        EL_new(2*(i-1)+2,2) = EL(i,3);
        EL_new(2*(i-1)+2,3) = EL(i,4);
        EL_new(2*(i-1)+2,4) = EL(i,9); 
        EL_new(2*(i-1)+2,5) = EL(i,7);
        EL_new(2*(i-1)+2,6) = EL(i,8);
        
    end
    
    EL = EL_new;
    
end



if isequal(element_type,'D2TR3N')

    for i = 1:m*(1+inc) 
        
        for j = 1:4*p
            
            if j == 1
                
                EL((i-1)*(4*p)+j,1) = (i-1)*(4*p) + 1;
                EL((i-1)*(4*p)+j,2) = EL((i-1)*(4*p)+j,1) + 1;
                EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j,1) + 4*p;
                EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j,4) + 1;
                
            elseif j == 4*p
                
                EL((i-1)*(4*p)+j,1) = (i)*(4*p);
                EL((i-1)*(4*p)+j,2) = (i-1)*(4*p) + 1;
                EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j,1) + 1;
                EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j,1) + 4*p;
                
            else
                
                EL((i-1)*(4*p)+j,1) = EL((i-1)*(4*p)+j-1,2);
                EL((i-1)*(4*p)+j,4) = EL((i-1)*(4*p)+j-1,3);
                EL((i-1)*(4*p)+j,3) = EL((i-1)*(4*p)+j,4) + 1;
                EL((i-1)*(4*p)+j,2) = EL((i-1)*(4*p)+j,1) + 1;
                
            end
        end
    end
    
    
    
    if inc
        n=8*p*m; % Number of elements up to now
        NN = 2*(p+1)*(2*m+1) + 2*(p-1)*(2*m+1) - p*4 + 1; % Node number of the down-left corner
        
        for i = 1:p
            for j = 1:p
                
                if (i == 1)&&(j == 1)
                    EL((i-1)*p+j+n,1) =                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         NN;
                    EL((i-1)*p+j+n,2) = NN+1;
                    EL((i-1)*p+j+n,4) = NN+4*p-1;
                    EL((i-1)*p+j+n,3) = NN+4*p;
                elseif (i == 1)&&(j == p)
                    EL((i-1)*p+j+n,1) = NN+p-1;
                    EL((i-1)*p+j+n,2) = NN+p;
                    EL((i-1)*p+j+n,4) = NN+4*p+(p-2);
                    EL((i-1)*p+j+n,3) = NN+p+1;
                elseif (j == 1)&&(i == p)
                    EL((i-1)*p+j+n,1) = NN+3*p+1;
                    EL((i-1)*p+j+n,2) = NoN-(p-2);
                    EL((i-1)*p+j+n,4) = NN+3*p;
                    EL((i-1)*p+j+n,3) = NN+3*p-1;
                elseif (j == p)&&(i == p) 
                    EL((i-1)*p+j+n,1) = NoN;
                    EL((i-1)*p+j+n,2) = NN+2*p-1;
                    EL((i-1)*p+j+n,4) = NN+2*p+1;
                    EL((i-1)*p+j+n,3) = NN+2*p;
                
                elseif i==1
                    EL((i-1)*p+j+n,1) = NN+j-1;
                    EL((i-1)*p+j+n,2) = NN+j;
                    EL((i-1)*p+j+n,4) = NN+4*p-2+j;
                    EL((i-1)*p+j+n,3) = NN+4*p-1+j;
                elseif i==p
                    EL((i-1)*p+j+n,1) = NoN-(p-2)-2+j;
                    EL((i-1)*p+j+n,2) = NoN-(p-2)-1+j;
                    EL((i-1)*p+j+n,4) = NN+3*p+1-j;
                    EL((i-1)*p+j+n,3) = NN+3*p-j;
                elseif j == 1
                    EL((i-1)*p+j+n,1) = NN+4*p+1-i;
                    EL((i-1)*p+j+n,2) = NN+4*p+(p-1)*(i-2);
                    EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n,1) - 1;
                    EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2) + (p-1);
                elseif j == p
                    EL((i-1)*p+j+n,1) = NN+4*p+(p-2)+(p-1)*(i-2);
                    EL((i-1)*p+j+n,2) = NN+p-1+i;
                    EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n,1) + (p-1);
                    EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,2) + 1;
                else
                    
                    EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,2);
                    EL((i-1)*p+j+n,4) = EL((i-1)*p+j-1+n,3);
                    EL((i-1)*p+j+n,2) = EL((i-1)*p+j+n,1) + 1;
                    EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,4) + 1;
    %                 
    %                 if j == 1
    %                     
    %                     EL((i-1)*p+j+n,1) = (i-1)*(p+1)+j + NN;
    %                     EL((i-1)*p+j+n,2) = EL((i-1)*p+j+n,1) + 1;
    %                     EL((i-1)*p+j+n,4) = EL((i-1)*p+j+n,1) + (p+1);
    %                     EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,4) + 1;
    %                     
    %                 else
    %                     
    %                     EL((i-1)*p+j+n,1) = EL((i-1)*p+j-1+n,2);
    %                     EL((i-1)*p+j+n,4) = EL((i-1)*p+j-1+n,3);
    %                     EL((i-1)*p+j+n,2) = EL((i-1)*p+j+n,1) + 1;
    %                     EL((i-1)*p+j+n,3) = EL((i-1)*p+j+n,4) + 1;
    %                     
    %                 end
                end
                
            end
        end
        
    end

    
    NPE_new = 3;
    EL_new = zeros(2*NoE, NPE_new);
    
    for i = 1:NoE
        
        EL_new(2*(i-1)+1,1) = EL(i,1); % Lower triangle
        EL_new(2*(i-1)+1,2) = EL(i,2);
        EL_new(2*(i-1)+1,3) = EL(i,3);
        
        EL_new(2*(i-1)+2,1) = EL(i,1); % Upper triangle
        EL_new(2*(i-1)+2,2) = EL(i,3);
        EL_new(2*(i-1)+2,3) = EL(i,4);
        
    end
    
    EL = EL_new;
    
end

end 





function [NL,NoN] = NL_Extender(NL,p,m,inc)

    NL_new = zeros(2*size(NL,1),2);

    for i = 1:size(NL,1)
        
        if i == size(NL,1)
            NL_new(2*i-1,:) = NL(i,:);
            NL_new(2*i,:) = (NL(i,:) + NL(i-(4*p)+1,:))/2;

        elseif rem(i,4*p) == 0
            NL_new(2*i-1,:) = NL(i,:);
            NL_new(2*i,:) = (NL(i,:) + NL(i-(4*p)+1,:))/2;

        else
            NL_new(2*i-1,:) = NL(i,:);
            NL_new(2*i,:) = (NL(i,:) + NL(i+1,:))/2;
        end
    end

    NL = NL_new;

    NL_new = zeros(2*size(NL,1)- 8*p,2);


    for j = 1:m*(1+inc)+1

        if j == m*(1+inc)+1
            NL_new(((8*p*(2*j-2))+1):8*p*(2*j-1),:) = NL(8*p*(j-1)+1:8*p*j,:);
        else

            NL_new(((8*p*(2*j-2))+1):8*p*(2*j-1),:) = NL(8*p*(j-1)+1:8*p*j,:);
    
            for k = 1:8*p

                if rem((2*j-1)*8*p+k,16*p) == 0

                    NL_new((2*j-1)*8*p+k,:) = (NL((j-1)*8*p+k+1,:) + NL((j-1)*8*p+k-1,:))/2;

                else
    
                    if rem((2*j-1)*8*p+k,2) == 0
    
                        NL_new((2*j-1)*8*p+k,:) = (NL(8*p*(j-1)+k-1,:) + NL(8*p*(j)+k+1,:))/2;
    
                    else
        
                        NL_new((2*j-1)*8*p+k,:) = (NL(8*p*(j-1)+k,:) + NL(8*p*(j)+k,:))/2;
                    end
                end
            end
        end
    end


    NL = NL_new;
    NoN = size(NL,1);

end
