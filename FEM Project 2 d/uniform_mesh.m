function [NL, EL] = uniform_mesh(d1, d2, p, m, element_type)

PD = 2;

q = [0 0; d1 0; 0 d2; d1 d2];

NoN = (p+1)*(m+1);

if isequal(element_type, 'D2TR3N')
    NPE = 3;
    NoE = p*m;
elseif isequal(element_type, 'D2TR6N')
    NPE = 6;
    NoE = (p/2)*(m/2);
elseif isequal(element_type, 'D2QU8N')
    NPE = 8;
    NoE = (p/2)*(m/2);
elseif isequal(element_type, 'D2QU9N')
    NPE = 9;
    NoE = (p/2)*(m/2);
elseif isequal(element_type, 'D2QU4N')
    NPE = 4;
    NoE = p*m;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Nodes   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NL = zeros(NoN, PD);

a = d1/p;
b = d2/m;

% if isequal(element_type, 'D2QU8N') 
%     n=1;
%     NoN = ((p+1)*((m/2)+1))+((p-1)*(m/2));
%     NL = zeros(NoN, PD);
% 
%     for i = 1:m+1
% 
%         if rem(i,2) == 1
%         
%            for j = 1:p+1
%                
%                NL(n,1) = q(1,1) + (j-1)*a;
%                NL(n,2) = q(1,2) + (i-1)*b;
%                n = n+1;
%                
%            end
% 
%         else
% 
%             for j = 1:2:p+1
%                
%                NL(n,1) = q(1,1) + (j-1)*a;
%                NL(n,2) = q(1,2) + (i-1)*b;
%                n = n+1;
%                
%            end
% 
%         end
%        
%     end
% 
% else
% 
%     n=1;
%     
%     for i = 1:m+1
%         
%        for j = 1:p+1
%            
%            NL(n,1) = q(1,1) + (j-1)*a;
%            NL(n,2) = q(1,2) + (i-1)*b;
%            n = n+1;
%            
%        end
%        
%     end
% end
    n=1;

    for i = 1:m+1
        
       for j = 1:p+1
           
           NL(n,1) = q(1,1) + (j-1)*a;
           NL(n,2) = q(1,2) + (i-1)*b;
           n = n+1;
           
       end
       
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Elements   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       
EL = zeros(NoE, NPE);

if isequal(element_type, 'D2QU4N')
    for i = 1:m
        
        for j = 1:p
            
            if j == 1
                
                EL((i-1)*p+j,1) = (i-1)*(p+1)+j;
                EL((i-1)*p+j,2) = EL((i-1)*p+j,1) + 1;
                EL((i-1)*p+j,4) = EL((i-1)*p+j,1) + (p+1);
                EL((i-1)*p+j,3) = EL((i-1)*p+j,4) + 1;
                
            else
                
                EL((i-1)*p+j,1) = EL((i-1)*p+j-1,2);
                EL((i-1)*p+j,4) = EL((i-1)*p+j-1,3);
                EL((i-1)*p+j,2) = EL((i-1)*p+j,1) + 1;
                EL((i-1)*p+j,3) = EL((i-1)*p+j,4) + 1;
                
            end
            
        end
        
    end

elseif isequal(element_type, 'D2QU9N')
    for i = 1:(m/2)
        
        for j = 1:(p/2)
            
            if j == 1
                
                EL((i-1)*(p/2)+j,1) = (i-1)*(p+1)*2 + j;
                EL((i-1)*(p/2)+j,2) = EL((i-1)*(p/2)+j,1) + 2;
                EL((i-1)*(p/2)+j,4) = EL((i-1)*(p/2)+j,1) + (2*p+2);
                EL((i-1)*(p/2)+j,3) = EL((i-1)*(p/2)+j,4) + 2;
                EL((i-1)*(p/2)+j,5) = EL((i-1)*(p/2)+j,1) + 1;
                EL((i-1)*(p/2)+j,6) = EL((i-1)*(p/2)+j,2) + p+1;
                EL((i-1)*(p/2)+j,7) = EL((i-1)*(p/2)+j,4) + 1;
                EL((i-1)*(p/2)+j,8) = EL((i-1)*(p/2)+j,1) + p+1;
                EL((i-1)*(p/2)+j,9) = EL((i-1)*(p/2)+j,8) + 1;
                
            else
                
                EL((i-1)*(p/2)+j,1) = EL((i-1)*(p/2)+j-1,2);
                EL((i-1)*(p/2)+j,4) = EL((i-1)*(p/2)+j-1,3);
                EL((i-1)*(p/2)+j,2) = EL((i-1)*(p/2)+j,1) + 2;
                EL((i-1)*(p/2)+j,3) = EL((i-1)*(p/2)+j,4) + 2;
                EL((i-1)*(p/2)+j,5) = EL((i-1)*(p/2)+j,1) + 1;
                EL((i-1)*(p/2)+j,6) = EL((i-1)*(p/2)+j,2) + p+1;
                EL((i-1)*(p/2)+j,7) = EL((i-1)*(p/2)+j,4) + 1;
                EL((i-1)*(p/2)+j,8) = EL((i-1)*(p/2)+j,1) + p+1;
                EL((i-1)*(p/2)+j,9) = EL((i-1)*(p/2)+j,8) + 1;
                
            end
            
        end
        
    end

elseif isequal(element_type, 'D2QU8N')
%     for i = 1:(m/2)
%         
%         for j = 1:(p/2)
%             
%             if j == 1
%                 
%                 EL((i-1)*(p/2)+j,1) = (i-1)*p*2 + j;
%                 EL((i-1)*(p/2)+j,2) = EL((i-1)*(p/2)+j,1) + 2;
%                 EL((i-1)*(p/2)+j,4) = EL((i-1)*(p/2)+j,1) + (2*p);
%                 EL((i-1)*(p/2)+j,3) = EL((i-1)*(p/2)+j,4) + 2;
%                 EL((i-1)*(p/2)+j,5) = EL((i-1)*(p/2)+j,1) + 1;
%                 EL((i-1)*(p/2)+j,6) = EL((i-1)*(p/2)+j,2) + p-j+1;
%                 EL((i-1)*(p/2)+j,7) = EL((i-1)*(p/2)+j,4) + 1;
%                 EL((i-1)*(p/2)+j,8) = EL((i-1)*(p/2)+j,1) + p-j+2;
%                 
%                 
%             else
%                 
%                 EL((i-1)*(p/2)+j,1) = EL((i-1)*(p/2)+j-1,2);
%                 EL((i-1)*(p/2)+j,4) = EL((i-1)*(p/2)+j-1,3);
%                 EL((i-1)*(p/2)+j,2) = EL((i-1)*(p/2)+j,1) + 2;
%                 EL((i-1)*(p/2)+j,3) = EL((i-1)*(p/2)+j,4) + 2;
%                 EL((i-1)*(p/2)+j,5) = EL((i-1)*(p/2)+j,1) + 1;
%                 EL((i-1)*(p/2)+j,6) = EL((i-1)*(p/2)+j,2) + p-j+1;
%                 EL((i-1)*(p/2)+j,7) = EL((i-1)*(p/2)+j,4) + 1;
%                 EL((i-1)*(p/2)+j,8) = EL((i-1)*(p/2)+j,1) + p-j+2;
%                 
%                 
%             end
%             
%         end
%         
%     end

    for i = 1:(m/2)
        
        for j = 1:(p/2)
            
            if j == 1
                
                EL((i-1)*(p/2)+j,1) = (i-1)*(p+1)*2 + j;
                EL((i-1)*(p/2)+j,2) = EL((i-1)*(p/2)+j,1) + 2;
                EL((i-1)*(p/2)+j,4) = EL((i-1)*(p/2)+j,1) + (2*p+2);
                EL((i-1)*(p/2)+j,3) = EL((i-1)*(p/2)+j,4) + 2;
                EL((i-1)*(p/2)+j,5) = EL((i-1)*(p/2)+j,1) + 1;
                EL((i-1)*(p/2)+j,6) = EL((i-1)*(p/2)+j,2) + p+1;
                EL((i-1)*(p/2)+j,7) = EL((i-1)*(p/2)+j,4) + 1;
                EL((i-1)*(p/2)+j,8) = EL((i-1)*(p/2)+j,1) + p+1;
                
                
            else
                
                EL((i-1)*(p/2)+j,1) = EL((i-1)*(p/2)+j-1,2);
                EL((i-1)*(p/2)+j,4) = EL((i-1)*(p/2)+j-1,3);
                EL((i-1)*(p/2)+j,2) = EL((i-1)*(p/2)+j,1) + 2;
                EL((i-1)*(p/2)+j,3) = EL((i-1)*(p/2)+j,4) + 2;
                EL((i-1)*(p/2)+j,5) = EL((i-1)*(p/2)+j,1) + 1;
                EL((i-1)*(p/2)+j,6) = EL((i-1)*(p/2)+j,2) + p+1;
                EL((i-1)*(p/2)+j,7) = EL((i-1)*(p/2)+j,4) + 1;
                EL((i-1)*(p/2)+j,8) = EL((i-1)*(p/2)+j,1) + p+1;
                
                
            end
            
        end
        
    end



elseif isequal(element_type, 'D2TR3N')

    for i = 1:m
        
        for j = 1:p
            
            if j == 1
                
                EL((i-1)*p+j,1) = (i-1)*(p+1)+j;
                EL((i-1)*p+j,2) = EL((i-1)*p+j,1) + 1;
                EL((i-1)*p+j,4) = EL((i-1)*p+j,1) + (p+1);
                EL((i-1)*p+j,3) = EL((i-1)*p+j,4) + 1;
                
            else
                
                EL((i-1)*p+j,1) = EL((i-1)*p+j-1,2);
                EL((i-1)*p+j,4) = EL((i-1)*p+j-1,3);
                EL((i-1)*p+j,2) = EL((i-1)*p+j,1) + 1;
                EL((i-1)*p+j,3) = EL((i-1)*p+j,4) + 1;
                
            end
            
        end
        
    end
    
    NPE_new = 3;
    EL_new = zeros(2*NoE, NPE_new);
    
    for i = 1:NoE
        
        EL_new(2*(i-1)+1,1) = EL(i,1); 
        EL_new(2*(i-1)+1,2) = EL(i,2);
        EL_new(2*(i-1)+1,3) = EL(i,3);
        
        EL_new(2*(i-1)+2,1) = EL(i,1); 
        EL_new(2*(i-1)+2,2) = EL(i,3);
        EL_new(2*(i-1)+2,3) = EL(i,4);
        
    end
    
    EL = EL_new;

elseif isequal(element_type, 'D2TR6N')
    
    for i = 1:(m/2)
        
        for j = 1:(p/2)
            
            if j == 1
                
                EL((i-1)*(p/2)+j,1) = (i-1)*(p+1)*2 + j;
                EL((i-1)*(p/2)+j,2) = EL((i-1)*(p/2)+j,1) + 2;
                EL((i-1)*(p/2)+j,4) = EL((i-1)*(p/2)+j,1) + (2*p+2);
                EL((i-1)*(p/2)+j,3) = EL((i-1)*(p/2)+j,4) + 2;
                EL((i-1)*(p/2)+j,5) = EL((i-1)*(p/2)+j,1) + 1;
                EL((i-1)*(p/2)+j,6) = EL((i-1)*(p/2)+j,2) + p+1;
                EL((i-1)*(p/2)+j,7) = EL((i-1)*(p/2)+j,4) + 1;
                EL((i-1)*(p/2)+j,8) = EL((i-1)*(p/2)+j,1) + p+1;
                EL((i-1)*(p/2)+j,9) = EL((i-1)*(p/2)+j,8) + 1;
                
            else
                
                EL((i-1)*(p/2)+j,1) = EL((i-1)*(p/2)+j-1,2);
                EL((i-1)*(p/2)+j,4) = EL((i-1)*(p/2)+j-1,3);
                EL((i-1)*(p/2)+j,2) = EL((i-1)*(p/2)+j,1) + 2;
                EL((i-1)*(p/2)+j,3) = EL((i-1)*(p/2)+j,4) + 2;
                EL((i-1)*(p/2)+j,5) = EL((i-1)*(p/2)+j,1) + 1;
                EL((i-1)*(p/2)+j,6) = EL((i-1)*(p/2)+j,2) + p+1;
                EL((i-1)*(p/2)+j,7) = EL((i-1)*(p/2)+j,4) + 1;
                EL((i-1)*(p/2)+j,8) = EL((i-1)*(p/2)+j,1) + p+1;
                EL((i-1)*(p/2)+j,9) = EL((i-1)*(p/2)+j,8) + 1;
                
            end
            
        end
        
    end

    EL_new = zeros(NoE, NPE);
    
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

    
    
end

