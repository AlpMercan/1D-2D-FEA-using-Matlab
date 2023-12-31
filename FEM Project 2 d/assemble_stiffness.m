function [K] = assemble_stiffness(ENL, EL, NL, PR, p, m, element_type)

NoE = size(EL,1); % Number of Elements
NPE = size(EL,2); % Nodes per Element

NoN = size(NL,1);
PD = size(NL,2);

K = zeros(NoN*PD, NoN*PD);

if NPE == 3  
    
    GPE = 1;
    
elseif NPE == 4
    
    GPE = 4; 

elseif NPE == 6
    
    GPE = 3;

elseif NPE == 8
    
    GPE = 4; 

else
    GPE = 9;
    
end

for i = 1:NoE
    
    nl = EL(i,1:NPE); % Node List of each element
    
    x = zeros(NPE, PD);
    
    for j = 1:NPE
        
        x(j,:) = ENL(nl(j),1:PD);

    end
    
    if i > 4*p*m
        E = PR(2,1);
        v = PR(2,2);
        
    else
        E = PR(1,1);
        v = PR(1,2);
    end
    
    k = Stiffness(x, GPE, element_type);  % Stiffness calculates element stiffness
    
    for r = 1:NPE %d1
        
        for p = 1:PD
           
            for q = 1:NPE %d2
                
                for s = 1:PD
                    row = ENL(nl(r), p+3*PD); %Orders wrt the global DOFs
                    column = ENL(nl(q),s+3*PD);
                    
                    value = k((r-1)*PD+p, (q-1)*PD+s);
                    K(row,column) = K(row,column) + value;
                end
            end
        end
    end
end