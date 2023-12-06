function [] = post_process(NL, EL, ENL,PR,p,m,element_type, stressswitch, figuree, colormapp, edgec)

    PD = size(NL,2);
    NPE = size(EL,2);
    
    NoN = size(ENL,1);
    NoE = size(EL,1);
    
    scale = 2; % Magnifies deformation
    
    [disp, stress, strain] = element_post_process(NL, EL, ENL,PR,p,m,element_type);
    
    switch NPE
        
        case 3
            
            for i = 1:NoE
                
                nl = EL(i,1:NPE); % node numbers connected to element
                
                for j = 1:NPE
                    
                    X(j,i) = ENL(nl(j),1) + scale*ENL(nl(j),4*PD+1);
                    Y(j,i) = ENL(nl(j),2) + scale*ENL(nl(j),4*PD+2);
                    
                end
                
                for j = 1:NPE
                    stress_xx(j,i) = stress(i,:,1,1);
                end
                
                for j = 1:NPE
                    stress_xy(j,i) = stress(i,:,1,2);
                end
                
                for j = 1:NPE
                    stress_yx(j,i) = stress(i,:,2,1);
                end
                
                for j = 1:NPE
                    stress_yy(j,i) = stress(i,:,2,2);
                end
                
                
                for j = 1:NPE
                    strain_xx(j,i) = stress(i,:,1,1);
                end
                
                for j = 1:NPE
                    strain_xy(j,i) = stress(i,:,1,2);
                end
                
                for j = 1:NPE
                    strain_yx(j,i) = stress(i,:,2,1);
                end
                
                for j = 1:NPE
                    strain_yy(j,i) = stress(i,:,2,2);
                end
                
                
                for j = 1:NPE
                    disp_x(j,i) = disp(i,j,1,1);
                end
                
                for j = 1:NPE
                    disp_y(j,i) = disp(i,j,2,1);
                end
            end
        

        case 4
            
            for i = 1:NoE
                
                nl = EL(i,1:NPE); % node numbers connected to element
                
                for j = 1:NPE
                    
                    X(j,i) = ENL(nl(j),1) + scale*ENL(nl(j),4*PD+1);
                    Y(j,i) = ENL(nl(j),2) + scale*ENL(nl(j),4*PD+2);
                    
                end
                
                for j = 1:NPE
                    val = stress(i,:,1,1);
                    stress_xx(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = stress(i,:,1,2);
                    stress_xy(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = stress(i,:,2,1);
                    stress_yx(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = stress(i,:,2,2);
                    stress_yy(j,i) = val(1,j);
                end
                
                
                for j = 1:NPE
                    val = strain(i,:,1,1);
                    strain_xx(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = strain(i,:,1,2);
                    strain_xy(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = strain(i,:,2,1);
                    strain_yx(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = strain(i,:,2,2);
                    strain_yy(j,i) = val(1,j);
                end
                
                
                for j = 1:NPE
                    disp_x(j,i) = disp(i,j,1,1);
                end
                
                for j = 1:NPE
                    disp_y(j,i) = disp(i,j,2,1);
                end
            end


        case 6  

%             stress = stress(:,[2 3 1],:,:);      conseptually working but
%             not in code
% 
%             strain = strain(:,[2 3 1],:,:);

            for i = 1:NoE
                NPE =3;
                nl = EL(i,1:NPE);
                
                for j  =1:NPE 
                    
                    X(j,i) = ENL(nl(j),1) +scale*ENL(nl(j),4*PD+1);
                    Y(j,i) = ENL(nl(j),2) +scale*ENL(nl(j),4*PD+2);
                end
                
                for j = 1:NPE
                    val = stress(i,:,1,1);
                    stress_xx(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = stress(i,:,1,2);
                    stress_xy(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = stress(i,:,2,1);
                    stress_yx(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = stress(i,:,2,2);
                    stress_yy(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = strain(i,:,1,1);
                    strain_xx(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = strain(i,:,1,2);
                    strain_xy(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = strain(i,:,2,1);
                    strain_yx(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = strain(i,:,2,2);
                    strain_yy(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    disp_x(j,i) = disp(i,j,1,1);
                end
                
                for j = 1:NPE
                    disp_y(j,i) = disp(i,j,2,1);
                end
                
            end


        case 9

%                 stress = stress(:,[1 5 2 6 3 7 4 8 9],:,:);
%                 conceptually working but not in code
%                 
%                 strain = strain(:,[1 5 2 6 3 7 4 8 9],:,:);
      
            for i = 1:NoE
                nl = EL(i,1:NPE);
                
                for j  =1:NPE 
                    
                    X(j,i) = ENL(nl(j),1) + scale*ENL(nl(j),4*PD+1);
                    Y(j,i) = ENL(nl(j),2) + scale*ENL(nl(j),4*PD+2);
                end
                
                for j = 1:NPE
                    val = stress(i,:,1,1);
                    stress_xx(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = stress(i,:,1,2);
                    stress_xy(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = stress(i,:,2,1);
                    stress_yx(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = stress(i,:,2,2);
                    stress_yy(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = strain(i,:,1,1);
                    strain_xx(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = strain(i,:,1,2);
                    strain_xy(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = strain(i,:,2,1);
                    strain_yx(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    val = strain(i,:,2,2);
                    strain_yy(j,i) = val(1,j);
                end
                
                for j = 1:NPE
                    disp_x(j,i) = disp(i,j,1,1);
                end
                
                for j = 1:NPE
                    disp_y(j,i) = disp(i,j,2,1);
                end
                
            end
   
    end

    switch stressswitch

        case "Stress_xx"
                
            f1 = figuree;
            patch(figuree,X,Y,stress_xx,'FaceColor','interp',"EdgeColor", edgec);
            colormap (figuree, colormapp);
            title(figuree,'Distribution of stress xx', 'FontWeight', 'bold');
            axis (figuree, "equal");
            colorbar(figuree, "Color","w");

        case "Stress_xy"
            f2 = figuree;
            patch(figuree,X,Y,stress_xy,'FaceColor','interp',"EdgeColor", edgec);
            colormap (figuree, colormapp);
            title(figuree,'Distribution of stress xy', 'FontWeight', 'bold');
            axis (figuree, "equal");
            colorbar(figuree, "Color","w");
            
        case "Stress_yx"
            f3 = figuree;
            patch(figuree,X,Y,stress_yx,'FaceColor','interp',"EdgeColor", edgec);
            colormap (figuree, colormapp);
            title(figuree,'Distribution of stress yx', 'FontWeight', 'bold');
            axis (figuree, "equal");
            colorbar(figuree, "Color","w");

        case "Stress_yy"
            f4 = figuree;
            patch(figuree, X,Y,stress_yy,'FaceColor','interp',"EdgeColor", edgec);
            colormap (figuree, colormapp);
            title(figuree,'Distribution of stress yy', 'FontWeight', 'bold');
            axis (figuree, "equal");
            colorbar(figuree, "Color","w");

        case "Strain_xx"
            f5 = figuree;
            patch(figuree,X,Y,strain_xx,'FaceColor','interp',"EdgeColor", edgec);
            colormap (figuree, colormapp);
            title(figuree,'Distribution of strain xx', 'FontWeight', 'bold');
            axis (figuree, "equal");
            colorbar(figuree, "Color","w");

        case "Strain_xy"
            f6 = figuree;
            patch(figuree,X,Y,strain_xy,'FaceColor','interp',"EdgeColor", edgec);
            colormap (figuree, colormapp);
            title(figuree,'Distribution of strain xy', 'FontWeight', 'bold');
            axis (figuree, "equal");
            colorbar(figuree, "Color","w");

        case "Strain_yx"
            f7 = figuree;
            patch(figuree,X,Y,strain_yx,'FaceColor','interp',"EdgeColor", edgec);
            colormap (figuree, colormapp);
            title(figuree,'Distribution of strain yx', 'FontWeight', 'bold');
            axis (figuree, "equal");
            colorbar(figuree, "Color","w");

        case "Strain_yy"
            f8 = figuree;
            patch(figuree,X,Y,strain_yy,'FaceColor','interp',"EdgeColor", edgec);
            colormap (figuree, colormapp);
            title(figuree,'Distribution of strain yy', 'FontWeight', 'bold');
            axis (figuree, "equal");
            colorbar(figuree, "Color","w");

        case "Disp_x"
            f9 = figuree;
            patch(figuree,X,Y,disp_x,'FaceColor','interp',"EdgeColor", edgec);
            colormap (figuree, colormapp);
            title(figuree,'Distribution of displacement x', 'FontWeight', 'bold');
            axis (figuree, "equal");
            colorbar(figuree, "Color","w");

        case "Disp_y"
            f10 = figuree;
            patch(figuree,X,Y,disp_y,'FaceColor','interp',"EdgeColor", edgec);
            colormap (figuree, colormapp);
            title(figuree,'Distribution of displacement y', 'FontWeight', 'bold');
            axis (figuree, "equal");
            colorbar(figuree, "Color","w");
        
%         movegui(f1, 'northwest');
%         movegui(f2, 'north');
%         movegui(f3, 'north');
%         movegui(f4, 'northeast');
%         
%         movegui(f5, 'southwest');
%         movegui(f6, 'south');
%         movegui(f7, 'south');
%         movegui(f8, 'southeast');
%         
%         movegui(f9, 'west');
%         movegui(f10, 'east');
    end
end

function [disp, stress, strain] = element_post_process(NL, EL, ENL,PR,p,m,element_type)

    PD = size(NL,2);
    NPE = size(EL,2);
    
    NoN = size(ENL,1);
    NoE = size(EL,1);
    
    if NPE == 3
        GPE = 1;
    elseif NPE == 4
        GPE = 4;    
    elseif NPE == 6
        GPE = 3;
    elseif NPE == 8
        GPE = 4;
    elseif NPE == 9
        GPE = 9;    
    end
    
    disp = zeros(NoE, NPE, PD, 1);
    stress = zeros(NoE, GPE, PD, PD);
    strain = zeros(NoE, GPE, PD, PD);


if isequal(element_type, 'D2TR3N')
    NoE_inc = p*8*m;
    
elseif isequal(element_type, 'D2TR6N')
    NoE_inc = 8*p*m;
   
elseif isequal(element_type, 'D2QU8N')
    
    NoE_inc = 4*p*m;

elseif isequal(element_type, 'D2QU9N')
    NoE_inc = 4*p*m;

elseif isequal(element_type, 'D2QU4N')
    NoE_inc = 4*p*m;

end


    
    for e = 1:NoE
        
        nl = EL(e,1:NPE); % node numbers connected to each element
        
        if e > NoE_inc
            E = PR(2,1);
            v = PR(2,2);

        else
            E = PR(1,1);
            v = PR(1,2);
        end
        
        for i = 1:NPE
            for j = 1:PD
                disp(e,i,j,1) = ENL(nl(i),4*PD+j);
            end
        end
        
        x = zeros(NPE,PD);
        for i = 1:NPE
            for j = 1:PD
                x(i,j) = NL(nl(i),j);
            end
        end
        
        coor = x';
        
        u = zeros(PD,NPE);
        for i = 1:NPE
            for j = 1:PD
                u(j,i) = ENL(nl(i),4*PD+j);
            end
        end
            
        for gp = 1:GPE
            
            epsilon = zeros(PD,PD);
            
            for i = 1:NPE
                
                J = zeros(PD,PD);
                grad = zeros(PD,NPE);
                
                [xi, eta, alpha] = GaussPoint(NPE, GPE, gp);
                
                grad_nat = grad_N_nat(NPE, xi, eta);
                
                J = coor * grad_nat';
                
                grad = inv(J)' * grad_nat;
                
                epsilon = epsilon + 1/2 *( dyad(grad(:,i),u(:,i)) + dyad(u(:,i), grad(:,i)) );
                
            end
            
            sigma = zeros(PD,PD);
            
            for a = 1:PD
                for b = 1:PD
                    for c = 1:PD
                        for d = 1:PD
                            sigma(a,b) = sigma(a,b) + constitutive(a,b,c,d, E, v) * epsilon(c,d);
                        end
                    end
                end
            end
            
            for a = 1:PD
                for b = 1:PD
                    strain(e,gp,a,b) = epsilon(a,b);
                    stress(e,gp,a,b) = sigma(a,b);
                end
            end
            
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = grad_N_nat(NPE, xi, eta)

    PD = 2;
    
    result = zeros(PD,NPE);
    
    if (NPE == 3)
        
        result(1,1) = 1;
        result(1,2) = 0;
        result(1,3) = -1;
        
        result(2,1) = 0;
        result(2,2) = 1;
        result(2,3) = -1;
        
    elseif NPE == 4
        
        result(1,1) = -1/4*(1-eta);
        result(1,2) =  1/4*(1-eta);
        result(1,3) =  1/4*(1+eta);
        result(1,4) = -1/4*(1+eta);
        
        result(2,1) = -1/4*(1-xi);
        result(2,2) = -1/4*(1+xi);
        result(2,3) =  1/4*(1+xi);
        result(2,4) =  1/4*(1-xi);

    elseif NPE == 6
        
        result(1,1) = -1 + xi;
        result(1,2) =  0;
        result(1,3) = -3 + (4*xi) + (4*eta);
        result(1,4) =  4*eta;
        result(1,5) = -4*eta;
        result(1,6) = -4 - (4*eta) + (8*xi);
        
        result(2,1) =  0;
        result(2,2) = -1 + (4*eta);
        result(2,3) = -3 + (4*xi) + (4*eta);
        result(2,4) =  4*xi;
        result(2,5) =  4 - (8*eta) - (4*xi);
        result(2,6) = -4*xi;

    elseif NPE == 8
        
        result(1,1) =  1/4*(1-eta)*(2*xi+eta);
        result(1,2) =  1/4*(1-eta)*(2*xi-eta);
        result(1,3) =  1/4*(1+eta)*(2*xi+eta);
        result(1,4) =  1/4*(1+eta)*(2*xi-eta);
        result(1,5) = -xi*(1-eta);
        result(1,6) =  1/2 *(1-eta)*(1+eta);
        result(1,7) = -xi*(1-eta);
        result(1,8) = -1/2*(1-eta)*(1+eta);
        
        result(2,1) =  1/4*(1-xi)*(xi+2*eta);
        result(2,2) =  1/4*(1+xi)*(-xi+2*eta);
        result(2,3) =  1/4*(1+xi)*(xi+2*eta);
        result(2,4) =  1/4*(1-xi)*(-xi+2*eta);
        result(2,5) = -1/2*(1-xi)*(1+xi);
        result(2,6) = -(1+xi)*eta;
        result(2,7) =  1/2*(1-xi)*(1+xi);
        result(2,8) = -(1-xi)*eta;

    elseif NPE == 9
        
        result(1,1) =  1/4*(1-eta)*(-2*xi+1)*eta;
        result(1,2) = -1/4*(1-eta)*(2*xi+1)*eta;
        result(1,3) =  1/4*(1+eta)*(2*xi+1)*eta;
        result(1,4) = -1/4*(1+eta)*(-2*xi+1)*eta;
        result(1,5) =  xi*eta*(1-eta);
        result(1,6) =  1/2*(1+2*xi)*(1+eta)*(1-eta);
        result(1,7) =  -xi*eta*(1+eta);
        result(1,8) = -1/2*(1-eta)*(1+eta)*(1-2*xi);
        result(1,9) = -2*xi*(1+eta)*(1-eta);
        
        result(2,1) =  1/4*(1-xi)*xi*(1-2*eta);
        result(2,2) = -1/4*(1+xi)*xi*(1-2*eta);
        result(2,3) =  1/4*(1+xi)*xi*(1+2*eta);
        result(2,4) = -1/4*(1-xi)*xi*(1+2*eta);
        result(2,5) =  1/2*(1-xi)*(1+xi)*(2*eta-1);
        result(2,6) = -2*(1+xi)*xi*eta;
        result(2,7) =  1/2*(1-xi)*(1+xi)*(1+2*eta);
        result(2,8) =  (1-xi)*xi*eta;
        result(2,9) = -2*eta*(1+xi)*(1-xi);
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xi, eta, alpha] = GaussPoint(NPE, GPE, gp)

    if (NPE == 4)||(NPE == 8)||(NPE == 9)
        if (GPE == 1)
            
            xi = 0; eta = 0; alpha = 4;
            
        elseif (GPE == 4)
            
            if (gp == 1)
                xi = -1/sqrt(3); eta = -1/sqrt(3); alpha = 1;
            elseif(gp == 2)
                xi = 1/sqrt(3); eta = -1/sqrt(3); alpha = 1;
            elseif(gp == 3)
                xi = 1/sqrt(3); eta = 1/sqrt(3); alpha = 1;
            elseif(gp == 4)
                xi = -1/sqrt(3); eta = 1/sqrt(3); alpha = 1;
            end
        
        elseif (GPE == 9)
            
            if (gp == 1)
                xi = -sqrt(3/5); eta = -sqrt(3/5); alpha = 25/81;
            elseif(gp == 2)
                xi = sqrt(3/5); eta = -sqrt(3/5); alpha = 25/81;
            elseif(gp == 3)
                xi = sqrt(3/5); eta = sqrt(3/5); alpha = 25/81;
            elseif(gp == 4)
                xi = -sqrt(3/5); eta = sqrt(3/5); alpha = 25/81;
            elseif(gp == 5)
                xi = 0; eta = -sqrt(3/5); alpha = 40/81;
            elseif(gp == 6)
                xi = sqrt(3/5); eta = 0; alpha = 40/81;
            elseif(gp == 7)
                xi = 0; eta = sqrt(3/5); alpha = 40/81;
            elseif(gp == 8)
                xi = -sqrt(3/5); eta= 0; alpha = 40/81;
            elseif(gp == 9)
                xi = 0; eta = 0; alpha = 64/81;
            end
            
        end
    end
    
    if (NPE == 3)||(NPE == 6)
        
        if (GPE == 1)
            
            xi = 1/3; eta = 1/3; alpha = 1;
            
        elseif (GPE == 3)
            
            if (gp == 1)
                xi = 1/6; eta = 1/6; alpha = 1/3;
            elseif(gp == 2)
                xi = 4/6; eta = 1/6; alpha = 1/3;
            elseif(gp == 3)
                xi = 1/6; eta = 4/6; alpha = 1/3;
            end
            
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = constitutive(i,j,k,l,E,nu)
   
    C = (E/(2*(1+nu))) * (delta(i,l)*delta(j,k) + delta(i,k)*delta(j,l)) + ((E*nu)/(1-nu^2)) * (delta(i,j)*delta(k,l));
              
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function delta = delta(i,j)

    if (i == j)
        delta = 1;
    else
        delta = 0;
    end
    
end
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = dyad(a,b)

    PD = 2;
    
    A = zeros(PD,PD);
    
    for i = 1:PD
        for j = 1:PD
            A(i,j) = a(i)*b(j);
        end
    end
    
end
