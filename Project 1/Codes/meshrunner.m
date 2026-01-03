clc;
clear all;
clf;

d1 = 1;
d2 = 1;
p = 10;
m = 4;
R = 0;

% [NL, EL] = void_mesh_rhombus(d1, d2, p, m, R);
% [NL, EL] = void_mesh_square(d1, d2, p, m, R);
% [NL, EL] = void_mesh(d1, d2, p, m, R);

NoN = size(NL,1);
NoE = size(EL,1);



%         for i = 1:NoE
%             hold on;
%             plot([NL(EL(i,1),1), NL(EL(i,2),1)],'m');
            
           
%         end

        %nodeları plot ediyor
        
        for i = 1:NoN
            hold on;
            plot(NL(i,1),NL(i,2),'o','MarkerSize',15, 'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1])
            text(NL(i,1),NL(i,2), num2str(i), 'Color','w','FontSize',12,'HorizontalAlignment','center')
        end
        
        axis equal      

        
        

        
      
       
        
        
        


