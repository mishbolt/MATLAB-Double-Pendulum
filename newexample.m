%syms a b c d e f

S=[0,0;
    -1,-3;
    2,-2;
    2,0;
    1,3;
    5,1];
SV=['a';'b';'c';'d';'e';'f'];
r = 1.6;


for i=1:length(S)
    grid on
    labels=[['(',num2str(S(i,1)),','],[num2str(S(i,2)),')'],[', ',SV(i)]];
    %labels=[num2str(S(i,1)),num2str(S(i,2))];
    plot(S(i,1),S(i,2),'d')
    text(S(i,1),S(i,2),labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
    hold on
end

%Axes Location
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

C1=[];

for i=1:length(S)
    for j= i+1:length(S)
        if norm(S(i,:)-S(j,:))<2*r
            %C2(i,j)= 1;
            C1 =[C1;S(i,:),S(j,:)];
            plot([S(i,1),S(j,1)],[S(i,2),S(j,2)],'-r')
            
            hold on
               
        end
    end
end

del_1=zeros(length(S),length(C1));
for i=1:length(S)
    for j=1:length(C1)
        if all(S(i,:)== C1(j,1:2)) || all([S(i,:)==C1(j,3:4)])
            del_1(i,j)=1;
        end
    end
end
C2=[];
for i=1:length(S)
    for j= i+1:length(S)
        for k=j+1:length(S)
            
            
            if norm(S(i,:)-S(j,:))<2*r && norm(S(i,:)-S(k,:))<2*r && norm(S(j,:)-S(k,:))<2*r
                C2 =[C2;S(i,:),S(j,:),S(k,:)];
                pgon = polyshape([S(i,1),S(j,1),S(k,1)],[S(i,2),S(j,2),S(k,2)]);
                plot(pgon)
                hold on
            end
        end
    end
end
del_2=zeros(length(C1),length(C2));
for i=1:length(C1)
    for j=1:length(C2(:,1))
        if all(C1(i,:)== C2(j,1:4)) || all(C1(i,:)==C2(j,3:6)) || all(C1(i,:)==C2(j,[1,2,5,6]))
            del_2(i,j)=1;
        end
    end
end

