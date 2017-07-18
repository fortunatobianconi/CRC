%  intersection is a recursive function. It takes in input the array of the selected distance function values returned by
%  the upperLowerSet_Nr function and intersects each single set with the others. In this way, it finds
%  all the samples of the parameter vector that give simultaneously the desired behavior for all the measured output variables.

function region=intersection(distance)
 
    % final results
    region={};    
    region_name={};
    
    % base case
       if(length(distance)==4)
           region{1,1}=intersect(distance{1,1},distance{1,3});
           region_name{1,1}={'UU'};

           region{1,2}=intersect(distance{1,1},distance{1,4});
           region_name{1,2}={'UL'};

           region{1,3}=intersect(distance{1,2},distance{1,3});
           region_name{1,3}={'LU'};

           region{1,4}=intersect(distance{1,2},distance{1,4});
           region_name{1,4}={'LL'};
     
     % recursive case       
       else
           partial_region=intersection(distance(1:length(distance)-2));

           for i=1:size(partial_region,2)
             region=[region intersect(partial_region{1,i},distance{1,length(distance)-1})];  %Upper
             if ~isempty(intersect(partial_region{1,i},distance{1,length(distance)-1}))
                 region_name{1,end+1}=strcat(partial_region{2,i},'U');
             end
             region=[region intersect(partial_region{1,i},distance{length(distance)})];   %Lower
             if ~isempty(intersect(partial_region{1,i},distance{1,length(distance)}))
                 region_name{1,end+1}=strcat(partial_region{2,i},'L');
             end
           end
       end
       
       region=[region;region_name];
end       