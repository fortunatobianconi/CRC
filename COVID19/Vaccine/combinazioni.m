%the function computes the intersection between the sets of indexes of the
%Cloud where values are high or low
%The input is a seires of cells C={U1 L1 U2 L2 ...} where U1 is the set of
%indexes where the chosen variable has high values

function risultati2=combinazioni(C)
    risultati2={};    
    risultati={};
    ris={};
    %caso base
       if(length(C)==4)
           %ris{1,1}={'UU'};
           risultati{1,1}=intersect(C{1,1},C{1,3});
           ris{1,1}={'TT'};
           %ris{1,2}={'UL'};
           risultati{1,2}=intersect(C{1,1},C{1,4});
           ris{1,2}={'TF'};
           %ris{1,3}={'LU'};
           risultati{1,3}=intersect(C{1,2},C{1,3});
           ris{1,3}={'FT'};
           %ris{1,4}={'LL'};
           risultati{1,4}=intersect(C{1,2},C{1,4});
           ris{1,4}={'FF'};
           
       else
           risultatiProvvisori=combinazioni(C(1:length(C)-2));

           for i=1:size(risultatiProvvisori,2)
            risultati=[risultati intersect(risultatiProvvisori{1,i},C{1,length(C)-1})];  %Upper
             if ~isempty(intersect(risultatiProvvisori{1,i},C{1,length(C)-1}))
                 ris{1,end+1}=strcat(risultatiProvvisori{2,i},'T');
             end
             risultati=[risultati intersect(risultatiProvvisori{1,i},C{length(C)})];   %Lower
             if ~isempty(intersect(risultatiProvvisori{1,i},C{1,length(C)}))
                 ris{1,end+1}=strcat(risultatiProvvisori{2,i},'F');
             end
           end
       end
       risultati2=[risultati;ris];
end       