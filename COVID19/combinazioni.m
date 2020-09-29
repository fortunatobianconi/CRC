%la funzione calcola l'intersezione fra gli insiemi di indici delle Cloud
%in cui le proteine hanno valore alto o basso
%L'ingresso � un insieme di celle di tipo C={U1 L1 U2 L2 ...} dove U1 �
%l'insieme degli indici in cui la prima proteina ha valori alti dell'area

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
%            for i=1:length(risultatiProvvisori)
%                risultati{1,i}=intersect(risultatiProvvisori{1,i},C{1,length(C)-1});
%                risultati{2,i}=intersect(risultatiProvvisori{1,i},C{length(C)});
%            end
           for i=1:size(risultatiProvvisori,2)
            risultati=[risultati intersect(risultatiProvvisori{1,i},C{1,length(C)-1})];  %Upper
             if ~isempty(intersect(risultatiProvvisori{1,i},C{1,length(C)-1}))
                 ris{1,end+1}=strcat(risultatiProvvisori{2,i},'T');
             end
             risultati=[risultati intersect(risultatiProvvisori{1,i},C{length(C)})];   %Lower
             if ~isempty(intersect(risultatiProvvisori{1,i},C{1,length(C)}))
                 ris{1,end+1}=strcat(risultatiProvvisori{2,i},'F');
             end
%            A=intersect(risultatiProvvisori{1,i},C{1,length(C)-1});
%            B=intersect(risultatiProvvisori{1,i},C{1,length(C)});
%            risultati

           end
       end
       risultati2=[risultati;ris];
end       