

function [s0,s1,s2]=compute_intervention(t,Tlock,s00,s11,s22,region_name)


 if(strcmp(region_name,'Italy'))
       s0 = (t<Tlock(1)) + s00(1) * (t>=Tlock(1)) * (t<Tlock(2)) + s00(2)*(t>=Tlock(2)) * (t<Tlock(3)) + s00(3)*(t>=Tlock(3))*(t<Tlock(4)) + s00(4)*(t>=Tlock(4));%*(t<Tlock(5)) + s00(5)*(t>=Tlock(5))*(t<Tlock(6)) + s00(6)*(t>=Tlock(6));
       s1 = (t<Tlock(4)) + s11(1) * (t>=Tlock(4));
       s2 = (t<Tlock(1)) + s22 * (t>=Tlock(1));        
    end
    
    
    
    if(strcmp(region_name,'Umbria'))
       s0 = (t<Tlock(2)) + s00(1) * (t>=Tlock(2)) * (t<Tlock(3)) + s00(2)*(t>=Tlock(3));
       s1 = (t<Tlock(3)) + s11(1) * (t>=Tlock(3));
       s2 = (t<Tlock(1)) + s22 * (t>=Tlock(1));        
    end
    
    
end
