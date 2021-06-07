

function [s0,s1,s2]=compute_intervention_vaccine(t,Tlock,s00,s11,s22,region_name)
    
    if(strcmp(region_name,'Umbria_second_third_wave'))
       s0 = (t<Tlock(1)) + s00(1) * (t>=Tlock(1)) * (t<Tlock(2)) + s00(2)*(t>=Tlock(2))*(t<Tlock(3)) + s00(3)*(t>=Tlock(3))*(t<Tlock(4))+ s00(4)*(t>=Tlock(4))*(t<Tlock(5))+s00(5)*(t>=Tlock(5))*(t<Tlock(6))+s00(6)*(t>=Tlock(6))*(t<Tlock(7))+s00(7)*(t>=Tlock(7)); %*(t<Tlock(8))+s00(8)*(t>=Tlock(8))*(t<Tlock(9)) + s00(9)*(t>=Tlock(9))*(t<Tlock(10)) + s00(10)*(t>=Tlock(10));
       s1 = (t<Tlock(1)) + s11(1) * (t>=Tlock(1));
       s2 = (t<Tlock(1)) + s22 * (t>=Tlock(1));        
    end  
    
    if(strcmp(region_name,'Italy_second_third_wave'))
       s0 = (t<Tlock(1)) + s00(1) * (t>=Tlock(1)) * (t<Tlock(2)) + s00(2)*(t>=Tlock(2))*(t<Tlock(3)) + s00(3)*(t>=Tlock(3))*(t<Tlock(4))+ s00(4)*(t>=Tlock(4))*(t<Tlock(5))+s00(5)*(t>=Tlock(5)); %*(t<Tlock(6))+s00(6)*(t>=Tlock(6))*(t<Tlock(7))+s00(7)*(t>=Tlock(7))*(t<Tlock(8))+s00(8)*(t>=Tlock(8)); 
       s1 = (t<Tlock(1)) + s11(1) * (t>=Tlock(1));
       s2 = (t<Tlock(1)) + s22(1) * (t>=Tlock(1));        
    end 
    
end
