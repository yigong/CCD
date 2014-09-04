function EventTemp = FindContainedElectrons1(Event)
%Takes a loaded 'zvent{1,}' and saves all contained electron events  

EventTest = {}; 
 n = 1; 
 for i = length(Event); 
     %EventTemp = {};
     %EventTemp= Event{out.eventNumImg(i)}; 
     %EventTest{i} = EventTemp;
     

       if ((Event{1,i}.Edep >= 200))
            EventTemp = {};
            try 
                if Event{1,i}.out.E > 10
                    EventTemp= Event{1,i}; 
                    EventTest{n} = EventTemp;
                    n = n+1;
                end
            catch e
            end
            
       end
  end
