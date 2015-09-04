 EventTest = {}; 
 error = [];
 for i = length(Event); 
     %EventTemp = {};
     %EventTemp= Event{out.eventNumImg(i)}; 
     %EventTest{i} = EventTemp;
     

       if ((Event{1,i}.Edep >= 1000))
            EventTemp = {};
            try 
                if Event{1,i}.out.E > 1000
                    EventTemp= Event{out.eventNumImg(i)}; 
                    EventTest{i} = EventTemp;
                    
                end
            catch e
                error = [ error e];
            end
                
       end
  end
