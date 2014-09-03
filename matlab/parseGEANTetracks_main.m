for idx = 10:5:40
    folder_in = '/global/home/users/ygzhang/scratch_link/g4_output/ePar_out/';
    tmp_in = ['ePar_', str(idx), '.out'];
    file_in = [folder_in, tmp_in];
    
    folder_out = folder_in;
    tmp_out = ['ePar_', str(idx), '.mat'];
    file_out = folder_out + tmp_out;
    parseGEANTetracks_func(file_in, file_out);
    


end