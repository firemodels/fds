def write_data_to_csv(data_set):
    import csv
    
    filename = 'workfile.csv'
    f = open(filename, 'w')
    fileWriter = csv.writer(f,delimiter=',')
    
    list_of_lengths = [[len(column) for column in record_set] for record_set in data_set]
    max_set = max(list_of_lengths)
    
    temp_list = []
    
    for group in range(len(data_set)):
        for col in range(len(data_set[group])):
            row_list=[]
            
            for value in data_set[group][col]:
                row_list.append(str(value))
                
            if len(row_list) < max_set[0]:
                for i in range(max_set[0]-len(row_list)):
                    row_list.append('')
                temp_list.append(row_list)
                
            elif len(row_list) == max_set[0]:
                temp_list.append(row_list)
                
    final_list = zip(*temp_list)
    
    for line in final_list:
        fileWriter.writerow(line)
        
    f.close()
    print 'Data written to:', filename

list_of_lists = [[['a','b','c','d','e','f'],['g','h','i','j','k','l']],[[1,2,3,4,5,6,7,8,9,10,11,12],[1,2,3,4,5,6,7,8,9,10,11,12]],[[13.5,14.5,15.5,16.5],[17.5,18.5,19.5,20.5]],[["word","to"],["your","momma"]]]

write_data_to_csv(list_of_lists)