row_median_subtraction = np.zeros_like(image_stack)
for image_num in range(len(image_stack)):
    image_temp = image_stack[image_num]
    for row_num in range(len(image_temp)):
        row_median_subtraction[image_num, row_num] = image_temp[row_num] - \
        np.median(image_temp[row_num])
        1/0