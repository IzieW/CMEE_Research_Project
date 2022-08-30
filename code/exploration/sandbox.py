

count120, count240, count360 = 0,0,0
for i in range(240):
    if i < 120:
        i =i
    else:
        i -= 120

    if i < 40:
        count120 +=1
    elif i < 80:
        count240 +=1
    else:
        count360 +=1
