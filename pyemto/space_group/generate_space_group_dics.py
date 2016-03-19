file = open('space_groups.txt','r')
file_lines = file.readlines()

number_symbol = []

for i in range(len(file_lines)):
    line_split = file_lines[i].split()
    for j in range(len(line_split)):
        number_symbol.append(line_split[j])

number = number_symbol[::2]
symbol = number_symbol[1::2]


ibz_list = []
bl_list = []

for i in range(1,len(number)+1):
    # Trigonal
    if i >= 1 and i <= 2:
        #print(symbol[i][0])
        if symbol[i-1][0] == 'P':
            ibz_list.append(14)
            bl_list.append('simple triclinic')
        else:
            print('spage group {0} is missing!'.format(i))
    # Monoclinic
    elif i >= 3 and i <= 15:
        if symbol[i-1][0] == 'P':
            ibz_list.append(12)
            bl_list.append('simple monoclinic')
        elif symbol[i-1][0] == 'C':
            ibz_list.append(13)
            bl_list.append('base-centered monoclinic')
        else:
            print('spage group {0} is missing!'.format(i))
    # Orthorhombic
    elif i >= 16 and i <= 74:
        if symbol[i-1][0] == 'P':
            ibz_list.append(8)
            bl_list.append('simple orthorhombic')
        elif symbol[i-1][0] == 'I':
            ibz_list.append(10)
            bl_list.append('body-centered orthorhombic')
        elif symbol[i-1][0] == 'A' or symbol[i-1][0] == 'B' or symbol[i-1][0] == 'C':
            ibz_list.append(9)
            bl_list.append('base-centered orthorhombic')
        elif symbol[i-1][0] == 'F':
            ibz_list.append(11)
            bl_list.append('face-centered orthorhombic')
        else:
            print('spage group {0} is missing!'.format(i))
    # Tetragonal
    elif i >= 75 and i <= 142:
        if symbol[i-1][0] == 'P':
            ibz_list.append(5)
            bl_list.append('simple tetragonal')
        elif symbol[i-1][0] == 'I':
            ibz_list.append(6)
            bl_list.append('body-centered tetragonal')
        else:
            print('spage group {0} is missing!'.format(i))
    # Trigonal/Rhombohedral
    elif i >= 143 and i <= 167:
        if symbol[i-1][0] == 'P':
            ibz_list.append(4)
            bl_list.append('hexagonal')
        elif symbol[i-1][0] == 'R':
            ibz_list.append(7)
            bl_list.append('rhombohedral')
        else:
            print('spage group {0} is missing!'.format(i))
    # Hexagonal
    elif i >= 168 and i <= 194:
        if symbol[i-1][0] == 'P':
            ibz_list.append(4)
            bl_list.append('hexagonal')
        else:
            print('spage group {0} is missing!'.format(i))
    # Cubic
    elif i >= 195 and i <= 230:
        if symbol[i-1][0] == 'P':
            ibz_list.append(1)
            bl_list.append('simple cubic')
        elif symbol[i-1][0] == 'I':
            ibz_list.append(3)
            bl_list.append('body-centered cubic')
        elif symbol[i-1][0] == 'F':
            ibz_list.append(2)
            bl_list.append('face-centered cubic')
        else:
            print('spage group {0} is missing!'.format(i))

# Generate dictionaries using some reasonable formatting.
line =  '        self.sg2ibz = {'
for i in range(23):
    if i == 0:
        line += '{0}:{1}, {2}:{3}, {4}:{5}, {6}:{7}, {8}:{9}, {10}:{11}, {12}:{13}, {14}:{15}, {16}:{17}, {18}:{19},\n'.format(
        number[10*i+0],ibz_list[10*i+0],number[10*i+1],ibz_list[10*i+1],number[10*i+2],ibz_list[10*i+2],number[10*i+3],ibz_list[10*i+3],number[10*i+4],ibz_list[10*i+4],
        number[10*i+5],ibz_list[10*i+5],number[10*i+6],ibz_list[10*i+6],number[10*i+7],ibz_list[10*i+7],number[10*i+8],ibz_list[10*i+8],number[10*i+9],ibz_list[10*i+9])
    else:
        line += '                       {0}:{1}, {2}:{3}, {4}:{5}, {6}:{7}, {8}:{9}, {10}:{11}, {12}:{13}, {14}:{15}, {16}:{17}, {18}:{19},\n'.format(
        number[10*i+0],ibz_list[10*i+0],number[10*i+1],ibz_list[10*i+1],number[10*i+2],ibz_list[10*i+2],number[10*i+3],ibz_list[10*i+3],number[10*i+4],ibz_list[10*i+4],
        number[10*i+5],ibz_list[10*i+5],number[10*i+6],ibz_list[10*i+6],number[10*i+7],ibz_list[10*i+7],number[10*i+8],ibz_list[10*i+8],number[10*i+9],ibz_list[10*i+9])

print(line)
print('\n')

line =  '        self.sg2bl = {'
for i in range(46):
    if i == 0:
        line += '{0}:\'{1}\', {2}:\'{3}\', {4}:\'{5}\', {6}:\'{7}\', {8}:\'{9}\',\n'.format(
        number[5*i+0],bl_list[5*i+0],number[5*i+1],bl_list[5*i+1],number[5*i+2],bl_list[5*i+2],number[5*i+3],bl_list[5*i+3],number[5*i+4],bl_list[5*i+4])
    else:
        line += '                      {0}:\'{1}\', {2}:\'{3}\', {4}:\'{5}\', {6}:\'{7}\', {8}:\'{9}\',\n'.format(
        number[5*i+0],bl_list[5*i+0],number[5*i+1],bl_list[5*i+1],number[5*i+2],bl_list[5*i+2],number[5*i+3],bl_list[5*i+3],number[5*i+4],bl_list[5*i+4])

print(line)
