import sys


def write_file(out_file, data, separator, header):
    line = ""
    with open(out_file, 'w') as outfile:
        #outfile.write(str('#') + str(header) + "\n")
        for row in range(len(data)):
            for col in range(len(data[row])):
                if col != len(data[row]) - 1:
                    line = line + str(data[row][col]) + str(separator)
                elif col == len(data[row]) - 1:
                    line = line + str(data[row][col])
            outfile.write( line + "\n" )
            line = ""

def read_file(in_file, separator):
    with open(in_file, 'r') as infile:
        if separator == 'single_line':
            g_list = []
            row = []
            while(True):
                line = infile.readline()
                if not line: break
                if line[0] != '#'and line[0] != '\\' and line[0] != '|':
                    g_list.append(float(line))
        else:
            g_list = []
            row = []
            while(True):
                line = infile.readline()
                if not line: break
                if line[0] != '#' and line[0] != '\\' and line[0] != '|':
                    if separator == ' ':
                        string = line.split()
                    else:
                        string = line.split(separator)
                    for col in range(0,len(string)):
                        row.append(string[col])
                    g_list.append(row)
                    row = []
    return g_list

def filter_interpolation(lw):
    """
    """
    if lw > 370.0 and lw < 410.0:
        return (lw - 370.0)*(0.365 - 0.35)/(410.0 - 370.0) + 0.35
    elif lw >= 410.0 and lw < 455.0:
        return (lw - 410.0)*(0.544 - 0.365)/(455.0 - 410.0) + 0.365
    elif lw >= 455.0 and lw < 470.0:
        return (lw - 455.0)*(0.68 - 0.544)/(470.0 - 455.0) + 0.544
    elif lw >= 470.0 and lw < 600.0:
        return (lw - 470.0)*(0.68 - 0.68)/(600.0 - 470.0) + 0.68
    elif lw >= 600.0 and lw < 650.0:
        return (lw - 600.0)*(0.49 - 0.68)/(650.0 - 600.0) + 0.68
    elif lw >= 650.0 and lw < 700.0:
        return (lw - 650.0)*(0.31 - 0.49)/(700.0 - 650.0) + 0.49
    elif lw >= 700.0 and lw < 800.0:
        return (lw - 700.0)*(0.0 - 0.31)/(800.0 - 700.0) + 0.31
    else:
        return 0.0

if __name__ == '__main__':
    """
    """
    filter_data = read_file(sys.argv[1], ' ')

    phot_filter = list()
    for i in range(len(filter_data)):
        phot_filter.append([float(filter_data[i][0]), filter_interpolation(1000.0*float(filter_data[i][0]))])

    write_file('most_hires.dat', phot_filter, ' ', '')

    print "Visca!"