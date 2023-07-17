import sys

sums = []
r = [0] * 4
with open("res.txt", "w+") as write_file:
    for i in range(len(sys.argv) - 1):
        curr_file = sys.argv[i + 1]
        with open(curr_file, "r") as open_file:
            open_file.readline()  # ignore
            open_file.readline()  # ignore
            open_file.readline()  # ignore
            open_file.readline()  # ignore
            open_file.readline()  # ignore
            for j in range(4):
                line = open_file.readline().rstrip().split(' ')
                # print(line)
                line = list(filter(lambda st: st != '', line))
                r[j] = int(line[0].replace(',', ''))
            sums.append(r[0] + r[1] * 2 + r[2] * 4 + r[3] * 8)
            write_file.write(f'{r[0]}, {r[1]}, {r[2]}, {r[3]}\n')
    write_file.write('\n')
    write_file.write(','.join([str(i) for i in sums]))
