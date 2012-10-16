import sys

flag = False

g = open('recos.txt', 'w')

for line in open(sys.argv[1]):
	if line.strip() == 'EndReco': flag = False
	if flag: g.write(line)
	if line.strip() == 'BeginReco': flag = True
g.close()

#----------------------------

f = open('recos.txt')
lines = [l[1:-1] for l in f.readlines()]
f.close()

curr_event = ''
events = []
for line in lines:
	if len(line) > 10:
		curr_event += line
		# print curr_event

		if line[-1] in ['T', 'F']:
			events.append(curr_event)
			curr_event = ''

# for event in events:
# 	print event

# print '-----------------------'

energies = []
for event in events:
	entries = [x for x in event.split(' ') if x != '']
	# print entries
	if entries[2] == '4':
		energies.append(float(entries[6]))

# print '-----------------------'

g = open('energies.txt', 'w')
for e in energies:
	g.write(str(e))
	g.write('\n')

print max(energies)
