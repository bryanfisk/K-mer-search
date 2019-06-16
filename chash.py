import random
import string
random.seed(42)	

def numsplit(string, n):
	if len(string) < n: return string
	numsplits = len(string) // n
	return [string[x : x + n] for x in range(0, numsplits * n, n)]
		
def pad(string, n):
	if len(string) < n:
		numzeroes = n - len(string)
		string = string + [0 for _ in range(numzeroes)]
	elif len(string) > n:
		totallength = (len(string) // n) * n + n
		numzeroes = totallength - len(string)
		string = string + [0 for _ in range(numzeroes)]
	return string

def chash(string):
	A0 = int(random.random() * 65536)
	B0 = int(random.random() * 65536)
	C0 = int(random.random() * 65536)
	D0 = int(random.random() * 65536)
	K = [int(random.random() * 65536) for _ in range(64)]
	s = [int(random.random() * 32) for _ in range(64)]
	string = list(map(ord, string))
	string = pad(string, 16)
	random.shuffle(string)
	string = numsplit(string, 16)
	for i in range(len(string)):
		for j in range(16):
			A = A0
			B = B0
			C = C0
			D = D0
			if j >= 0 and j < 4:
				F = (B & C) | (~B & D)
			if j >= 4 and j < 8:
				F = (D & B) | (~D & C)
			if j >= 8 and j < 12:
				F = B ^ C ^ D
			if j >= 12 and j < 16:
				F = C ^ (B | ~D)
			F = F + A + K[j * i] + string[i][j]
			A = D
			D = C
			C = B
			B = B + F << s[j * i]
		A0 = A0 + A
		B0 = B0 + B
		C0 = C0 + C
		D0 = D0 + D

	return str(abs(A)) + str(abs(B)) + str(abs(C)) + str(abs(D)) 

collisions = {}
random.seed(None)
for _ in range(10000000):
	temp = ''.join([random.choice(string.ascii_letters) for _ in range(15)])
	temp = int(chash(temp)) % 8999999999999999 + 1000000000000000
	if temp not in collisions:
		collisions[temp] = 1
	else:
		collisions[temp] += 1

print(len([k for k, v in collisions.items() if v > 1]))