import math
p = (
151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,36,103,
30,69,142,8,99,37,240,21,10,23,190,6,148,247,120,234,75,0,26,197,
62,94,252,219,203,117,35,11,32,57,177,33,88,237,149,56,87,174,20,
125,136,171,168,68,175,74,165,71,134,139,48,27,166,77,146,158,231,
83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,102,
143,54,65,25,63,161,1,216,80,73,209,76,132,187,208,89,18,169,200,
196,135,130,116,188,159,86,164,100,109,198,173,186,3,64,52,217,226,
250,124,123,5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,
58,17,182,189,28,42,223,183,170,213,119,248,152,2,44,154,163,70,
221,153,101,155,167,43,172,9,129,22,39,253,19,98,108,110,79,113,
224,232,178,185,112,104,218,246,97,228,251,34,242,193,238,210,144,
12,191,179,162,241,81,51,145,235,249,14,239,107,49,192,214,31,181,
199,106,157,184,84,204,176,115,121,50,45,127,4,150,254,138,236,
205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180,
151,160,137,91,90,15,131,13,201,95,96,53,194,233,7,225,140,36,103,
30,69,142,8,99,37,240,21,10,23,190,6,148,247,120,234,75,0,26,197,
62,94,252,219,203,117,35,11,32,57,177,33,88,237,149,56,87,174,20,
125,136,171,168,68,175,74,165,71,134,139,48,27,166,77,146,158,231,
83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,102,
143,54,65,25,63,161,1,216,80,73,209,76,132,187,208,89,18,169,200,
196,135,130,116,188,159,86,164,100,109,198,173,186,3,64,52,217,226,
250,124,123,5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,
58,17,182,189,28,42,223,183,170,213,119,248,152,2,44,154,163,70,
221,153,101,155,167,43,172,9,129,22,39,253,19,98,108,110,79,113,
224,232,178,185,112,104,218,246,97,228,251,34,242,193,238,210,144,
12,191,179,162,241,81,51,145,235,249,14,239,107,49,192,214,31,181,
199,106,157,184,84,204,176,115,121,50,45,127,4,150,254,138,236,
205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180)
  
def lerp(t, a, b):
    return a + t * (b - a)
  
def fade(t):
    return t * t * t * (t * (t * 6 - 15) + 10)

def grad1D(hash, x):
    h = hash & 15
    if h < 8:
        u = x
    else:
        u = 0
    if h < 4:
        v = 0
    elif h == 12 or h == 14:
        v = x
    else:
        v = 0
    if h & 1 != 0:
        u = -u
    if h & 2 != 0:
        v = -v
    return u + v

def grad2D(hash, x, y):
    h = hash & 15
    if h < 8:
        u = x
    else:
        u = y
    if h < 4:
        v = y
    elif h == 12 or h == 14:
        v = x
    else:
        v = 0
    if h & 1 != 0:
        u = -u
    if h & 2 != 0:
        v = -v
    return u + v
  
def grad3D(hash, x, y, z):
    h = hash & 15
    if h < 8:
        u = x
    else:
        u = y
    if h < 4:
        v = y
    elif h == 12 or h == 14:
        v = x
    else:
        v = z
    if h & 1 != 0:
        u = -u
    if h & 2 != 0:
        v = -v
    return u + v

def pnoise1D(x, y):
    global p
    X = int(math.floor(x)) & 255
    Y = int(math.floor(y)) & 255
    x -= math.floor(x)
    y -= math.floor(y)
    
    u = fade(x)
    v = fade(y)
    
    A =  p[X]
    B =  p[(X + 1) & 255]
    
    gradA = grad1D(A, x)
    gradB = grad1D(B, x-1)
    return lerp(u, gradA, gradB)
 
def pnoise2D(x, y):
    global p
    X = int(math.floor(x)) & 255
    Y = int(math.floor(y)) & 255
    x -= math.floor(x)
    y -= math.floor(y)
    
    u = fade(x)
    v = fade(y)
    
    A = p[X] + Y
    B = p[(X + 1) & 255] + Y
    
    AA = p[A]
    AB = p[(A + 1) & 255]
    BA = p[B]
    BB = p[(B + 1) & 255]
    
    gradAA = grad2D(AA, x,   y)
    gradBA = grad2D(BA, x-1, y)
    gradAB = grad2D(AB, x,   y-1)
    gradBB = grad2D(BB, x-1, y-1)
    return lerp(v, lerp(u, gradAA, gradBA), lerp(u, gradAB, gradBB)) 

def pnoise3D(x, y, z):
    global p
    X = int(math.floor(x)) & 255
    Y = int(math.floor(y)) & 255
    Z = int(math.floor(z)) & 255
    x -= math.floor(x)
    y -= math.floor(y)
    z -= math.floor(z)
    
    u = fade(x)
    v = fade(y)
    w = fade(z)
    
    A = p[X] + Y
    B = p[(X + 1) & 255] + Y
    
    AA = p[A] + Z
    AB = p[(A + 1) & 255] + Z
    BA = p[B] + Z
    BB = p[(B + 1) & 255] + Z
    
    AAA = p[AA]
    ABA = p[AB]
    BAA = p[BA]
    BBA = p[BB]
    AAB = p[AA + 1]
    BAB = p[BA + 1]
    ABB = p[AB + 1]
    BBB = p[BB + 1]
    
    gradAAA = grad3D(AAA, x,   y,   z)
    gradBAA = grad3D(BAA, x-1, y,   z)
    gradABA = grad3D(ABA, x,   y-1, z)
    gradBBA = grad3D(BBA, x-1, y-1, z)
    gradAAB = grad3D(AAB,x,   y,   z-1)
    gradBAB = grad3D(BAB,x-1, y,   z-1)
    gradABB = grad3D(ABB,x,   y-1, z-1)
    gradBBB = grad3D(BBB,x-1, y-1, z-1)
    return lerp(w, 
    lerp(v, lerp(u, gradAAA, gradBAA), lerp(u, gradABA, gradBBA)),
    lerp(v, lerp(u, gradAAB, gradBAB), lerp(u, gradABB, gradBBB)))

def noise1D(x, depth, roughness):
	f = 1.0
	n = 0.0
	sum_f = 0.0
	
	for i in range(0, depth + 1):
		n += f * pnoise1D(x)
		sum_f += f
		f *= roughness
		x *= 2.0
	
	return n / sum_f

def noise2D(x, y, depth, roughness):
	f = 1.0
	n = 0.0
	sum_f = 0.0
	
	for i in range(0, depth + 1):
		n += f * pnoise2D(x, y)
		sum_f += f
		f *= roughness
		x *= 2.0
		y *= 2.0
	
	return n / sum_f

def noise3D(x, y, z, depth, roughness):
	f = 1.0
	n = 0.0
	sum_f = 0.0
	
	for i in range(0, depth + 1):
		n += f * pnoise3D(x, y, z)
		sum_f += f
		f *= roughness
		x *= 2.0
		y *= 2.0
		z *= 2.0
	
	return n / sum_f
