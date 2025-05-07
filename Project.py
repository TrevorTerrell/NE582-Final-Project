import numpy as np
import matplotlib.pyplot as plt
ran = np.random.default_rng()

class crossSection():
    def __init__(self, num_mat, file):
        self.s = np.zeros((num_mat, 3))
        self.c = np.zeros((num_mat, 3))
        self.f = np.zeros((num_mat, 3))
        self.nu = np.zeros((num_mat))
        for i in np.arange(num_mat):
            line = file.readline() #number of material

            line = file.readline() #scattering
            line = line.replace("    Scatter: ", "")
            temp = line.split(",")
            for j in np.arange(3):
                self.s[i, j] = float(temp[j])
            
            line = file.readline() #capture
            line = line.replace("    Capture: ", "")
            temp = line.split(",")
            for j in np.arange(3):
                self.c[i, j] = float(temp[j])
            
            line = file.readline() #fission
            line = line.replace("    Fission: ", "")
            temp = line.split(",")
            for j in np.arange(3):
                self.f[i, j] = float(temp[j])

            line = file.readline() #nu
            line = line.replace("    Nu: ", "")
            self.nu[i] = float(line)

            line = file.readline()

        assert np.sum(self.s, 0)[2] == 0, "Scattering cross sections in group 3 must be 0"

        self.tot = np.sum(self.s + self.c + self.f, 0)

    def get_nu(self, material):
        nu = self.nu[material]
        val = np.floor(nu)
        chance = nu - val
        if ran.random() < chance:
            val += 1
        return val
            
class neutron():
    def __init__(self, x, y, xsec: type[crossSection]):
        self.e = 0
        self.x = x
        self.y = y
        self.r = -np.log(1 - ran.random())/xsec.tot[self.e]
        self.mu = (1 - ran.random()) * (ran.integers(0,2) * 2 - 1)
        self.eta = (1 - ran.random()) * (ran.integers(0,2) * 2 - 1)
        norm_factor = np.sqrt(np.pow(self.mu, 2) + np.pow(self.eta, 2))
        self.mu /= norm_factor
        self.eta /= norm_factor
        del norm_factor

    def new_r(self, xsec: type[crossSection]):
        self.r = -np.log(1 - ran.random())/xsec.tot[self.e]
    
    def new_dir(self):
        self.mu = (1 - ran.random()) * (ran.integers(0,2) * 2 - 1)
        self.eta = (1 - ran.random()) * (ran.integers(0,2) * 2 - 1)
        norm_factor = np.sqrt(np.pow(self.mu, 2) + np.pow(self.eta, 2))
        self.mu /= norm_factor
        self.eta /= norm_factor
        del norm_factor

    def reflect_x(self):#reflecting on x = right boundary
        self.mu = - self.mu

    def reflect_y(self):
        self.eta = - self.eta

    def collType(self, xsec):
        val = ran.random() * xsec.tot[self.e]
        j = 0
        if val < np.sum(xsec.s[:, self.e]):
            for i in np.arange(len(xsec.s[:, self.e])):
                if val > np.sum(xsec.s[:(i + 1), self.e]):
                    j += 1
            return 0, j
        elif val < np.sum(xsec.s[:, self.e]) + np.sum(xsec.c[:, self.e]):
            val -= np.sum(xsec.s[:, self.e])
            for i in np.arange(len(xsec.c[:, self.e])):
                if val > np.sum(xsec.c[:(i + 1), self.e]):
                    j += 1
            return 1, j
        elif val < np.sum(xsec.s[:, self.e]) + np.sum(xsec.c[:, self.e]) + np.sum(xsec.f[:, self.e]):
            val -= np.sum(xsec.s[:, self.e]) + np.sum(xsec.c[:, self.e])
            for i in np.arange(len(xsec.f[:, self.e])):
                if val > np.sum(xsec.f[:(i + 1), self.e]):
                    j += 1
            return 2, j

class tally():
    def __init__(self, bins, xmax, ymax):
        self.flux = np.zeros((3, bins, bins))#energy group, xbin, y bin
        self.fis = np.zeros((3, bins, bins))
        x_ref = np.ones(bins) * xmax / bins
        x_ref = np.cumsum(x_ref)
        self.X_REF = x_ref
        y_ref = np.ones(bins) * ymax / bins
        y_ref = np.cumsum(y_ref)
        self.Y_REF = y_ref

    def normalize(self, num_hist):
        self.flux /= num_hist * self.X_REF[0] * self.Y_REF[0]
        self.fis /= num_hist * self.X_REF[0] * self.Y_REF[0]

def MCTMa(n: type[neutron], gen_n: type[list], xsec: type[crossSection], tallies: type[tally]):
    n.x += n.r * n.mu
    n.y += n.r * n.eta

    if n.x > tallies.X_REF[-1]:
        #neutron hits the reflection boundary
        n.x = 2 * tallies.X_REF[-1] - n.x
        n.reflect_x()
    if n.y > tallies.Y_REF[-1]:
        #neutron hits the reflection boundary
        n.y = 2 * tallies.Y_REF[-1] - n.y
        n.reflect_y()

    #if leaked into vacuum, neutron dies
    if (n.x < 0) or (n.y < 0):
        return gen_n

    #find bin number
    i = 0
    j = 0
    while n.x > tallies.X_REF[i]:
        i += 1
    while n.y > tallies.Y_REF[j]:
        j +=1

    #collision
    c, m = n.collType(xsec)
    match c:
        case 0:
            #scattering in material m
            #update tallies
            tallies.flux[n.e, i, j] += 1/xsec.s[m, n.e]
            #resample distance and direction
            n.new_r(xsec)
            n.new_dir()
            n.e += 1
            #rerun
            return MCTMa(n, gen_n, xsec, tallies)
        case 1:
            #capture in material m
            #update tallies
            tallies.flux[n.e, i, j] += 1/xsec.c[m, n.e]
            return gen_n
        case 2:
            #fission in material m
            #update tallies
            tallies.flux[n.e, i, j] += 1/xsec.f[m, n.e]
            #save neutron birth positions
            array = np.array([[n.x], [n.y]])
            for k in np.arange(xsec.get_nu(m)):
                gen_n = np.append(gen_n, array, 1)
            return gen_n
        
def MCTMi(n: type[neutron], gen_n: type[list], xsec: type[crossSection], tallies: type[tally]):
    n.x += n.r * n.mu
    n.y += n.r * n.eta

    if n.x > tallies.X_REF[-1]:
        #neutron hits the reflection boundary
        n.x = 2 * tallies.X_REF[-1] - n.x
        n.reflect_x()
    if n.y > tallies.Y_REF[-1]:
        #neutron hits the reflection boundary
        n.y = 2 * tallies.Y_REF[-1] - n.y
        n.reflect_y()

    #if leaked into vacuum, neutron dies
    if (n.x < 0) or (n.y < 0):
        return gen_n

    #find bin number
    i = 0
    j = 0
    while n.x > tallies.X_REF[i]:
        i += 1
    while n.y > tallies.Y_REF[j]:
        j +=1

    #collision
    c, m = n.collType(xsec)
    match c:
        case 0:
            #scattering in material m
            #resample distance and direction
            n.new_r(xsec)
            n.new_dir()
            n.e += 1
            #rerun
            return MCTMa(n, gen_n, xsec, tallies)
        case 1:
            #capture in material m
            return gen_n
        case 2:
            #fission in material m
            #save neutron birth positions
            array = np.array([[n.x], [n.y]])
            for k in np.arange(xsec.get_nu(m)):
                gen_n = np.append(gen_n, array, 1)
            return gen_n

def simulate(num_hist, num_gen, num_bins, source_x, source_y, x_edge, y_edge, xsec):
    xs = np.ones(num_hist) * source_x
    ys = np.ones(num_hist) * source_y
    pop_this_gen = np.array((xs,ys))
    del xs
    del ys
    pop_next_gen = np.zeros((2,0))

    tallies = tally(num_bins, x_edge, y_edge)

    #generation loop
    for g in np.arange(num_gen):
        if g == num_gen - 1:#active
            for i in np.arange(len(pop_this_gen[0,:])):
                n = neutron(pop_this_gen[0,i], pop_this_gen[1,i], xsec)
                pop_next_gen = MCTMa(n, pop_next_gen, xsec, tallies)
            tallies.normalize(len(pop_this_gen[0,:]))

        else:#Inactive
            for i in np.arange(len(pop_this_gen[0,:])):
                n = neutron(pop_this_gen[0,i], pop_this_gen[1,i], xsec)
                pop_next_gen = MCTMi(n, pop_next_gen, xsec, tallies)
            tallies.normalize(len(pop_this_gen[0,:]))

        print(f'{g}: k = {len(pop_next_gen[0,:])/len(pop_this_gen[0,:])}, N = {len(pop_this_gen[0,:])}')
        if len(pop_this_gen[0,:]) > num_hist * 2:#lower the population to prevent long simulation times
            print(f'\tCulling Population...')
            temp = np.flip(np.arange(int(len(pop_next_gen[0,:])/2)) * 2, 0)
            for i in temp:
                pop_next_gen = np.delete(pop_next_gen, i, 1)
        if len(pop_this_gen[0,:]) < num_hist / 2:#raise the population to prevent botteming out
            print(f'Cloning Population...')
            pop_next_gen = np.append(pop_next_gen, pop_next_gen, 1)
        del pop_this_gen
        pop_this_gen = np.copy(pop_next_gen)
        pop_next_gen = np.zeros((2,0))

    return tallies

print("Reading File")
#Step 1: get info from input file
file = open("ProjectInput.txt", "r")
line = file.readline() #Number of histories
line = line.replace("N: ", "")
hist = int(line)
line = file.readline() #Number of generations
line = line.replace("G: ", "")
gen = int(line)
line = file.readline() #Number of bins
line = line.replace("Bins: ", "")
bins = int(line)
line = file.readline() #Width
line = line.replace("Width: ", "")
d_x = float(line)
line = file.readline() #Height
line = line.replace("Height: ", "")
d_y = float(line)
line = file.readline() #Source x
line = line.replace("Source x: ", "")
x0 = float(line)
line = file.readline() #Source y
line = line.replace("Source y: ", "")
y0 = float(line)
line = file.readline() #Number of Materials
line = line.replace("Number of Materials: ", "")
num_mat = int(line)

xsec = crossSection(num_mat, file)

del file
del line

assert d_x > 0, "System width must be grater than 0."
assert d_y > 0, "System height must be grater than 0."
assert (x0 >= 0) & (x0 <= d_x) & (y0 >= 0) & (y0 <= d_y), "Initial source location must be in bounds."


print('Simulating...')
tallies = simulate(hist, gen, bins, x0, y0, d_x, d_y, xsec)
total_flux = np.sum(tallies.flux, 0)

plt.figure(0)
plt.imshow(total_flux, norm = 'log', vmin = np.min(total_flux), vmax = np.max(total_flux), origin = 'lower')
locs, labels = plt.xticks()
locs = np.delete(locs, 0)
labels = np.delete(labels, 0)
locs = np.delete(locs, -1)
labels = np.delete(labels, -1)
labels = [float(i)*d_x/bins for i in locs]
plt.xticks(locs, labels)
locs, labels = plt.yticks()
locs = np.delete(locs, 0)
labels = np.delete(labels, 0)
locs = np.delete(locs, -1)
labels = np.delete(labels, -1)
labels = [float(i)*d_y/bins for i in locs]
plt.yticks(locs, labels)
plt.colorbar()
plt.show()