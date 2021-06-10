import numpy as np
import itertools

class analysis:
    def __init__(self, groFile, traj=None):
        self.gro = groFile
        self._traj = traj
        self._box = None
        self.g = None
        self.rho = None
        self.conversion = 0.376
        self.expTraj = None
        
    @property
    def traj(self):
        if self._traj is None:
            self._traj = self.get_traj()
        return self._traj
    
    @property
    def box(self):
        if self._box is None:
            self._box = self.get_box()
        return self._box
    
    def get_traj(self):
        with open(self.gro, 'r') as f:
            lines = f.readlines()
            numAtoms = int(lines[1])
            numFrames = int(len(lines)/(numAtoms +3))
            traj = np.zeros((numFrames,numAtoms,3))
            frame = -1
            atom = 0
            for line in lines:
                if line=='\n' and frame<numFrames-1: frame+=1
                elif line=='\n': break
                elif line==lines[1]: pass 
                elif line==lines[-1]: pass
                else: 
                    traj[frame][atom] = self.parseAtomLine(line)
                    #print(f"traj[{frame}][{atom}] = {traj[frame][atom]}")
                    if atom<numAtoms-1: atom += 1
                    else: atom = 0
        return traj

    def get_box(self):
        with open(self.gro, 'r') as f:
            lines = f.readlines()
            box = np.fromstring(lines[-1], sep=' ')/self.conversion
        return box
        
    def parseAtomLine(self, line):
        pos = line.find('.')
        x = line[pos-1:pos+4]
        y = line[pos+7:pos+12]
        z = line[pos+15:pos+20]
        xyz = [float(x)/self.conversion, float(y)/self.conversion, float(z)/self.conversion]
        #print(xyz)
        return xyz
    
    def calcDistances(self, frame):
        distances = []
        #print(frame)
        for i in range(0,len(frame)-1):
            for e in range(i+1, len(frame)):
                distances.append(np.sqrt(pow(frame[i][0]-frame[e][0], 2)+pow(frame[i][1]-frame[e][1], 2)+pow(frame[i][2]-frame[e][2], 2)))
        return distances
    
    def calcDistancesPBC(self, frame, xframes):
        distances = []
        #print(frame)
        for i in range(0,len(frame)-1):
            for e in range(i+1, len(frame)):
                distances.append(np.sqrt(pow(frame[i][0]-frame[e][0], 2)+pow(frame[i][1]-frame[e][1], 2)+pow(frame[i][2]-frame[e][2], 2)))
            for e in range(len(xframes)):
                distances.append(np.sqrt(pow(frame[i][0]-xframes[e][0], 2)+pow(frame[i][1]-xframes[e][1], 2)+pow(frame[i][2]-xframes[e][2], 2)))
        return distances

    def generateHist(self, traj, rMax, step, PBC, n):
        bins = np.arange(0,rMax+step, step=step)
        hist = np.array([bins[:-1], np.zeros(len(bins)-1)])
        if PBC:
            expTraj = self.expandTraj(traj, n)
            for frame in range(len(traj)):
                frameHist = np.array(np.histogram(self.calcDistancesPBC(traj[frame], expTraj[frame], ), bins=bins))
                hist[1] = hist[1] + frameHist[0]
        else:          
            for frame in traj:
                frameHist = np.array(np.histogram(self.calcDistances(frame), bins=bins))
                hist[1] = hist[1] + frameHist[0]
            #mid_bins = (bins[:-1] + bins[1:])/2
        return hist

    def calculate_g(self, rMax, step, first, last, n,PBC=False):
        self.calc_rho() 
        hist = self.generateHist(self.traj[first:last], rMax, step, PBC, n)
        hist[1] = hist[1]/len(self.traj[first:last])
        hist[1] = hist[1]/len(self.traj[0])
        hist[1] = hist[1]/(self.rho*4*np.pi*step)
        #hist[1] = hist[1]/(self.rho*4/3*np.pi*(np.power(hist[0]+step, 3)-np.power(hist[0], 3)))
        hist[1] = hist[1]/np.power(hist[0],2)
        self.g = hist
        
    def calc_rho(self):
        v = self.box[0]*self.box[1]*self.box[2]
        self.rho = len(self.traj[0])/v
    
    def expandTraj(self,traj, n):
        newTraj = []
        for frame in traj:
            newTraj.append(self.expandFrame(frame, n))
        return np.array(newTraj)
            
            
    def block(self, frame, xyz):
        xyz=np.array(xyz)
        addition = xyz*self.box
        #print(addition)
        newBlock = frame + addition[np.newaxis]
        '''newBlock = np.zeros(frame.shape)
        for atom in range(len(frame)):
            newBlock[atom] = frame[atom]+ addition[np.newaxis]'''
        #print(frame.shape)
        #print(newBlock.shape)
        return newBlock
    
    def expandFrame(self, frame, n):
        disps = np.arange(-n, n+1, 1)
        newFrame = []
        for p in itertools.product(disps, disps, disps):
            if p != [0,0,0]:
                newFrame.append(self.block(frame, p))
        return np.concatenate(newFrame)
