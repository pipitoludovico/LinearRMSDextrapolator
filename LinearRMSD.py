import MDAnalysis as Mda
import matplotlib.pyplot as plt
from sys import argv


class Getter:
    def __init__(self, file, traj, topology):
        self.file = file
        self.traj = traj
        self.topology = topology
        self.rmsd_values: list = []
        self.val: dict = {}
        self.sorted_keys = None
        self.RMSDS, self.FRAMES = [], []
        self.Limit = 100
        self.decimals = 0
        self.GetMetric()
        self.GetUniquesRMSDS()
        self.GetRMSDandFRAME()
        self.WriteFrames()
        self.GetTrajectory()
        if '-plot' in [str(arg.lower()) for arg in argv]:
            try:
                sliceStr = argv[5]
                start, end = map(int, sliceStr.split(':'))
                self.Limit = slice(start, end)
            except:
                print('\nNo slice defined in the -plot arguments: using the default 0:50')
                self.Limit = slice(0, 50)
            self.Plot(self.Limit)

    def __str__(self):
        for key in self.sorted_keys:
            print("RMSD ", key, "Frame ", self.val[key])
        return ""

    def WriteFrames(self):
        with open("frame_list.csv", 'w') as flist:
            flist.write("FRAME,RMSD\n")
            for key in self.sorted_keys:
                flist.write(f"{str(self.val[key])},{int(key)}\n")

    def GetMetric(self) -> None:
        with open(self.file, 'r') as outVmd:
            for line in outVmd.readlines():
                self.rmsd_values.append(float(line.split("\t")[1]))

    def GetUniquesRMSDS(self) -> None:
        for idx, rmsd in enumerate(self.rmsd_values):
            rmsd_ = round(self.rmsd_values[idx], self.decimals)
            if rmsd_ not in self.val:
                self.val[rmsd_] = idx
        self.sorted_keys = sorted(self.val, key=lambda k: float(k))

    def GetRMSDandFRAME(self) -> (list, list):
        for key in self.sorted_keys:
            self.RMSDS.append(key)
            self.FRAMES.append(self.val[key])

    def GetTrajectory(self):
        Universe = Mda.Universe(self.topology, self.traj)
        atomsel = Universe.select_atoms('all')
        frameSelection = [ts for ts in self.FRAMES]
        atomsel.write(f'merged.xtc', frames=Universe.trajectory[frameSelection])

    def Plot(self, Limit):
        plt.plot(range(len(self.RMSDS[Limit])), self.RMSDS[Limit], marker='o', linestyle='-')
        plt.xlabel('Data Point')
        plt.ylabel('RMSD')
        plt.title('Frame/RMSD extrapolated from the trajectory')
        plt.grid(True)
        plt.savefig(f'RMSD_DATAPOINT.png')


if len(argv) < 4:
    print("example usage: python LinearRMSD.py RMSDs_from_VMD.dat traj.xtc structure.psf -plot 0:61")
    exit()
else:
    getter = Getter(argv[1], argv[2], argv[3])
    print(getter)
