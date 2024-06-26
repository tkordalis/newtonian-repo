####2
import unv
import binarytecplot as bt

mesh = unv.Reader        ("./remesh_39.3100.unv")
sol  = bt.LoadTecplotFile("./time_39.4000.plt")

time = sol.getZone().getName()



##<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
## Define Functions
##<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
def readBoundary(filename):

    Snodes = []
    for node in mesh.getNodeIDsThatBelongToGroup(filename):
        z = sol.getZone()["Z"][node]
        r = sol.getZone()["R"][node]
        Snodes.append( "{0:.8f} {1:.8f}".format(z,r) )
    return Snodes
##<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
##<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

with open("Bubble1_{}.dat".format(time) , "w") as f: [f.write(line + "\n" ) for line in readBoundary("Bubble1")]
with open("Bubble2_{}.dat".format(time) , "w") as f: [f.write(line + "\n" ) for line in readBoundary("Bubble2")]

