def readBoundaryNodes(filename:str):
    nodes = []
    try:
        with open(filename,'r') as f:
            lines = f.readlines()
            for line in lines:
                ll = line.split()
                vertex_coordinate = (float(ll[0]), float(ll[1]))
                if vertex_coordinate in nodes:
                    continue
                nodes.append(vertex_coordinate)
    except Exception as e:
        print(e)
        raise FileNotFoundError(" The file does not exist.")
    return nodes
