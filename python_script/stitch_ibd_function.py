def stitch_ibd(ts):
	Pair = []																
	Left = np.array([])
	Right = np.array([])
	Node = np.array([])
	Time = np.array([])
	for k, v in segments.items():
		kLeft = v.left
		kRight = v.right
		kNode = v.node
		order = kLeft.argsort()
		l = kLeft[order]
		n = kNode[order]
		r = kRight[order]
		node = (np.convolve(n + (n[0] == 0), [-1, 1]) != 0)[:-1]	
		discontinuities = np.concatenate(([True], l[1:] != r[:-1]))
		i = np.where(node | discontinuities)[0]
		kPair = [k, ]*len(i)
		kLeft = l[i]
		kRight = np.append(r[i[1:]-1], r[len(r)-1])
		kNode = n[i] 
		Pair = Pair + kPair
		Left = np.append(Left, kLeft)
		Right = np.append(Right, kRight)	
		Node = np.append(Node, kNode)	
		Time = np.append(Time, ts.tables.nodes.time[kNode])
	Length = np.subtract(Right,Left)	
	Pair = [(int(p[0]), int(p[1])) for p in Pair]
	dfIBD = pd.DataFrame(columns = ['Pair', 'Hap1', 'Hap2', 'Ind1', 'Ind2', 'Left', 'Right', 'Length', 'Node', 'Time'])
	dfIBD['Pair'] = Pair                                  														
	dfIBD['Left'] = Left
	dfIBD['Right'] = Right 
	dfIBD['Length'] = Length 
	dfIBD['Node'] = Node
	dfIBD['Time'] = Time
	return dfIBD
