
package fuzzyLZ;

import common.*;

class Matches_Sparse extends FuzzyLZ.Matches {
    Sparse b,e;
    ExactMatches hash;
    int winSize;
    int computeWin;
    boolean fullN2, plotActive;
    int debug;
    boolean fwd;
    long hashHits;
    double cutML;
    double plotBrightness;

    double beginCost;

    class Sparse_Fsm implements Sparse.Obj {
	Mutation_FSM.With_Counts fsm;
	Sparse_Fsm(Mutation_FSM.With_Counts fsm) {
	    this.fsm = fsm;
	}

	public double get_val() {
	    return fsm.get_val();
	}

	public Object clone() {
	    return new Sparse_Fsm((Mutation_FSM.With_Counts)fsm.clone());
	}
    }

    Matches_Sparse(boolean fwd, Mutation_FSM.With_Counts fsm, Params p, int countIndex, char[] seq) {
        super(fsm,p,countIndex,seq);

	this.fwd = fwd;

	Sparse_Fsm sparse_fsm = new Sparse_Fsm(fsmType);
	b = new Sparse(sparse_fsm);
	e = new Sparse(sparse_fsm);
	hash       = null;

	winSize    = 25;
	computeWin = 10;
	cutML      = 4;

	plotBrightness = 0.5;
	
	fullN2     = false;
	plotActive = false;

	debug = 0;
	hashHits = 0;
    }

    void constructHash() {
	hash = new ExactMatches(sequence, winSize);
    }

    void setHash(ExactMatches h) {
	hash = h;
    }

    ExactMatches getHash() {
	return hash;
    }


    public void beginLinks(int i, double startLen, Counts startCounts) 
    {
	beginCost = startLen;

	// Do begin links into b,  (startLen & startCounts)
	if (debug>1) System.err.print("Active at ");
	for (Sparse.Iterate iter = b.moveFwd(null); iter.o != null; iter = b.moveFwd(iter)) {
	    Mutation_FSM.With_Counts cell = (Mutation_FSM.With_Counts)iter.o;
	    cell.or(startLen, startCounts);

	    if (plotActive) plot.putMax(i, iter.i, (fwd ? 0.5 : 0), (fwd ? 0 : 0.5), 0);

	    if (debug>1) System.err.print("("+i+","+iter.i+") ");
	}
	if (debug>1) System.err.println();
    }

    public double update(char aChar, int i, Counts retCounts) {
	if (hash==null && !fullN2)
	    constructHash();

	// Duplicate 'b' into 'e'.  Discard cells with bad msgLen, and chop where appropriate.
	if (debug>1) System.err.println("Before copy_cull_cut. Sparse = " + b);
	e.copy_cull_cut(b, (fwd ? 1 : -1), winSize, beginCost-cutML, computeWin);
	if (debug>1) System.err.println("After  copy_cull_cut. Sparse = " + e);

	// Add any new hash hits into 'e'.  'e' must not have useful data in it for this to work.
	if (i+1+winSize <= seqLen && !fullN2) {
	    ExactMatches.MyList l;
	    if (fwd)
		l = hash.get( sequence, i+1 );
	    else {
		String str = new String(sequence, i+1, winSize);
		l = hash.get( ExactMatches.Reverse_Complement_DNA.conv(str) );
	    }
	    if (l!=null)
		for (ExactMatches.MyList.L l2=l.start; l2!=null && l2.val<i; l2=l2.next) {
		    int hit = l2.val;
		    int n1  = MyMath.max2(0, hit-computeWin);
		    int n2  = MyMath.min2(i, hit+1+computeWin);
		    e.add(n1, n2);

		    if (plotActive) 
			for (int j=n1; j<n1+1; j++) 
			    plot.putMax(i, j, (fwd ? 1 : 0), (fwd ? 0 : 1), 0);
		    if (debug>0) System.err.println("New hash "+i+" = "+hit+" ("+n1+","+n2+")");

		    hashHits++;
		}
	}

	if (fullN2)
	    e.add(0, i);
	
	if (debug>1) System.err.println("Before join. Sparse = " + e);
	e.join();
	if (debug>1) System.err.println("After  join. Sparse = " + e);
	e.checkAlloc();

	Misc.assert(e.tail==null || e.tail.end <= i, "Cheating! Sparse end past current character");
	Misc.assert(e.head==null || e.head.start >= 0, "Bugger. Sparse start < 0");

	// Reset all 'e' cells.
	for (Sparse.Iterate eIter=e.moveFwd(0, null); eIter.o!=null; eIter=e.moveFwd(eIter))
	    ((Mutation_FSM.With_Counts)eIter.o).reset();

	{
	    // Compute 'e' from 'b' using FSM
	    Sparse.Iterate eIter = null;
	    Sparse.Iterate bIter = fwd ? b.moveFwd(null) : b.moveRev(null);
	    while (bIter.o != null) {
		Mutation_FSM.With_Counts cell,hcell,vcell,dcell;

		cell  = (Mutation_FSM.With_Counts)bIter.o;
		hcell = (Mutation_FSM.With_Counts)(fwd ? b.getNext(bIter) : b.getPrev(bIter));

		eIter = fwd ? e.moveFwd(bIter.i, eIter, true) : e.moveRev(bIter.i, eIter, true);
		vcell = (Mutation_FSM.With_Counts)eIter.o;
		dcell = (Mutation_FSM.With_Counts)(fwd ? e.getNext(eIter) : e.getPrev(eIter));
		
		int j = bIter.i;
		char bChar;
		if (fwd)
		    bChar = sequence[j]; 
		else 
		    bChar = (j>0 ? sequence[j-1] : '-');
		cell.calc(hcell, vcell, dcell, aChar, bChar, i, j);

		bIter = fwd ? b.moveFwd(bIter) : b.moveRev(bIter);
	    }
	}
	// Compute return base & Counts from 'e' and END_COPY
	double ret_msgLen = Double.POSITIVE_INFINITY;
	retCounts.zero();
	for (Sparse.Iterate eIter=e.moveFwd(0, null); eIter.o!=null; eIter=e.moveFwd(eIter)) {
	    Mutation_FSM.With_Counts cell  = (Mutation_FSM.With_Counts)eIter.o;

	    double endLen = cell.get_val() + encEnd;

	    // Hack: Increment the endCopy count in the _cell_  (will undo after next statement)
	    cell.get_counts().inc(countIndex + endIndex, 1);
	    retCounts.combine_with_lens(ret_msgLen, cell.get_counts(), endLen);
	    cell.get_counts().inc(countIndex + endIndex, -1);

	    ret_msgLen = MyMath.logplus(ret_msgLen, endLen);
	}

	// Compute 'e' to include CONT_COPY
	for (Sparse.Iterate eIter=e.moveFwd(0, null); eIter.o!=null; eIter=e.moveFwd(eIter)) {
	    Mutation_FSM.With_Counts cell  = (Mutation_FSM.With_Counts)eIter.o;

	    cell.add(encContinue, countIndex+contIndex);
	}

	// Increase the age of all active elements.
	e.incAge();

	// Swap 'e' and 'b'
	Sparse t = b;
	b = e;
	e = t;

	return ret_msgLen;
    }

    public double msgLen() {
	double res = Double.POSITIVE_INFINITY;
	for (Sparse.Iterate bIter=b.moveFwd(0, null); bIter.o!=null; bIter=b.moveFwd(bIter)) {
	    Mutation_FSM.With_Counts cell  = (Mutation_FSM.With_Counts)bIter.o;

	    res = MyMath.logplus(res, cell.get_val());
	}
	return res;
    }

    public double normalise(double base) {
	// Note: base must be less than _all_ vals in the cells.
	//       Use the result of msgLen() or something less
	double res = Double.POSITIVE_INFINITY;
	for (Sparse.Iterate bIter=b.moveFwd(0, null); bIter.o!=null; bIter=b.moveFwd(bIter)) {
	    Mutation_FSM.With_Counts cell  = (Mutation_FSM.With_Counts)bIter.o;
	    cell.normalise(base);
	    res = MyMath.logplus(res, cell.get_val());
	}
	return res;
    }


    public void plotVals(int i, double base) {
	if (plotActive) return;
	for (Sparse.Iterate bIter=b.moveFwd(0, null); bIter.o!=null; bIter=b.moveFwd(bIter)) {
	    Mutation_FSM.With_Counts cell  = (Mutation_FSM.With_Counts)bIter.o;
	    double v = MyMath.exp2(base-cell.get_val());
	    v = Math.pow(v, plotBrightness);
 	    plot.putMax(i, bIter.i, (fwd ? v : 0), (fwd ? 0 : v), 0);
	}
    }


    public String toString() {
	StringBuffer r = new StringBuffer();
	r.append((fwd ? "fwd: " : "rev: "));
	for (Sparse.Iterate bIter=b.moveFwd(0, null); bIter.o!=null; bIter=b.moveFwd(bIter)) {
	    Mutation_FSM.With_Counts cell  = (Mutation_FSM.With_Counts)bIter.o; 
	    r.append(bIter.i + ": " + cell.get_val() + " ");
	}
	return r.toString();
    }

    public void display_stats() {
	if (debug==0) return;
	System.err.print(fwd ? "fwd: " : "rev: ");
	System.err.println("Number of hash hits = "+hashHits);
	b.display_stats();
    }

}
