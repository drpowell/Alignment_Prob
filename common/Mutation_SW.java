
package common;

import java.io.*;

import com.braju.format.*; // in hb15.zip  provides the Format.printf routines.

/**
   Implement Smith-Waterman costs.
*/

public abstract class Mutation_SW extends Mutation_FSM {

    private class FSM_Params implements Serializable {
	Two_Seq_Model_Counts s;

	double diag_fromD, start_fromD;
	double diag_fromI, start_fromI, cont_fromI;
	/*
	FSM_Params(Two_Seq_Model_Counts s, double diag_cost, double start_gap, double cont_gap) {
	    this.s = s;

	    normalize_costs(diag_cost, start_gap, cont_gap);
	}
	*/

	FSM_Params(Two_Seq_Model_Counts s, Params p) {
	    this.s = s;

	    if (!p.exists("diag_fromD")) {
		set_default_costs();
	    } else {
		diag_fromD  = p.get("diag_fromD");
		start_fromD = p.get("start_fromD");
		diag_fromI  = p.get("diag_fromI");
		start_fromI = p.get("start_fromI");
		cont_fromI  = p.get("cont_fromI");
		
		normalize_costs();
	    }
	}

	void set_default_costs() {
	    diag_fromD  = 0;
	    start_fromD = 3;

	    diag_fromI  = 0;
	    start_fromI = 3;
	    cont_fromI  = 1;

	    normalize_costs();
	}

	void normalize_costs() {
	    double sum;

	    sum = MyMath.exp2(-diag_fromD) + 2*MyMath.exp2(-start_fromD);
	    diag_fromD  = diag_fromD  + MyMath.log2( sum );
	    start_fromD = start_fromD + MyMath.log2( sum );

	    sum = MyMath.exp2(-diag_fromI) + MyMath.exp2(-start_fromI) + MyMath.exp2(-cont_fromI);
	    diag_fromI  = diag_fromI  + MyMath.log2( sum );
	    start_fromI = start_fromI + MyMath.log2( sum );
	    cont_fromI  = cont_fromI  + MyMath.log2( sum );
	}

        public String toString() {
	    return "diag_fromD="+diag_fromD+" start_fromD="+start_fromD+"\n"+
		"diag_fromI="+diag_fromI+" start_fromI="+start_fromI+" cont_fromI="+cont_fromI;
	}
    }

    FSM_Params p;

    double dval, hval, vval;
    Counts d_counts, h_counts, v_counts;

    final protected int diag_fromD=0, start_fromD=1,
	diag_fromI=2, start_fromI=3, cont_fromI=4;

    public Mutation_SW(Two_Seq_Model_Counts s, Params par, int numCounts, int countIndex) {
	super(numCounts, countIndex);
	p = new FSM_Params(s, par);
	reset();
    }
    
    public Mutation_SW(FSM_Params p, int numCounts, int countIndex) { 
	super(numCounts, countIndex);
	this.p = p; 
	reset();
    }

    public void reset() {
	super.reset();
	dval = MyMath.Big_Double;
	hval = MyMath.Big_Double;
	vval = MyMath.Big_Double;

	if (d_counts==null) d_counts = new Counts(numCounts);
	if (h_counts==null) h_counts = new Counts(numCounts);
	if (v_counts==null) v_counts = new Counts(numCounts);
	d_counts.zero();
	h_counts.zero();
	v_counts.zero();

	//counts = null;
    }


    public static int required_counts() { return 5; };

    // counts_to_params - convert counts into parameters.  
    //   Call the Two_Seq_Model to convert its own.
    //   Note the 0.5* in the start_fromD computation.  This is cause there are 2 start_fromD 
    //   arcs out of each state, but we only have one count for them combined.
    public Params counts_to_params(Counts c) { 
	Params par = new Params();
	double sum1 = c.counts[countIndex+diag_fromD] + c.counts[countIndex+start_fromD];
	double sum2 = c.counts[countIndex+diag_fromI] + c.counts[countIndex+start_fromI] +
	    c.counts[countIndex+cont_fromI];

	par.put("diag_fromD", -MyMath.log2(c.counts[countIndex+diag_fromD]/sum1));
	par.put("start_fromD", -MyMath.log2(0.5*c.counts[countIndex+start_fromD]/sum1));
	par.put("diag_fromI", -MyMath.log2(c.counts[countIndex+diag_fromI]/sum2));
	par.put("start_fromI", -MyMath.log2(c.counts[countIndex+start_fromI]/sum2));
	par.put("cont_fromI", -MyMath.log2(c.counts[countIndex+cont_fromI]/sum2));
	par.join( p.s.counts_to_params(c) );

	return par;
    };

    public double encode_params() {
	Counts c = get_counts();
	double[] probs;
	double dataLen;

	// First encode match/change paramters. Pass the number of these events to encode_params
	double len = p.s.encode_params(c.counts[countIndex+diag_fromD] + c.counts[countIndex+diag_fromI]);

	// 'probs' is the paramaters of the multinomial distribution.
	// 'dataLen' is the number of things in this multinomial.

	probs = new double[] { MyMath.exp2(-p.diag_fromD), 2*MyMath.exp2(-p.start_fromD) };
	dataLen = c.counts[countIndex+diag_fromD] + c.counts[countIndex+start_fromD];

	len += Multinomial.MMLparameter_cost(probs, dataLen);

	probs = new double[] { MyMath.exp2(-p.diag_fromI), MyMath.exp2(-p.start_fromI), 
			       MyMath.exp2(-p.cont_fromI) };
	dataLen = c.counts[countIndex+diag_fromI] + c.counts[countIndex+start_fromI] + 
	    c.counts[countIndex+start_fromI];

	len += Multinomial.MMLparameter_cost(probs, dataLen);

	return len;
    }

    public double alignmentLength() {
	Counts c = get_counts();
	return c.counts[countIndex+diag_fromD] + 
	       c.counts[countIndex+start_fromD] + 
	       c.counts[countIndex+diag_fromI] + 
	       c.counts[countIndex+start_fromI] + 
	       c.counts[countIndex+cont_fromI];
    }

    public void init_val(double v) { dval=v; hval = vval = MyMath.Big_Double; };

    public void normalise(double v) { dval -= v; hval -= v; vval -= v; }

    /** Initialise diag counts */
    public void init_counts(Counts c) {
	d_counts.duplicate(c);
    }


    public String paramsToString() {
	return this.getClass() + ": " + p + "\n";
    }

    public String toString() {
	return 
	    "dval="+dval+" hval="+hval+" vval="+vval+"\n"+
	    "d_counts="+d_counts+"\n"+
	    "h_counts="+h_counts+"\n"+
	    "v_counts="+v_counts+"\n";
    }

    /** add() - something extra must be encoded in this state.  
	Update all three states
    */
    public void add(double v, int cIndex)  { 
	dval += v;
	hval += v;
	vval += v;
	if (cIndex>=0) {
	    d_counts.inc(cIndex, 1);
	    h_counts.inc(cIndex, 1);
	    v_counts.inc(cIndex, 1);
	}
    }

    public static class One extends Mutation_SW implements TraceBack_Info {
	
	private static class TB_Data extends TraceBack_Data {
	    int state;
	    TB_Data() {
		super();
		state=-1;
	    }
	    TB_Data(TraceBack_Data td, int s) {
		super();
		i = td.i;
		j = td.j;
		state=s;
	    }
	}

	TB_Data d_id, v_id, h_id;
	TB_Data d_from, v_from, h_from;

	final static double startGap = 16;
	final static double contGap = 4;
	final static double match = -5;
	final static double mismatch = 4;

	public One(FSM_Params p, int numCounts, int countIndex) { 
	    super(p, numCounts, countIndex);
	    d_id = v_id = h_id = null;
	    d_from = v_from = h_from = null;
	}
	
	public One(Two_Seq_Model_Counts s, Params par, int numCounts, int countIndex) {
	    super(s, par, numCounts, countIndex);
	    d_id = v_id = h_id = null;
	    d_from = v_from = h_from = null;
	}

	public void set_tbdata(TraceBack_Data id) {  
	    d_id = new TB_Data(id, 0); 
	    v_id = new TB_Data(id, 1); 
	    h_id = new TB_Data(id, 2); 
	}

	public TraceBack_Data get_tbdata()        {  
	    if (dval <= hval && dval <= vval) {
		return d_id;
	    } else if (vval <= dval && vval <= hval) {
		return v_id;
	    } else if (hval <= dval && hval <= vval) {
		return h_id;
	    }
	    Misc.assert(false, "Bad vals");
	    return null;    
	}
	
	public TraceBack_Data get_from(TraceBack_Data td) { 
	    TB_Data t = (TB_Data)td;
	    Misc.assert(t.i==d_id.i && t.j==d_id.j, "Not my traceback data!");
	    switch(t.state) {
	    case 0: return d_from;
	    case 1: return v_from;
	    case 2: return h_from;
	    default:
		Misc.assert(false, "Bad trackback data. t.state="+t.state);
	    }
	    return null;
	}

	public Object clone() {
	    return new Mutation_SW.One(p, numCounts, countIndex);
	}

	public double get_val() { 
	    return MyMath.min3(dval, vval, hval);
	};

	public Counts get_counts() { 
	    // Get counts from state with smallest val
	    if (dval <= hval && dval <= vval)
		return d_counts;

	    if (vval <= dval && vval <= hval)
		return v_counts;

	    if (hval <= dval && hval <= vval)
		return h_counts;

	    Misc.assert(false, "Bad vals!");
	    return null;
	};

	public void calc(Mutation_FSM h, Mutation_FSM v, Mutation_FSM d, 
			 char a, char b, int i, int j) {
	    Mutation_SW.One hcell = (Mutation_SW.One)h;
	    Mutation_SW.One vcell = (Mutation_SW.One)v;
	    Mutation_SW.One dcell = (Mutation_SW.One)d;

	    Counts tcounts = (Counts)counts.clone();
	    if (hcell != null) {
		hcell.or_h(dval + startGap, d_id);
		hcell.or_h(vval + startGap, v_id);
		hcell.or_h(hval + contGap,  h_id);
	    }

	    if (vcell != null) {
		vcell.or_v(dval + startGap, d_id);
		vcell.or_v(vval + contGap,  v_id);
		vcell.or_v(hval + startGap, h_id);
	    }

	    if (dcell != null) {
		double char_cost = (a==b ? match : mismatch);
		dcell.or_d(dval + char_cost, d_id);
		dcell.or_d(vval + char_cost, v_id);
		dcell.or_d(hval + char_cost,  h_id);
	    }
	}

	public void or(double d, Counts c) {
	    or(d, c, null);
	}

	/** or() - a new transition into this cell.
	    All new transitions start in the diag state, so just call or_d
	*/
	public void or(double d, Counts c, Mutation_FSM from) {
	    TB_Data f = null;
	    if (from!=null) f = (TB_Data)(((One)from).get_tbdata());
	    or_d(d, f);
	}


	public void or_d(double d, TB_Data from_id) {
	    if (d<dval) {
		dval = d;
		d_from = from_id;
	    }
	}

	public void or_h(double d, TB_Data from_id) {
	    if (d<hval) {
		hval = d;
		h_from = from_id;
	    }
	}

	public void or_v(double d, TB_Data from_id) {
	    if (d<vval) {
		vval = d;
		v_from = from_id;
	    }
	}
    }

}

