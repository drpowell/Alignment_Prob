
package common;

import java.io.*;

import com.braju.format.*; // in hb15.zip  provides the Format.printf routines.

/** Despite its name, this class does _not_ implement a full 3 state FSM.
    It implements a 3 state FSM for the purpose of Linear gapped costs for
    a DPA.  There are therefore 5 parameters: 
              given the last operation was a diagonal (match or change): diag, start gap
              given the last operation was a gap: diag, start new gap, continue current gap

    So, the start_fromD encodes the start of a gap _and_ the first character
         2**-diag_fromD + 2 * 2**-start_fromD = 1         (2 possible ways to start a gap: ins or del)
         2**-diag_fromI + 2**-start_fromI + 2**-cont_fromI = 1

    The match/mismatch costs are part of the Two_Seq_Model
*/

public class FSM_3State extends Mutation_FSM.With_Counts {

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
		"diag_fromI="+diag_fromI+" start_fromI="+start_fromI+" cont_fromI="+cont_fromI+"\n";
	}
    }

    FSM_Params p;

    double dval, hval, vval;
    Counts d_counts, h_counts, v_counts;

    final private int diag_fromD=0, start_fromD=1,
	diag_fromI=2, start_fromI=3, cont_fromI=4;

    /*
    public FSM_3State(Two_Seq_Model_Counts s, double diag_cost, double start_gap,
		    double cont_gap, int numCounts, int countIndex) {
	super(numCounts, countIndex);
	p = new FSM_Params(s, diag_cost, start_gap, cont_gap);
	reset();
    }
    */

    public FSM_3State(Two_Seq_Model_Counts s, Params par, int numCounts, int countIndex) {
	super(numCounts, countIndex);
	p = new FSM_Params(s, par);
	reset();
    }
    
    public FSM_3State(FSM_Params p, int numCounts, int countIndex) { 
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

	counts = null;
    }


    public static int required_counts() { return 5; };
    public Params counts_to_params(Counts c) { 
	Params par = new Params();
	par.put("diag_fromD", -MyMath.log2(c.counts[countIndex+diag_fromD]));
	par.put("start_fromD", -MyMath.log2(c.counts[countIndex+start_fromD]));
	par.put("diag_fromI", -MyMath.log2(c.counts[countIndex+diag_fromI]));
	par.put("start_fromI", -MyMath.log2(c.counts[countIndex+start_fromI]));
	par.put("cont_fromI", -MyMath.log2(c.counts[countIndex+cont_fromI]));
	par.join( p.s.counts_to_params(c) );

	// Normalise the costs, and replace them in 'par'
	FSM_Params new_p = new FSM_Params(p.s, par);
	par.put("diag_fromD",  new_p.diag_fromD);
	par.put("start_fromD", new_p.start_fromD); 
	par.put("diag_fromI",  new_p.diag_fromI);
	par.put("start_fromI", new_p.start_fromI);
	par.put("cont_fromI",  new_p.cont_fromI);

	return par;
    };
    public Counts get_counts() { 
	// Combine counts from each of the three states weighted by their val
	Counts c = (Counts)d_counts.clone();
	double val = dval;

	c.combine_with_lens(val, h_counts, hval);
	val = MyMath.logplus(val, hval);

	c.combine_with_lens(val, v_counts, vval);

	return c; 
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

    public Object clone() {
	return new FSM_3State(p, numCounts, countIndex);
    }

    public void init_val(double v) { dval=v; hval = vval = MyMath.Big_Double; };
    public double get_val() { 
	double v = MyMath.logplus(dval, hval); 
	v = MyMath.logplus(v, vval);
	return v;
    };
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

    public void calc(Mutation_FSM h, Mutation_FSM v, Mutation_FSM d, char a, char b, int i, int j) {
        FSM_3State hcell = (FSM_3State)h;
        FSM_3State vcell = (FSM_3State)v;
        FSM_3State dcell = (FSM_3State)d;

	if (hcell != null) {
	    double char_cost = p.s.encB(b,j);
	    double w;
	    double count_inc;

	    // From 'd' state
	    w = dval + char_cost + p.start_fromD;
	    count_inc = 1.0/(1.0+MyMath.exp2(w-hcell.hval));
	    if (Double.isNaN(count_inc)) count_inc=0;
	    hcell.or_h(w, d_counts);
	    hcell.p.s.update_count_encB(hcell.h_counts, count_inc, b, j);
	    hcell.h_counts.inc(countIndex+start_fromD, count_inc);

	    // From 'v' state
	    w = vval + char_cost + p.start_fromI;
	    count_inc = 1.0/(1.0+MyMath.exp2(w-hcell.hval));
	    if (Double.isNaN(count_inc)) count_inc=0;
	    hcell.or_h(w, v_counts);
	    hcell.p.s.update_count_encB(hcell.h_counts, count_inc, b, j);
	    hcell.h_counts.inc(countIndex+start_fromI, count_inc);

	    // From 'h' state
	    w = hval + char_cost + p.cont_fromI;
	    count_inc = 1.0/(1.0+MyMath.exp2(w-hcell.hval));
	    if (Double.isNaN(count_inc)) count_inc=0;
	    hcell.or_h(w, h_counts);
	    hcell.p.s.update_count_encB(hcell.h_counts, count_inc, b, j);
	    hcell.h_counts.inc(countIndex+cont_fromI,  count_inc);
	}

	if (vcell != null) {
	    double char_cost = p.s.encA(a,i);
	    double w;
	    double count_inc;

	    // From 'd' state
	    w = dval + char_cost + p.start_fromD;
	    count_inc = 1.0/(1.0+MyMath.exp2(w-vcell.vval));
	    if (Double.isNaN(count_inc)) count_inc=0;
	    vcell.or_v(w, d_counts);
	    vcell.p.s.update_count_encA(vcell.v_counts, count_inc, a, i);
	    vcell.v_counts.inc(countIndex+start_fromD,  count_inc);

	    // From 'v' state
	    w = vval + char_cost + p.cont_fromI;
	    count_inc = 1.0/(1.0+MyMath.exp2(w-vcell.vval));
	    if (Double.isNaN(count_inc)) count_inc=0;
	    vcell.or_v(w, v_counts);
	    vcell.p.s.update_count_encA(vcell.v_counts, count_inc, a, i);
	    vcell.v_counts.inc(countIndex+cont_fromI,  count_inc);

	    // From 'h' state
	    w = hval + char_cost + p.start_fromI;
	    count_inc = 1.0/(1.0+MyMath.exp2(w-vcell.vval));
	    if (Double.isNaN(count_inc)) count_inc=0;
	    vcell.or_v(w, h_counts);
	    vcell.p.s.update_count_encA(vcell.v_counts, count_inc, a, i);
	    vcell.v_counts.inc(countIndex+start_fromI,  count_inc);
	}

	if (dcell != null) {
	    double char_cost = p.s.encBoth(a,b,i,j);
	    double w;
	    double count_inc;

	    // From 'd' state
	    w = dval + char_cost + p.diag_fromD;
	    count_inc = 1.0/(1.0+MyMath.exp2(w-dcell.dval));
	    if (Double.isNaN(count_inc)) count_inc=0;
	    dcell.or_d(w, d_counts);
	    dcell.p.s.update_count_encBoth(dcell.d_counts, count_inc, a, b, i, j);
	    dcell.d_counts.inc(countIndex+diag_fromD,  count_inc);

	    // From 'v' state
	    w = vval + char_cost + p.diag_fromI;
	    count_inc = 1.0/(1.0+MyMath.exp2(w-dcell.dval));
	    if (Double.isNaN(count_inc)) count_inc=0;
	    dcell.or_d(w, v_counts);
	    dcell.p.s.update_count_encBoth(dcell.d_counts, count_inc, a, b, i, j);
	    dcell.d_counts.inc(countIndex+diag_fromI,  count_inc);

	    // From 'h' state
	    w = hval + char_cost + p.diag_fromI;
	    count_inc = 1.0/(1.0+MyMath.exp2(w-dcell.dval));
	    if (Double.isNaN(count_inc)) count_inc=0;
	    dcell.or_d(w, h_counts);
	    dcell.p.s.update_count_encBoth(dcell.d_counts, count_inc, a, b, i, j);
	    dcell.d_counts.inc(countIndex+diag_fromI,  count_inc);
	}

    };

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

    /** or() - a new transition into this cell.
	All new transitions start in the diag state, so just call or_d
    */
    public void or(double d, Counts c) {
	or_d(d, c);
    }

    public void or_d(double d, Counts c) {
	// Update the counts first
	d_counts.combine_with_lens(dval, c, d);
	dval = MyMath.logplus(dval, d);
    }

    public void or_h(double d, Counts c) {
	// Update the counts first
	h_counts.combine_with_lens(hval, c, d);
	hval = MyMath.logplus(hval, d);
    }

    public void or_v(double d, Counts c) {
	// Update the counts first
	v_counts.combine_with_lens(vval, c, d);
	vval = MyMath.logplus(vval, d);
    }
}

