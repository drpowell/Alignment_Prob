
package common;

import java.io.*;

import com.braju.format.*; // in hb15.zip  provides the Format.printf routines.

public class FSM_Prob extends Mutation_FSM.With_Counts {

    private class FSM_Prob_Params implements Serializable {
	double ins_cost, del_cost, diag_cost;
	Two_Seq_Model_Counts s;

	FSM_Prob_Params(Two_Seq_Model_Counts s, double ins_cost, double del_cost, double diag_cost) {
	    this.s = s;
	    this.ins_cost = ins_cost;
	    this.del_cost = del_cost;
	    this.diag_cost = diag_cost;

	    normalize_costs();
	}

	FSM_Prob_Params(Two_Seq_Model_Counts s, Params p) {
	    this.s = s;

	    if (!p.exists("ins_cost")) {
		set_default_costs();
	    } else {
		ins_cost = p.get("ins_cost");
		del_cost = p.get("del_cost");
		diag_cost = p.get("diag_cost");
	    }

	    normalize_costs();
	}

	void set_default_costs() {
	    ins_cost  = -MyMath.log2(0.05);
	    del_cost  = -MyMath.log2(0.05);
	    diag_cost = -MyMath.log2(0.9);
	}

	void normalize_costs() {
	    double sum = MyMath.exp2(-ins_cost)+MyMath.exp2(-del_cost)+MyMath.exp2(-diag_cost);
	    ins_cost  = ins_cost+MyMath.log2( sum );
	    del_cost  = del_cost+MyMath.log2( sum );
	    diag_cost = diag_cost+MyMath.log2( sum );
	    //	    Format.printf("ins_cost=%.5f del_cost=%.5f diag_cost=%.5f\n", 
	    //			  new Parameters(ins_cost).add(del_cost).add(diag_cost));
	}

        public String toString() {
	    return "diag_cost="+diag_cost+" ins_cost="+ins_cost+" del_cost="+del_cost+"\n";
	}
    }


    FSM_Prob_Params p;

    double val;

    final private int insIndex=0, delIndex=1, diagIndex=2;

    public FSM_Prob(Two_Seq_Model_Counts s, double ins_cost, double del_cost,
		    double diag_cost, int numCounts, int countIndex) {
	super(numCounts, countIndex);
	p = new FSM_Prob_Params(s, ins_cost, del_cost, diag_cost);
	reset();
    }

    public FSM_Prob(Two_Seq_Model_Counts s, Params par, int numCounts, int countIndex) {
	super(numCounts, countIndex);
	p = new FSM_Prob_Params(s, par);
	reset();
    }
    
    public FSM_Prob(FSM_Prob_Params p, int numCounts, int countIndex) { 
	super(numCounts, countIndex);
	this.p = p; 
	reset();
    }

    public void reset() {
	super.reset();
	val = MyMath.Big_Double;
	counts.zero();
    }


    public static int required_counts() { return 3; };
    public Params counts_to_params(Counts c) { 
	double sum = c.counts[countIndex+insIndex] + 
	    c.counts[countIndex+delIndex] + 
	    c.counts[countIndex+diagIndex];
	Params par = new Params();
	par.put("ins_cost", -MyMath.log2(c.counts[countIndex+insIndex]/sum));
	par.put("del_cost", -MyMath.log2(c.counts[countIndex+delIndex]/sum));
	par.put("diag_cost", -MyMath.log2(c.counts[countIndex+diagIndex]/sum));
	par.join( p.s.counts_to_params(c) );
	return par;
    };

    public double encode_params() {
	Counts c = get_counts();

	// First encode match/change paramters. Pass the number of these events to encode_params
	double len = p.s.encode_params(c.counts[countIndex+diagIndex]);

	double dataLen = c.counts[countIndex+diagIndex] + c.counts[countIndex+insIndex] + 
	    c.counts[countIndex+delIndex];

	len += Multinomial.MMLparameter_cost(new double[] {
	    MyMath.exp2(-p.ins_cost), MyMath.exp2(-p.del_cost), MyMath.exp2(-p.diag_cost) },
					     dataLen);
	return len;
    }

    public double alignmentLength() {
	Counts c = get_counts();
	return c.counts[countIndex+insIndex] + 
	       c.counts[countIndex+delIndex] + 
	       c.counts[countIndex+diagIndex];
    }

    public Object clone() {
	return new FSM_Prob(p, numCounts, countIndex);
    }

    public void init_val(double v) { val=v; };
    public double get_val() { return val; };
    public void normalise(double v) { val -= v; }

    public String paramsToString() {
	return this.getClass() + ": " + p + "\n";
    }

    public String toString() {
	return "val="+val+" counts="+counts+"\n";
    }

    public void calc(Mutation_FSM h, Mutation_FSM v, Mutation_FSM d, char a, char b, int i, int j) {
        FSM_Prob hcell = (FSM_Prob)h;
        FSM_Prob vcell = (FSM_Prob)v;
        FSM_Prob dcell = (FSM_Prob)d;

	if (hcell != null) {
	    double w = val + p.ins_cost  + p.s.encB(b,j);
            double count_inc = 1.0/(1.0+MyMath.exp2(w-hcell.val));
            if (Double.isNaN(count_inc)) count_inc=0;
	    hcell.or(w, counts);
	    hcell.p.s.update_count_encB(hcell.counts, count_inc, b, j);
	    hcell.counts.inc(countIndex+insIndex,  count_inc);
	}

	if (vcell != null) {
	    double w = val + p.del_cost  + p.s.encA(a,i);
            double count_inc = 1.0/(1.0+MyMath.exp2(w-vcell.val));
            if (Double.isNaN(count_inc)) count_inc=0;
	    vcell.or(w, counts);
	    vcell.p.s.update_count_encA(vcell.counts, count_inc, a, i);
	    vcell.counts.inc(countIndex+delIndex,  count_inc);
	}

	if (dcell != null) {
	    double w = val + p.diag_cost  + p.s.encBoth(a,b,i,j);
            double count_inc = 1.0/(1.0+MyMath.exp2(w-dcell.val));
            if (Double.isNaN(count_inc)) count_inc=0;
	    dcell.or(w, counts);
	    dcell.p.s.update_count_encBoth(dcell.counts, count_inc, a, b, i, j);
	    dcell.counts.inc(countIndex+diagIndex,  count_inc);
	}
    };

    public void add(double v, int cIndex) {
	val += v; 
	if (cIndex>=0)
	    counts.inc(cIndex, 1);
    };

    public void or(double d, Counts c) {
	// Update the counts first
	counts.combine_with_lens(val, c, d);
	val = MyMath.logplus(val, d);
    }
}

