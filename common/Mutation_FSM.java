
package common;

import java.io.*;
import java.util.*;


public abstract class Mutation_FSM implements Serializable {
    class FSM_Params {
	Two_Seq_Model s;
	FSM_Params(Two_Seq_Model s) { this.s = s; };
    }

    public abstract void init_val(double v);
    public abstract double get_val();
    public abstract void normalise(double v);

    public abstract void calc(Mutation_FSM h, Mutation_FSM v, Mutation_FSM d, char a, char b, int i, int j);

    //    public abstract double enc_del(char a, int i);
    //    public abstract double enc_ins(char a, int i);
    //    public abstract double enc_diag(char a, char b, int i, int j);

    public abstract void or(double v);
    public abstract void add(double v);

    public String paramsToString() {
	return this.getClass() + ": paramsToString() not defined" + "\n";
    }

    public double encode_params() { 
	System.err.println("WARNING: encode_params() not implemented in " + this.getClass());
	return 0; 
    }

    public abstract Object clone();

    public abstract static class With_Counts extends Mutation_FSM {
	Counts counts;
	int numCounts;
	int countIndex;
	
	With_Counts(int numCounts, int countIndex) {
	    this.numCounts = numCounts;
	    counts = new Counts(numCounts);
	    this.countIndex = countIndex;
	}

	public void reset() {
	    counts.zero();
	}

	/** Initialise all counts */
	public void init_counts(Counts c) {
	    counts.duplicate(c);
	}
	
	public static int required_counts() { return 0; };
	public abstract Params counts_to_params(Counts c);
	
	public Counts get_counts() { return counts; };
	public double alignmentLength() {
	    System.err.println("WARNING: alignmentLength() not implemented in " + this.getClass());
	    return 0; 
	};
	
	public void add(double d) { System.err.println("method add(double) has no impl"); System.exit(-1); }
	public void or(double d)  { System.err.println("method or(double) has no impl"); System.exit(-1); }
	
	public abstract void add(double d, int cIndex);
	public abstract void or(double d, Counts c);
    }
}
