
package common;

import java.io.*;

import com.braju.format.*; // in hb15.zip  provides the Format.printf routines.

public final class Counts implements Serializable {
    public double counts[];
    public int num;
    
    public Counts(int num) {
	counts = new double[num];
	this.num = num;
    }

    public double get(int i) { return counts[i]; };
    
    public void zero() {
	for (int i=0; i<num; i++) { counts[i] = 0; }
    }

    public void combine_with_lens(double myLen, Counts c, double otherLen) {
	double min = MyMath.min2(myLen, otherLen);
	double w1 = (min==myLen    ? 1 : MyMath.exp2(min-myLen) );
	double w2 = (min==otherLen ? 1 : MyMath.exp2(min-otherLen) );
	scale(w1);
	linearWeight(c, w2);
	scale(1.0/(w1+w2));
    }

    public void linearWeight(Counts c, double w) {
	for (int i=0; i<num; i++) {
	    counts[i] += w * c.counts[i];
	}
    }

    public void inc(int index, double w) {
	counts[index] += w;
    }

    public void scale(double w) {
	for (int i=0; i<num; i++)
	    counts[i] *= w;
    }

    public void duplicate(Counts c) {
	for (int i=0; i<num; i++)
	    counts[i] = c.counts[i];
    }

    public Object clone() {
	Counts c = new Counts(num);
	c.duplicate(this);
	return c;
    }

    public String toString() {
	String r = new String("");
	for (int i=0; i<num; i++) {
	    r = Format.sprintf("%s %d:%.3f", new Parameters(r).add(i).add(counts[i]));
	}
	return r;
    }
}
