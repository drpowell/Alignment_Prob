
package common;

/** 
    An arbitrary order Markov Model for an arbitrary alphabet of characters.
    This meant to be used when the transition probabilities are known.

    If order == -1 then we are a uniform model over the alphabet
*/

public class MarkovN_exact implements Seq_Model {
    char[] chars;
    double[] probs;
    int order;
    
    StringBuffer past;

    public MarkovN_exact(int order, char[] chars) {
        Misc.assert(order>=-1, "Bad order="+order);

        this.chars = chars;
        this.order = order;

        if (order>=0) {
	    probs = new double[ (int)Math.pow(chars.length,order+1) ];
            for (int i=0; i<probs.length; i++) probs[i] = 1.0/chars.length;
        }

        past = new StringBuffer();
    }

    protected String num2Chars(int n) {
	StringBuffer res = new StringBuffer();
	for (int i=0; i<order+1; i++) {
	    res.append( chars[ n%chars.length ] );
	    n = n/chars.length;
	}
	return res.reverse().toString();
    }

    protected int chars2Num(String c) {
        int res=0, i;
        for (i=0; i<c.length(); i++) {
            int j;
            for (j=0; j<chars.length; j++)
                if (chars[j] == c.charAt(i)) break;
            Misc.assert(j<chars.length, "Character '"+c.charAt(i)+"' is unexpected");
            res = (res*chars.length) + j;
        }
        //      System.out.println("char2Num("+c+")="+res);
        return res;
    }

    public void setProb(String c, double p) {
	Misc.assert(c.length()-1 == order, "Bad length of chars passed to setProb: c="+c);
	if (order<0) return;
	probs[ chars2Num( c ) ] = p;
    }

    public void normalise() {
	for (int n=0; n < probs.length; n+=chars.length) {
	    double sum=0;
	    for (int i=0; i<chars.length; i++) sum += probs[n+i];
	    for (int i=0; i<chars.length; i++) probs[n+i] /= sum;
	}
    }

    public double encodeLen(char a, int i) {
        if (past.length() < order || order<0) return -MyMath.log2((double)1.0/chars.length);
        double p = probs[ chars2Num(past.substring(i-order) + a) ];
        return -MyMath.log2(p);
    }

    public double update(char a, int i) {
        double res = encodeLen(a, i);
        past.append(a);
        return res;
    }

    public String toString() {
	StringBuffer res = new StringBuffer();
	res.append("Markov_exact:  probs:\n");
	if (order<0) return res.toString();
	for (int i=0; i < probs.length; i++) {
	    String s = num2Chars(i);
	    res.append("p[" + s.substring(order) + " | " + s.substring(0, order) + "]");
	    res.append(" = " + probs[i] + "\n");
	}
	return res.toString();
    }


    public static void main(String args[]) {
        String s = args[0];
        char[] a = {'a','t','g','c'};
        MarkovN_exact m = new MarkovN_exact(1, a);

	m.setProb("aa", 1);
	m.setProb("at", 1);
	m.setProb("ag", 2);
	m.setProb("ac", 1);
	m.normalise();

	System.out.println(m);

        double tot = 0;
        for (int i=0; i<s.length(); i++) {
            double r = m.update(s.charAt(i), i);
            tot += r;
            System.out.println(r);
        }
        System.out.println("Total entropy = "+tot+" bits/ch");
    }
 
}
