
package common;

/** 
    An arbitrary order Markov Model for an arbitrary alphabet of characters.
    The model is first fitted to the sequence.
    This is _not_ a proper model because the model cost is not accounted for.
    Thus the entropy of a sequence calculated by this class is below its real value.
    The discrepency become worse for higher order models.

    If order == -1 then we are a uniform model over the alphabet
*/
public class MarkovN_fitted implements Seq_Model {
    char[] chars;
    int[] charCounts;
    int[] countTotal;
    int order;
    
    StringBuffer past;

    public MarkovN_fitted(int order, char[] chars, String s) {
	Misc.assert(order>=-1, "Bad order="+order);

	this.chars = chars;
	this.order = order;

	if (order>=0) {
	    charCounts = new int[ (int)Math.pow(chars.length,order+1) ];
	    countTotal = new int[ (int)Math.pow(chars.length,order) ];
	    for (int i=0; i<charCounts.length; i++) charCounts[i]=1;
	    for (int i=0; i<countTotal.length; i++) countTotal[i]=chars.length;
	}

	past = new StringBuffer();

	// Iterate over sequence to fill in all counts
	for (int i=0; i<s.length(); i++) {
	    if (i<order || order<0) continue;
	    charCounts[ chars2Num(s.substring(i-order, i) + s.charAt(i)) ]++;
	    countTotal[ chars2Num(s.substring(i-order, i)) ]++;
	}
    }

    private int chars2Num(String c) {
	int res=0, i;
	for (i=0; i<c.length(); i++) {
	    int j;
	    for (j=0; j<chars.length; j++)
		if (chars[j] == c.charAt(i)) break;
	    Misc.assert(j<chars.length, "Character '"+c.charAt(i)+"' is unexpected");
	    res = (res*chars.length) + j;
	}
	//	System.out.println("char2Num("+c+")="+res);
	return res;
    }

    public double encodeLen(char a, int i) {
	if (past.length() < order || order<0) return -MyMath.log2((double)1.0/chars.length);
	int n = charCounts[ chars2Num(past.substring(i-order) + a) ];
	int d = countTotal[ chars2Num(past.substring(i-order)) ];
	return -MyMath.log2((double) n / d);
    }

    public double update(char a, int i) {
	double res = encodeLen(a, i);
	past.append(a);
	return res;
    }

    public static void main(String args[]) {
	String s = args[0];
	char[] a = {'a','t','g','c'};
	MarkovN_fitted m = new MarkovN_fitted(0, a, s);
	
	double tot = 0;
	for (int i=0; i<s.length(); i++) {
	    double r = m.update(s.charAt(i), i);
	    tot += r;
	    System.out.println(r);
	}
	System.out.println("Total entropy = "+tot+" bits/ch");
    }
}

