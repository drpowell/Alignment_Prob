
package common;

/** 
    An arbitrary order, adaptive Markov Model for an arbitrary alphabet of characters

    If order == -1 then we are a uniform model over the alphabet
*/
public class MarkovN implements Seq_Model {
    char[] chars;
    int[] charCounts;
    int[] countTotal;
    int order;
    
    StringBuffer past;

    public MarkovN(int order, char[] chars) {
	Misc.my_assert(order>=-1, "Bad order="+order);

	this.chars = chars;
	this.order = order;

	if (order>=0) {
	    charCounts = new int[ (int)Math.pow(chars.length,order+1) ];
	    countTotal = new int[ (int)Math.pow(chars.length,order) ];
	    for (int i=0; i<charCounts.length; i++) charCounts[i]=1;
	    for (int i=0; i<countTotal.length; i++) countTotal[i]=chars.length;
	}

	past = new StringBuffer();
    }

    int chars2Num(String c) {
	int res=0, i;
	for (i=0; i<c.length(); i++) {
	    int j;
	    for (j=0; j<chars.length; j++)
		if (chars[j] == c.charAt(i)) break;
	    Misc.my_assert(j<chars.length, "Character '"+c.charAt(i)+"' is unexpected");
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
	if (past.length() < order || order<0) {
		past.append(a);
		return -MyMath.log2((double)1.0/chars.length);
	}
	double res = encodeLen(a,i);
	charCounts[ chars2Num(past.substring(i-order) + a) ]++;
	countTotal[ chars2Num(past.substring(i-order)) ]++;
	past.append(a);
	return res;
    };

    public static void main(String args[]) {
	String s = args[0];
	char[] a = {'a','t','g','c'};
	MarkovN m = new MarkovN(0, a);

	double tot = 0;
	for (int i=0; i<s.length(); i++) {
	    double r = m.update(s.charAt(i), i);
	    tot += r;
	    System.out.println(r);
	}
	System.out.println("Total entropy = "+tot+" bits/ch");
    }
}

