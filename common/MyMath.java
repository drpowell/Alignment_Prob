
package common;

public final class MyMath {
    public static final double Big_Double = Double.POSITIVE_INFINITY;

    public static final double loge2 = Math.log(2);
    public static final double log2e = 1.0/loge2;

    public static double log2(double a) {
	//return Math.log(a) / loge2;
	return Math.log(a) * log2e;
    }

    public static double exp2(double a) {
	//return Math.exp(a);
	return Math.exp(a * loge2);
    }

    public static double logplus(double a, double b) {
	//  return -log(exp(-a) + exp(-b));    But this would lose accuracy
	if (Double.isInfinite(a)) return b;
	if (Double.isInfinite(b)) return a;
	if(a>b) {double t=b;b=a;a=t;};   // make b>a
	if (b>a+50) return a;          // Approx if b >> a
	return a-log2(1+exp2(a-b));
	//return a-Math.log(1+Math.exp(a-b));
    }

    // logstar(x) - a prior over integers 1..infinity
    //              returns length in bits for encoding that integer
    //              This is for integers only!
    public static double logstar_discrete(int x) {
	if (x > 1) {
	    x = (int)log2(x);
	    return 1 + x + logstar_discrete(x);
	} else
	    return 1;
    }

    // logstar(x) - a prior over integers 1..infinity
    //              returns length in bits for encoding that integer
    //              This is Rissanen's continuous approximation to logstar_discrete
    //   according to Rissanen (1983)
    public static double logstar_continuous(double x) {
	if (x > 1) {
	    x = log2(x);
	    return x + logstar_continuous(x);
	} else
	    return log2(2.865);
    }

    public static double factorial(int N) {
	Misc.assert(N>=0, "Bad paramater to factorial");
	double res = 1;
	for (int i=2; i<=N; i++) res *= i;
	return res;
    }

    public static int min2(int a, int b) {
	return (a<b ? a : b);
    }

    public static int max2(int a, int b) {
	return (a>b ? a : b);
    }

    public static double max2(double a, double b) {
	return (a>b ? a : b);
    }

    public static double min2(double a, double b) {
	return (a<b ? a : b);
    }

    public static double min3(double a, double b, double c) {
	return (a<b ? (a<c ? a : c) : (b<c ? b : c) );
    }

    public static double min4(double a, double b, double c, double d) {
	return min2(min3(a,b,c),d);
    }

    public static void main(String args[]) {
	for (double i=1; i<4; i+=1) {
	    double v1 = MyMath.logstar_continuous(i);
	    System.out.println("logstar("+i+") = "+v1);
	}


	    /*	double s = Double.POSITIVE_INFINITY;
		for (double i=1; i<400; i+=0.01) {
		double v1 = MyMath.logstar_continuous(i);
		double v2 = MyMath.logstar_discrete((int)i);
		System.out.println(i+" "+v1+" "+v2);
		//s = MyMath.logplus(v,s);
		//System.out.println("i="+i+" v="+v+" sum="+s+" sum_p="+MyMath.exp2(-s));
	    */
    }
}

