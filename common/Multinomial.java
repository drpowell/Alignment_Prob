
package common;

public class Multinomial {
    /** Compute the MML estimate for encoding the parameters to the multinomial.
	p[] are the probabilites (should sum to 1)
	N is the number of data items that will be encoded.

	Encoding length of the Multinomial parameters are returned in bits.
    */
    static public double MMLparameter_cost(double[] p, double N) {

	double h = MyMath.factorial(p.length-1); // h(theta) - prior probabilty density = (K-1)!
	double F = 1.0/p[0];	// F will be the Fischer = N^(K-1)/(p1*p2*...*pk)
	for (int i=1; i < p.length-1; i++) {
	    F *= N / p[i];
	}

	double cost = 0.5 * MyMath.log2(1 + F/(h*h*Math.pow(12, p.length-1)));

	cost += 0.5 * (p.length-1) * MyMath.log2e;

	return cost;
    }
}
