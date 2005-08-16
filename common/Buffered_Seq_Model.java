/*
 * Created on 21/12/2004 by powell
 */
package common;

/**
 * @author powell
 *
 * Buffered_Seq_Model is a simple implementation of
 * Seq_Model that just wraps another Seq_Model.  Its purpose
 * is to cache the results of the last call to encodeLen(),
 * and hopefully speed things up!
 */
public class Buffered_Seq_Model implements Seq_Model {
    
    Seq_Model mdl;
    int last_i = 0;
    char last_a = 0;
    double last_len;

    // For some stats
    long cached = 0;
    long notCached = 0;
    
    public Buffered_Seq_Model(Seq_Model mdl) {
        this.mdl = mdl;
    }
    
    public double encodeLen(char a, int i) {
        //System.err.println("encodeLen("+a+", "+i+")");
        if (a==last_a && i==last_i) {
            cached++;
            return last_len;
        }
        
        Misc.my_assert(last_i == i, "Bad use of Buffered_Seq_Model.encodeLen. i="+i+" last_i="+last_i);
        last_a = a;
        last_len = mdl.encodeLen(a, i); 
        notCached++;
        return last_len;
    }
    
    public double update(char a, int i) {
        //System.err.println("update() : i="+i);
        Misc.my_assert(last_i == i, "Bad use of Buffered_Seq_Model.update. i="+i+" last_i="+last_i);
        last_i = i+1;
        last_a = 0;
        return mdl.update(a, i); 
    }

    public void printStats() {
        System.out.println("Buffered_Seq_Model hits cached : "+cached+"   not cached : "+notCached);
    }
}
