
package common;

import java.io.*;

public interface Two_Seq_Model extends Serializable {
    public abstract double encA(char a, int i);
    public abstract double encB(char a, int i);
    public abstract double encBoth(char a, char b, int i, int j);

    public abstract double encode_params(double N);
}

