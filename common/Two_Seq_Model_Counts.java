package common;

public interface Two_Seq_Model_Counts extends Two_Seq_Model {
    public abstract void update_count_encA(Counts c, double w, char a, int i);
    public abstract void update_count_encB(Counts c, double w, char a, int i);
    public abstract void update_count_encBoth(Counts c, double w, char a, char b, int i, int j);

    public abstract Params counts_to_params(Counts c);
}
